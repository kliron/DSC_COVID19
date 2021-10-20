import SimpleITK as sitk
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# By default simpleitk will try to use the Fiji program to show images, change to slicer
os.environ['SITK_SHOW_COMMAND'] = 'H:/Slicer 4.11.20210226/Slicer.exe' if os.name == 'nt' else \
    '/Applications/Slicer.app/Contents/MacOS/Slicer'
default_dicom_tags_json_path = 'H:/Dokument/Projects/DSC_COVID19/dicom_tags.json' if os.name == 'nt' else \
    './dicom_tags.json'
ACQUISITION_NUMBER_DICOM_TAG = '0020|0012'


def get_dicom_tags(json_path: str = default_dicom_tags_json_path, invert=False) -> dict:
    """
    Reads a dictionary of DICOM tag to human readable tag description.
    if `invert` is True, switch key-value pairs
    """
    with open(json_path, 'r') as r:
        tags = json.loads(r.read())
        return tags if not invert else {val: key for key, val in tags.items()}


# Flip Angle -> '0018|1314'
# Echo Time -> '0018|0081'
# Repetition Time -> '0018|0080'
# Inversion Time -> '0018|0082'
# Sequence Name -> '0018|0024'
# Magnetic Field Strength -> '0018|0087'
# Echo Train Length -> '0018|0091'
# In-plane Phase Encoding Direction -> '0018|1312'
# Acquisition Matrix -> '0018|1310'
# Transmit Coil Name -> '0018|1251'
# Receive Coil Name -> '0018|1250'
# Content Time -> '0008|0033'
def get_tag_by_name(name: str) -> str:
    """
    :param name: Tag human readable description
    :return: DICOM tag
    """
    for key, val in get_dicom_tags().items():
        if val == name:
            return key


def read_3d_series(dsc_path: str) -> (sitk.Image, dict):
    """
    Reads a DICOM series and returns a tuple of the (multislice) image and a dictionary of metadata from the first slice
    :param dsc_path: Path to directory containing DICOM files
    :return: SimpleITK.Image
    """
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dsc_path)
    reader.SetFileNames(dicom_names)
    reader.LoadPrivateTagsOn()
    reader.MetaDataDictionaryArrayUpdateOn()
    series = reader.Execute()
    tags = get_dicom_tags()
    metadata = {}
    for key in reader.GetMetaDataKeys(slice=0):
        try:
            desc = tags[key]
        except KeyError:
            print(f'DICOM tag {key} is missing a descriptive name')
            desc = key
        metadata[desc] = reader.GetMetaData(0, key)
    return series, metadata


def read_4d_series(dsc_path: str) -> [(sitk.Image, dict)]:
    """
    There is no built-in way to read 4D series like a DSC perfusion time series in ITK. We do what is described
    here: https://github.com/SimpleITK/SimpleITK/issues/879 to get an ordered list of 3D volumes instead, each being one
     frame in the time series.
     Returns a list of tuples, each containing the 3D Volume of one frame and a dictionary of DICOM tags read from its
     first slice.
    :param dsc_path: Path to the directory containing all DICOM files
    :return: [(sitk.Image, dict)]
    """
    reader = sitk.ImageSeriesReader()
    all_dicom_names = reader.GetGDCMSeriesFileNames(dsc_path)
    file_reader = sitk.ImageFileReader()
    file_lists = []
    for dcm_name in all_dicom_names:
        file_reader.SetFileName(dcm_name)
        file_reader.ReadImageInformation()
        _ = file_reader.Execute()
        # Acquisition Number starts at 1 and strictly increases by 1
        acquisition_number = int(file_reader.GetMetaData(ACQUISITION_NUMBER_DICOM_TAG))
        if len(file_lists) < acquisition_number:
            file_lists.append([dcm_name])
        else:
            file_lists[acquisition_number-1].append(dcm_name)

    series_list = []
    for dcm_names_list in file_lists:
        reader.SetFileNames(dcm_names_list)
        reader.LoadPrivateTagsOn()
        reader.MetaDataDictionaryArrayUpdateOn()
        series = reader.Execute()
        metadata = {k: reader.GetMetaData(0, k) for k in reader.GetMetaDataKeys(slice=0)}
        series_list.append((series, metadata))

    return series_list


def main():
    """
    Segmentations are saved as labelmaps where each voxel belonging to the segmentation has the value 1 (default)
    and the rest the value 0. Thus, by using orig_img_arr[segm==1] we can index into the original image and get the
    voxels corresponding to the volume of the segmentation.

    :return: None
    """
    data_root = 'H:/Dokument' if os.name == 'nt' else '/Users/kliron'
    dicom_root = 'Data/dsc_covid19/Examinations'
    segmentations_dir = 'Data/dsc_covid19/Segmentations'
    segmentation_file = 'Plexus.nrrd'
    segmentation_label = 1  # default
    background_label = 0    # default
    with open('./paths.json') as r:
        js = r.read()
        dicom_dirs = json.loads(js)
    [(uid, params)] = dicom_dirs.items()
    dicom_path = os.path.join(data_root, dicom_root, uid, 'DICOM', params['dicom_dir'])
    data = read_4d_series(dicom_path)
    segm = sitk.ReadImage(os.path.join(data_root, segmentations_dir, uid, segmentation_file))
    # Get the average voxel values for the non-contrast series
    noncontrast_imgs = [data[frame][0] for frame in params['noncontrast_frames']]
    contrast_imgs = [data[frame][0] for frame in params['contrast_frames']]

    # To see the segmentation overlayed on the actual frame:
    # sitk.Show(sitk.LabelOverlay(contrast_img, segm, opacity=0.01), title='Overlay')

    #
    # From the  Bouzerar et al paper (Neuroradiology 2013, DOI 10.1007/s00234-013-1290-2):
    # "For each pixel, the pre-bolus baseline S0 was estimated by
    #  averaging 10 to 12 points before tB. For an individual voxel,
    #  significant enhancement occurred when the mean signal intensity
    #  over the final ten time points was greater than 1 standard
    #  deviation over the voxel average baseline S0. After comparison
    #  between the post-bolus signal and the pre-bolus baseline S0,
    #  pixels within the mask that did not demonstrate significant signal
    #  enhancement were categorized as non-enhancing pixels and
    #  included for the calculation of ΔR*(t)."
    #
    noncontrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in noncontrast_imgs])/len(noncontrast_imgs)).astype(int)
    contrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in contrast_imgs])/len(contrast_imgs)).astype(int)
    non_contrast_std = np.round(np.std(noncontrast_avg_array)).astype(int)

    # To edit the segmentation itself, we save the indexes of each one of its voxels:
    segm_arr = sitk.GetArrayFromImage(segm)
    mask_indexes = list(zip(*np.where(segm_arr == segmentation_label)))

    remove_condition = contrast_avg_array[segm_arr == 1] > np.std(noncontrast_avg_array[segm_arr == 1])

    for idx, cond in zip(mask_indexes, remove_condition):
        if cond:
            # Setting the voxel to 0 (the default label for background) removes it from the segmentation
            segm_arr[idx] = background_label

    final_segm = sitk.GetImageFromArray(segm_arr)
    final_segm.CopyInformation(segm)
    sitk.WriteImage(final_segm, os.path.join(data_root, segmentations_dir, uid, 'Plexus_final.nrrd'))


if __name__ == '__main__':
    main()

