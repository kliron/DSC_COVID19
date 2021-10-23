# -*- coding: utf-8 -*-
import SimpleITK as sitk
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import argparse
from typing import Any

# This is to make relative imports work both when this file is run "on the fly" in an IDE
# and when it is run as a script (as __main__). For the first case to work, the console
# must of course cd in the directory this file is in.
parent_module = sys.modules['.'.join(__name__.split('.')[:-1]) or '__main__']
if __name__ == '__main__' or parent_module.__name__ == '__main__':
    from paths import paths
else:
    from . import paths


os.environ['SITK_SHOW_COMMAND'] = 'H:/Slicer 4.11.20210226/Slicer.exe' if os.name == 'nt' else \
    '/Applications/Slicer.app/Contents/MacOS/Slicer'
default_dicom_tags_json_path = 'H:/Dokument/Projects/dsc_covid19/dicom_tags.json' if os.name == 'nt' else \
    './dicom_tags.json'

DATA_ROOT = 'H:/Dokument' if os.name == 'nt' else '/Users/kliron'
DICOM_ROOT = 'Data/dsc_covid19/Examinations'
SEGMENTATIONS_ROOT = 'Data/dsc_covid19/Segmentations'
SEGMENTATION_FILENAME = 'Plexus.nrrd'
DICOM_TEMPORAL_POSITION_IDENTIFIER = '0020|0100'
DICOM_ACQUISITION_NUMBER = '0020|0012'


def get_dicom_tags(json_path: str = default_dicom_tags_json_path, invert=False) -> dict:
    """
    Reads a dictionary of DICOM tag to human readable tag description.
    if `invert` is True, the key-value pairs are switched
    """
    with open(json_path, 'r') as r:
        tags = json.loads(r.read())
        return tags if not invert else {val: key for key, val in tags.items()}


def get_tag_by_name(name: str) -> str:
    """
    :param name: Human readable description of a DICOM tag
    :return: (Machine readable) DICOM tag
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


def read_4d_series(dcm_path: str) -> [(sitk.Image, dict)]:
    """
    There is no built-in way to read 4D series like a DSC perfusion time series in ITK. We do what is described
    here: https://github.com/SimpleITK/SimpleITK/issues/879 to get an ordered list of 3D volumes instead, each being one
     frame in the time series.
     Returns a list of tuples, each containing the 3D Volume of one frame and a dictionary of DICOM tags read from its
     first slice.
    :param dcm_path: Path to the directory containing all DICOM files
    :return: [(sitk.Image, dict)]
    """
    reader = sitk.ImageSeriesReader()
    all_dicom_names = reader.GetGDCMSeriesFileNames(dcm_path)
    file_reader = sitk.ImageFileReader()
    file_lists = []
    warnings = set([])
    for dcm_name in all_dicom_names:
        file_reader.SetFileName(dcm_name)
        file_reader.ReadImageInformation()
        _ = file_reader.Execute()
        # Prefer temporal position identifier if it is defined, fall back to acquisition number if it is not.
        try:
            position_identifier = int(file_reader.GetMetaData(DICOM_TEMPORAL_POSITION_IDENTIFIER))
        except RuntimeError as e:
            warnings.add(f'Warning: Temporal position identifier key is not defined, falling back to acquisition number')
            position_identifier = int(file_reader.GetMetaData(DICOM_ACQUISITION_NUMBER))
        if len(file_lists) < position_identifier:
            file_lists.append([dcm_name])
        else:
            file_lists[position_identifier-1].append(dcm_name)

    series_list = []
    for dcm_names_list in file_lists:
        reader.SetFileNames(dcm_names_list)
        reader.LoadPrivateTagsOn()
        reader.MetaDataDictionaryArrayUpdateOn()
        series = reader.Execute()
        metadata = {k: reader.GetMetaData(0, k) for k in reader.GetMetaDataKeys(slice=0)}
        series_list.append((series, metadata))

    if len(warnings):
        for w in warnings:
            print(w)

    return series_list


def img_histogram(img: sitk.Image, bins=None) -> None:
    """
    Get a histogram of voxel intensities for an image
    """
    plt.hist(sitk.GetArrayFromImage(img).flatten(), bins=bins)


def show_overlay(bg_img: sitk.Image, overlayed_img: sitk.Image, opacity: float = 0.01) -> None:
    sitk.Show(sitk.LabelOverlay(bg_img, overlayed_img, opacity=opacity), title='Overlay')


def edit_segmentation(uid: str,
                      params: dict[str, Any],
                      foreground_label: int = 1,
                      background_label: int = 0,
                      show: bool = False) -> None:
    """
    :param uid: Unique identification number for an exam.
    :param params: Path and slice number parameters defined in `paths` module.
    :param foreground_label: Foreground integer value of the segmentation LabelMap (Default 1).
    :param background_label: Background integer value of the segmentation LabelMap (Default 0).
    :param show: If true, show an overlay of the final segmentation and the first slice of the first frame of the series.
    :return: None
    """
    print(f'Edit segmentation for {uid}')
    dcm_path = os.path.join(DATA_ROOT, DICOM_ROOT, uid, 'DICOM', params['dicom_dir'])
    data = read_4d_series(dcm_path)
    segm = sitk.ReadImage(os.path.join(DATA_ROOT, SEGMENTATIONS_ROOT, uid, SEGMENTATION_FILENAME))

    # Get the average voxel values for the non-contrast series
    noncontrast_imgs = [data[frame][0] for frame in params['noncontrast_frames']]
    contrast_imgs = [data[frame][0] for frame in params['contrast_frames']]
    noncontrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in noncontrast_imgs])/len(noncontrast_imgs)).astype(int)
    contrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in contrast_imgs])/len(contrast_imgs)).astype(int)

    # To edit the segmentation itself, we save the indexes of each one of its voxels
    segm_arr = sitk.GetArrayFromImage(segm)
    mask_indexes = list(zip(*np.where(segm_arr == foreground_label)))

    #
    # Segmentations are saved as labelmaps where each voxel belonging to the segmentation has the value 1 (default)
    # and the rest the value 0. Thus, by using orig_img_arr[segm == 1] we can index into the original image and get the
    # voxels corresponding to the volume of the segmentation.
    #
    # Remove voxels whose noncontrast minus contrast average intensity is lower than one standard deviation of the
    # average noncontrast intensity.
    remove_condition = (noncontrast_avg_array[segm_arr == 1] - contrast_avg_array[segm_arr == 1]) < np.std(noncontrast_avg_array[segm_arr == 1])

    # Setting the voxel to 0 (the default label for background) removes it from the segmentation
    for idx, cond in zip(mask_indexes, remove_condition):
        if cond:
            segm_arr[idx] = background_label

    # The final segmentation contains voxels that take up contrast and don't contain calcifications
    final_segm = sitk.GetImageFromArray(segm_arr)
    final_segm.CopyInformation(segm)

    # Erosion does not work for these images as they are just 1-2 voxels wide and the erosion converts it all to bg.
    # Leaving it here for reference.
    # f = sitk.BinaryErodeImageFilter()
    # f.SetKernelType(sitk.sitkBall)
    # f.SetKernelRadius(1)
    # f.SetBackgroundValue(0)
    # f.SetForegroundValue(1)
    # eroded_segm: sitk.Image = f.Execute(final_segm)

    if show:
        show_overlay(noncontrast_imgs[0], final_segm)

    sitk.WriteImage(final_segm, os.path.join(DATA_ROOT, SEGMENTATIONS_ROOT, uid, f'final_{SEGMENTATION_FILENAME}'))


def run_perfusion(uid: str, params: dict[str, Any], show: bool = False) -> None:
    """
    Uses DSCMRIAnalysis slicer CLI module to extract perfusion maps for all examinations using the final segmentations
    created by the step above as ROI.
    :param uid: Unique exam ID string
    :param params: Path and slice number parameters defined in `paths` module.
    :param show: If true, shows the resulting perfusion maps.
    :return: None
    """


def main():
    parser = argparse.ArgumentParser(description='DSC perfusion analysis.',
                                     usage=f'{sys.argv[0]} [--uid UID] [--segmentation | --perfusion]\n'
                                           f'If neither --segmentation nor --perfusion argument is provided, '
                                           f'first segmentation and then perfusion pipeline will execute')
    parser.add_argument('-u', '--uid', help='Specify a single UID to process', nargs=1)
    parser.add_argument('-s', '--segmentation', action='store_true', help='Process segmentations')
    parser.add_argument('-p', '--perfusion', action='store_true', help='Run perfusion analysis')
    args = parser.parse_args()

    exams = paths.items()
    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)

    if args.uid:
        uid = args.uid[0]
        params = paths[uid]
        if args.segmentation:
            edit_segmentation(uid, params, show=True)
        elif args.perfusion:
            run_perfusion(uid, params)
        else:
            edit_segmentation(uid, params, show=True)
            run_perfusion(uid, params)
    else:
        if args.segmentation:
            for uid, params in exams:
                edit_segmentation(uid, params, show=True)
        elif args.perfusion:
            for uid, params in exams:
                run_perfusion(uid, params, show=True)
        else:
            for uid, params in exams:
                edit_segmentation(uid, params, show=True)
                run_perfusion(uid, params, show=True)


if __name__ == '__main__':
    main()
