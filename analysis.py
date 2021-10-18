import SimpleITK as sitk
import os
import sys
import json
import numpy

# By default simpleitk will try to use the Fiji program to show images, change to slicer
os.environ['SITK_SHOW_COMMAND'] = '/Applications/Slicer.app/Contents/MacOS/Slicer'

default_dicom_tags_json_path = './dicom_tags.json'


def get_dicom_tags(json_path: str = default_dicom_tags_json_path) -> dict:
    """
    Reads a dictionary of DICOM tag to human readable tag description
    """
    with open(json_path, 'r') as r:
        js = r.read()
        return json.loads(js)


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
    tags = get_dicom_tags()
    for k, v in tags.items():
        if v == name:
            return k


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
    print(f'Read series with size {series.GetSize()}')

    tags = get_dicom_tags()
    metadata = {}
    for k in reader.GetMetaDataKeys(slice=0):
        try:
            desc = tags[k]
        except KeyError:
            print(f'DICOM tag {k} is missing a descriptive name')
            desc = k
        metadata[desc] = reader.GetMetaData(0, k)

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
        # Acquisition Number starts at 1 and is strictly increasing by 1
        acquisition_number = int(file_reader.GetMetaData('0020|0012'))
        if len(file_lists) < acquisition_number:
            file_lists.append([dcm_name])
        else:
            file_lists[acquisition_number-1].append(dcm_name)

    tags = get_dicom_tags()

    series_list = []
    for dcm_names_list in file_lists:
        reader.SetFileNames(dcm_names_list)
        metadata = {}
        reader.LoadPrivateTagsOn()
        reader.MetaDataDictionaryArrayUpdateOn()
        series = reader.Execute()
        for k in reader.GetMetaDataKeys(slice=0):
            try:
                desc = tags[k]
            except KeyError:
                print(f'DICOM tag {k} is missing a descriptive name')
                desc = k
            metadata[desc] = reader.GetMetaData(0, k)

        series_list.append((series, metadata))

    return series_list


def register(ser1: sitk.Image, ser2: sitk.Image) -> sitk.Image:
    """

    :param ser1:
    :param ser2:
    :return:
    """


if __name__ == '__main__':
    # dicom_dir = '/Users/kliron/Downloads/tmp/data/DICOM/0000931B/AAC70388/AA080C75'
    dicom_dir = sys.argv[1]
    # ser, meta = read_3d_series(os.path.join(dicom_dir, '000055DE'))
    data = read_4d_series(os.path.join(dicom_dir, '000055DE'))
    sitk.Show(data[7][0], title="First frame of series", debugOn=False)

    print((f'Sequence: {meta["Sequence Name"]}\n'
           f'FA: {meta["Flip Angle"]}\n'
           f'TE: {meta["Echo Time"]}\n'
           f'TR: {meta["Repetition Time"]}\n'
           f'MFS: {meta["Magnetic Field Strength"]}\n'
           f'Matrix: {meta["Acquisition Matrix"]}'))

    get_tag_by_name('Content Time')
    # Convert to numpy
    # arr = sitk.GetArrayFromImage(ser)

