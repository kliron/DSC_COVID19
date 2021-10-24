import os
import json
import SimpleITK as sitk
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE

DICOM_TAGS_JSON = 'H:/Dokument/Projects/dsc_covid19/dicom_tags.json' if os.name == 'nt' else './dicom_tags.json'
DICOM_TEMPORAL_POSITION_IDENTIFIER = '0020|0100'
DICOM_ACQUISITION_NUMBER = '0020|0012'


def img_histogram(img: sitk.Image, bins=None) -> None:
    """
    Get a histogram of voxel intensities for an image
    """
    plt.hist(sitk.GetArrayFromImage(img).flatten(), bins=bins)


def show_overlay(bg_img: sitk.Image, overlayed_img: sitk.Image, opacity: float = 0.01) -> None:
    sitk.Show(sitk.LabelOverlay(bg_img, overlayed_img, opacity=opacity), title='Overlay')


def get_dicom_tags(json_path: str = DICOM_TAGS_JSON, invert=False) -> dict:
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
            warnings.add(f'Warning: `Temporal position identifier` DICOM tag is not defined, falling back to acquisition number')
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


def run_n_subprocesses(commands: [[str]]):
    # We pass shell=True to execute in a shell environment and read FSL related environment variables
    procs = [Popen(c,
                   stdout=PIPE,
                   stderr=PIPE,
                   shell=True,
                   env=os.environ.copy(),
                   universal_newlines=True) for c in commands]
    for p in procs:
        # `communicate()` call will wait for process to complete.
        stdin, stderr = p.communicate()
        print(stdin)
        print(stderr)
        print(f'\nProcess returned exit code {p.returncode}')
