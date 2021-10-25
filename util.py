import os
import json
import SimpleITK as sitk
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from typing import Dict


DICOM_TAGS_JSON = 'H:/Dokument/Projects/dsc_covid19/dicom_tags.json' if os.name == 'nt' else './dicom_tags.json'
TEMPORAL_POSITION_IDENTIFIER_TAG = '0020|0100'
ACQUISITION_NUMBER_TAG = '0020|0012'


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


def get_dicom_tag_by_name(name: str) -> str:
    """
    :param name: Human readable description of a DICOM tag
    :return: (Machine readable) DICOM tag
    """
    for key, val in get_dicom_tags().items():
        if val == name:
            return key


def read_3d_series(dcm_path: str) -> (sitk.Image, Dict[str, str]):
    """
    Reads a DICOM series in a SimpleITK.Image
    :param dcm_path: Path to directory containing DICOM files
    :return: SimpleITK.Image and Dict of DICOM tags read from 1st slice
    """
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dcm_path)
    reader.SetFileNames(dicom_names)
    reader.LoadPrivateTagsOn()
    reader.MetaDataDictionaryArrayUpdateOn()
    series = reader.Execute()
    tags = get_dicom_tags()
    dicom_metadata = {}
    for key in reader.GetMetaDataKeys(slice=0):
        try:
            desc = tags[key]
        except KeyError:
            print(f'DICOM tag {key} is missing a descriptive name')
            desc = key
        dicom_metadata[desc] = reader.GetMetaData(0, key)
    return series, dicom_metadata


def read_4d_series(dcm_path: str) -> ([sitk.Image], Dict[str, str]):
    """
    There is no built-in way to read 4D series like a DSC perfusion time series in ITK. We do what is described
    here: https://github.com/SimpleITK/SimpleITK/issues/879 to get an ordered list of 3D volumes instead, each being one
     frame in the time series.
    :param dcm_path: Path to the directory containing all DICOM files
    :return: A list of sitk.Image and dictionary of DICOM tags read from 1st slice of 1st frame
    """
    reader = sitk.ImageSeriesReader()
    reader.LoadPrivateTagsOn()
    reader.MetaDataDictionaryArrayUpdateOn()
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
            position_identifier = int(file_reader.GetMetaData(TEMPORAL_POSITION_IDENTIFIER_TAG))
        except RuntimeError as e:
            warnings.add(f'Warning: `Temporal position identifier` DICOM tag is not defined, falling back to acquisition number')
            position_identifier = int(file_reader.GetMetaData(ACQUISITION_NUMBER_TAG))
        if len(file_lists) < position_identifier:
            file_lists.append([dcm_name])
        else:
            file_lists[position_identifier-1].append(dcm_name)

    series_list = []
    dicom_metadata = {}
    tags = get_dicom_tags()
    for idx, dcm_names_list in enumerate(file_lists):
        reader.SetFileNames(dcm_names_list)
        series_list.append(reader.Execute())
        if idx == 0:
            for key in reader.GetMetaDataKeys(slice=0):
                try:
                    desc = tags[key]
                except KeyError:
                    print(f'DICOM tag {key} is missing a descriptive name')
                    desc = key
                dicom_metadata[desc] = reader.GetMetaData(0, key)

    if len(warnings):
        for w in warnings:
            print(w)

    return series_list, dicom_metadata


def run_n_subprocesses(commands: [str]) -> [int]:
    """
    Runs N commands in N subprocesses. Because we need the programs to be able to read environment variables, we pass
    shell=True. When shell is True, the commands need to be raw strings, not lists of arguments.
    :param commands: The raw command string to execute in a shell.
    :return: A list of exit codes from each subprocess
    """
    # We pass shell=True to execute in a shell environment and read FSL related environment variables
    procs = [Popen(c,
                   stdout=PIPE,
                   stderr=PIPE,
                   shell=True,
                   env=os.environ.copy(),
                   universal_newlines=True) for c in commands]
    for p in procs:
        # `communicate()` call will wait for process to complete.
        stdout, stderr = p.communicate()
        if p.returncode:
            print(stdout)
            print(stderr)

    return [p.returncode for p in procs]
