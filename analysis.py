#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os.path
import sys
import argparse
from typing import Dict, Any
import numpy as np

# This is to make relative imports work both when this file is run "on the fly" in an IDE
# and when it is run as a script (as __main__). For the first case to work, the console
# must of course cd in the directory this file is in.
parent_module = sys.modules['.'.join(__name__.split('.')[:-1]) or '__main__']
if __name__ == '__main__' or parent_module.__name__ == '__main__':
    from paths import paths
    from util import *
else:
    from .paths import paths
    from .util import *

os.environ[
    'SITK_SHOW_COMMAND'] = 'H:/Slicer 4.11.20210226/Slicer.exe' if os.name == 'nt' else '/Applications/Slicer.app/Contents/MacOS/Slicer'
os.environ['FSLOUTPUTTYPE'] = 'NIFTI'  # this is to make bet2 tool not bitch when we execute it with subprocess

DATA_ROOT = 'H:/Dokument' if os.name == 'nt' else '/Users/kliron'
DICOM_ROOT = 'Data/dsc_covid19/Examinations'
DERIVED_ROOT = 'Data/dsc_covid19/Derived'
SEGMENTATION_FILENAME = 'Plexus.nrrd'
# BET (brain extraction tool) does not support 4D images. We need to convert each 4D series to a list of 3D images,
# save it to a temporary nifti file, call BET to process it, and convert the output back to nrrd which is the format the
# DSC Slicer module expects
PERFUSION_3D_TMP_NIFTI_FILE = '_Perfusion3D.nii'
PERFUSION_SKULL_STRIPPED_3D_FILE = '_stripped_Perfusion3d.nii'
BET_EXECUTABLE = '/usr/local/fsl/bin/bet2'
# This file is created by importing the original DICOM series as a MultiVolume in Slicer and saving as a nrrd file
# We need this to copy the metadata tags back to the final (skull-stripped) MultiVolume that the DSC module is going to
# work with
ORIGINAL_PERFUSION_MULTIVOLUME = 'Perfusion.nrrd'
# This is the file DSCMRIAnalysis module will work on
FINAL_PERFUSION_MULTIVOLUME = 'final_Perfusion.nrrd'
DSC_EXECUTABLE = '/Applications/Slicer.app/Contents/Extensions-30329/DSCMRIAnalysis/lib/Slicer-4.13/cli-modules/DSCMRIAnalysis'


def edit_segmentation(uid: str,
                      params: Dict[str, Any],
                      foreground_label: int = 1,
                      background_label: int = 0,
                      show: bool = False) -> [sitk.Image]:
    """
    :param uid: Unique identification number for an exam.
    :param params: Path and slice number parameters defined in `paths` module.
    :param foreground_label: Foreground integer value of the segmentation LabelMap (Default 1).
    :param background_label: Background integer value of the segmentation LabelMap (Default 0).
    :param show: If true, show an overlay of the final segmentation and the first slice of the first frame of the series.
    :return: The 3D image list to process in other functions downstream
    """
    print(f'Editing segmentation for {uid}')
    derived_root = os.path.join(DATA_ROOT, DERIVED_ROOT, uid)

    dcm_path = os.path.join(DATA_ROOT, DICOM_ROOT, uid, 'DICOM', params['dicom_dir'])
    series, dicom_tags = read_4d_series(dcm_path)
    segm = sitk.ReadImage(os.path.join(derived_root, SEGMENTATION_FILENAME))

    # Get the average voxel values for the non-contrast series
    noncontrast_imgs = [series[frame] for frame in params['noncontrast_frames']]
    contrast_imgs = [series[frame] for frame in params['contrast_frames']]
    noncontrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in noncontrast_imgs]) / len(noncontrast_imgs)).astype(int)
    contrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in contrast_imgs]) / len(contrast_imgs)).astype(int)

    # To edit the segmentation itself, we save the indexes of each one of its voxels
    segm_arr = sitk.GetArrayFromImage(segm)
    mask_indexes = list(zip(*np.where(segm_arr == foreground_label)))

    # Segmentations are saved as labelmaps where each voxel belonging to the segmentation has the value 1 (default)
    # and the rest the value 0. Thus, by using orig_img_arr[segm == 1] we can index into the original image and get the
    # voxels corresponding to the volume of the segmentation.
    #
    # Remove voxels whose noncontrast minus contrast average intensity is lower than one standard deviation of the
    # average noncontrast intensity.
    remove_condition = (noncontrast_avg_array[segm_arr == 1] - contrast_avg_array[segm_arr == 1]) < np.std(
        noncontrast_avg_array[segm_arr == 1])

    # Setting the voxel to 0 (the default label for background) removes it from the segmentation
    for idx, cond in zip(mask_indexes, remove_condition):
        if cond:
            segm_arr[idx] = background_label

    # The final segmentation contains voxels that take up contrast and don't contain calcifications
    final_segm = sitk.GetImageFromArray(segm_arr)
    final_segm.CopyInformation(segm)
    show_overlay(noncontrast_imgs[0], final_segm) if show else ()
    sitk.WriteImage(final_segm, os.path.join(derived_root, f'final_{SEGMENTATION_FILENAME}'))
    # with open(os.path.join(derived_root, 'dicom_tags.json'), 'w') as w:
    #    w.write(json.dumps(dicom_tags, indent=4))

    return series


def extract_brain(uid: str, series: [sitk.Image], nii_dir: str = 'nii') -> None:
    """
    For each of N 3D volumes in `series` a NIFTI file will be saved under DATA_ROOT/DERIVED_ROOT/${uid}/tmpdir/ directory.
    BET will be run on N parallel subprocesses to extract the brain and save a new NIFTI file. Each of the new files
    will be joined in a 4D image with JoinSeriesImageFilter which will be saved as a .nrrd file.
    :param uid: Unique exam id
    :param series: A list of 3D volume images representing a 4D DSC series
    :param nii_dir: A subdirectory under which we save all temporary NIFTI files
    :return: Paths to the final .nrrd file that was created
    """
    nii_path = os.path.join(DATA_ROOT, DERIVED_ROOT, uid, nii_dir)
    if not os.path.exists(nii_path):
        os.mkdir(nii_path)
    else:
        if not os.path.isdir(nii_path) or not os.access(nii_path, os.W_OK):
            raise Exception(f'ERROR {nii_path} is not a writeable directory')

    in_files = []
    out_files = []
    # Write series as a bunch of nifti files
    for idx, img in enumerate(series):
        file = str(idx) + PERFUSION_3D_TMP_NIFTI_FILE
        sitk.WriteImage(img, os.path.join(nii_path, file))
        in_files.append(os.path.join(nii_path, file))
        out_files.append(os.path.join(nii_path, f'{idx}_stripped_{file}'))

    # Call BET on N subprocesses to extract brain
    commands = [' '.join([BET_EXECUTABLE, i, o]) for i, o in zip(in_files, out_files)]
    exit_codes = run_n_subprocesses(commands)
    if any(exit_codes):
        print('Warning: some of the subprocesses returned nonzero exit status')

    # Delete temporary nifti files
    for f in in_files:
        os.unlink(f)


def copy_dicom_tags_to_final_multivolume(uid: str) -> None:
    """
    Copies the image metadata needed for DSC module to be able to understand the MultiVolume from the original perfusion
    file
    :param uid:
    :return:
    """
    orig_path = os.path.join(DATA_ROOT, DERIVED_ROOT, uid, ORIGINAL_PERFUSION_MULTIVOLUME)
    final_path = os.path.join(DATA_ROOT, DERIVED_ROOT, uid, FINAL_PERFUSION_MULTIVOLUME)
    orig = sitk.ReadImage(orig_path)
    meta = {k: orig.GetMetaData(k) for k in orig.GetMetaDataKeys()}
    final = sitk.ReadImage(final_path)
    for k, v in meta.items():
        final.SetMetaData(k, v)
    # Overwrites final_path
    sitk.WriteImage(final, final_path)


def run_perfusion(uid: str) -> None:
    """
    Uses DSCMRIAnalysis Slicer CLI module (https://www.slicer.org/w/index.php/Documentation/Nightly/Modules/DSC_MRI_Analysis)
    to extract perfusion maps for all examinations using the final segmentations created by the step above as ROI.

    This is an example of the command that will run in a subprocess:
    /Applications/Slicer.app/Contents/Extensions-30329/DSCMRIAnalysis/lib/Slicer-4.13/cli-modules/DSCMRIAnalysis \
            --roiMask /Users/kliron/Data/dsc_covid19/Segmentations/SEKBF00012197064/final_Plexus.nrrd \
            --aifMask /Users/kliron/Data/dsc_covid19/Segmentations/SEKBF00012197064/Artery.nrrd  \
            --outputK1 path_to_K1.nrrd \
            --outputK2 path_to_K2.nrrd \
            --outputMTT path_to_MTT.nrrd \
            --outputCBF path_to_CBF.nrrd \
            /Users/kliron/Downloads/multivolume.nrrd

    :param uid: Unique exam ID string
    :return: None
    """
    print(f'Computing perfusion maps for {uid}')
    root_path = os.path.join(DATA_ROOT, DERIVED_ROOT, uid)
    args = ['--roiMask', os.path.join(root_path, 'final_Plexus.nrrd'),
            '--aifMask', os.path.join(root_path, 'Artery.nrrd'),
            '--outputK1', os.path.join(root_path, 'K1.nrrd'),
            '--outputK2', os.path.join(root_path, 'K2.nrrd'),
            '--outputMTT', os.path.join(root_path, 'MTT.nrrd'),
            '--outputCBF', os.path.join(root_path, 'CBF.nrrd'),
            os.path.join(root_path, FINAL_PERFUSION_MULTIVOLUME)]

    print('Executing DSC CLI module...')
    print(DSC_EXECUTABLE + ' \\ \n' +
          ' \\ \n    '.join([f'{args[i]} {args[j]}' for i, j in zip(range(0, len(args), 2), range(1, len(args), 2))]) +
          ' \\ \n    ' + args[len(args) - 1])

    exit_codes = run_n_subprocesses([DSC_EXECUTABLE + ' ' + ' '.join(args)])
    if any(exit_codes):
        print('Warning: some of the subprocesses returned nonzero exit status')


def main():
    parser = argparse.ArgumentParser(description='DSC perfusion analysis.',
                                     usage=f'{sys.argv[0]} [--uid UID] --segmentation | --perfusion\n')
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
            run_perfusion(uid)
        else:
            edit_segmentation(uid, params)
            run_perfusion(uid)
    else:
        if args.segmentation:
            for uid, params in exams:
                edit_segmentation(uid, params)

        if args.perfusion:
            for uid, params in exams:
                run_perfusion(uid)


if __name__ == '__main__':
    main()
