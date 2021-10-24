#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from typing import Dict, Any

import SimpleITK
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


os.environ['SITK_SHOW_COMMAND'] = 'H:/Slicer 4.11.20210226/Slicer.exe' if os.name == 'nt' else '/Applications/Slicer.app/Contents/MacOS/Slicer'
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
# Each 4D perfusion DICOM series is written to a single nrrd file to be readable from DSC Slicer module
PERFUSION_4D_FILE = 'Perfusion4D.nrrd'
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
    dcm_path = os.path.join(DATA_ROOT, DICOM_ROOT, uid, 'DICOM', params['dicom_dir'])
    data = read_4d_series(dcm_path)
    segm = sitk.ReadImage(os.path.join(DATA_ROOT, DERIVED_ROOT, uid, SEGMENTATION_FILENAME))

    # Get the average voxel values for the non-contrast series
    noncontrast_imgs = [data[frame][0] for frame in params['noncontrast_frames']]
    contrast_imgs = [data[frame][0] for frame in params['contrast_frames']]
    noncontrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in noncontrast_imgs])/len(noncontrast_imgs)).astype(int)
    contrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in contrast_imgs])/len(contrast_imgs)).astype(int)

    # To edit the segmentation itself, we save the indexes of each one of its voxels
    segm_arr = sitk.GetArrayFromImage(segm)
    mask_indexes = list(zip(*np.where(segm_arr == foreground_label)))

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
    show_overlay(noncontrast_imgs[0], final_segm) if show else ()
    sitk.WriteImage(final_segm, os.path.join(DATA_ROOT, DERIVED_ROOT, uid, f'final_{SEGMENTATION_FILENAME}'))

    return [it[0] for it in data]


def extract_brain(uid: str, series: [sitk.Image]) -> None:
    # Save as nifti
    in_files = []
    out_files = []
    for idx, img in enumerate(series):
        file = str(idx) + PERFUSION_3D_TMP_NIFTI_FILE
        sitk.WriteImage(img, os.path.join(DATA_ROOT, DERIVED_ROOT, uid, file))
        in_files.append(os.path.join(DATA_ROOT, DERIVED_ROOT, uid, file))
        out_files.append(os.path.join(DATA_ROOT, DERIVED_ROOT, uid, f'{idx}_stripped_{file}'))

    commands = [' '.join([BET_EXECUTABLE, i, o]) for i, o in zip(in_files, out_files)]
    run_n_subprocesses(commands)
    # Read back nifti and save 4D .nrrd image instead
    # TODO
    reader = sitk.ImageFileReader()
    reader.SetFileName()


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
            os.path.join(root_path, PERFUSION_4D_FILE)]

    print('Executing DSC CLI module...')
    print(DSC_EXECUTABLE + ' \\ \n' +
          ' \\ \n    '.join([f'{args[i]} {args[j]}' for i, j in zip(range(0, len(args), 2), range(1, len(args), 2))]) +
          ' \\ \n    ' + args[len(args)-1])

    run_n_subprocesses([DSC_EXECUTABLE + ' ' + ' '.join(args)])


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
