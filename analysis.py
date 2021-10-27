#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os.path
import sys
import argparse
import numpy as np
import pandas as pd
from statsmodels.stats.weightstats import ttest_ind
import shutil

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
DICOM_ROOT = os.path.join(DATA_ROOT, 'Data/dsc_covid19/Examinations')
DERIVED_ROOT = os.path.join(DATA_ROOT, 'Data/dsc_covid19/Derived')
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
NII_DIR = 'nii'
# This is the file produced by reading the skull stripped nii files
# This is the file DSCMRIAnalysis module will work on
FINAL_PERFUSION_MULTIVOLUME = 'final_Perfusion.nrrd'
DSC_EXECUTABLE = '/Applications/Slicer.app/Contents/Extensions-30329/DSCMRIAnalysis/lib/Slicer-4.13/cli-modules/DSCMRIAnalysis'
PLEXUS_STATS_FILE = os.path.join(DERIVED_ROOT, 'plexus_stats.xlsx')
PERFUSION_DATA_FILE = os.path.join(DERIVED_ROOT, 'perfusion_data.xlsx')


# Step 1
def edit_segmentation(uid: str,
                      foreground_label: int = 1,
                      background_label: int = 0,
                      show: bool = False) -> [sitk.Image]:
    """
    :param uid: Unique identification number for an exam.
    :param foreground_label: Foreground integer value of the segmentation LabelMap (Default 1).
    :param background_label: Background integer value of the segmentation LabelMap (Default 0).
    :param show: If true, show an overlay of the final segmentation and the first slice of the first frame of the series.
    :return: The 3D image list to process in other functions downstream
    """
    print(f'Editing segmentation for {uid}')
    derived_root = os.path.join(DERIVED_ROOT, uid)
    params = paths[uid]
    dcm_path = os.path.join(DICOM_ROOT, uid, 'DICOM', params['dicom_dir'])
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


# Step 2
def extract_brain(uid: str, series: [sitk.Image]) -> None:
    """
    For each of N 3D volumes in `series` a NIFTI file will be saved under DERIVED_ROOT/${uid}/tmpdir/ directory.
    BET will be run on N parallel subprocesses to extract the brain and save a new NIFTI file. Each of the new files
    will be joined in a 4D image with JoinSeriesImageFilter which will be saved as a .nrrd file.
    :param uid: Unique exam id
    :param series: A list of 3D volume images representing a 4D DSC series
    :return: Paths to the final .nrrd file that was created
    """
    nii_path = os.path.join(DERIVED_ROOT, uid, NII_DIR)
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

    # Delete temporary nifti files, but don't delete the skull-stripped ones
    # (we'll need those to convert to a MultiVolume in Slicer)
    for f in in_files:
        os.unlink(f)

###################################################################################################
#                                                                                                 #
# STEP 3: In slicer, import nifti as MultiVolume and save the file as FINAL_PERFUSION_MULTIVOLUME #
#                                                                                                 #
###################################################################################################


# Step 4
def correct_and_copy_dicom_tags_to_final_multivolume(uid: str) -> None:
    """
    Copies the image metadata needed for DSC module to be able to understand the MultiVolume from the original perfusion
    file
    :param uid: Exam unique id
    :return: None
    """
    print(f'Copying image metadata from original MultiVolume for {uid}...')
    orig = sitk.ReadImage(os.path.join(DERIVED_ROOT, uid, ORIGINAL_PERFUSION_MULTIVOLUME))
    final = sitk.ReadImage(os.path.join(DERIVED_ROOT, uid, f'{uid}.nrrd'))

    for k in orig.GetMetaDataKeys():
        final.SetMetaData(k, orig.GetMetaData(k))

    # For some reason stripping slightly fucks up the image origin which causes DSC module to error out with a
    # "Inputs do not occupy the same physical space!" message. Set origin to original image.
    print('Correcting origin, spacing and direction...')
    final.SetOrigin(orig.GetOrigin())
    final.SetSpacing(orig.GetSpacing())
    final.SetDirection(orig.GetDirection())

    sitk.WriteImage(final, os.path.join(DERIVED_ROOT, uid, FINAL_PERFUSION_MULTIVOLUME))


def run_perfusion(uid: str) -> None:
    """
    Uses DSCMRIAnalysis Slicer CLI module (https://www.slicer.org/w/index.php/Documentation/Nightly/Modules/DSC_MRI_Analysis)
    to extract perfusion maps for all examinations using the final segmentations created by the step above as ROI.

    This is an example of the command that will run in a subprocess:
    /Applications/Slicer.app/Contents/Extensions-30329/DSCMRIAnalysis/lib/Slicer-4.13/cli-modules/DSCMRIAnalysis \
            --aifMask /Users/kliron/Data/dsc_covid19/Derived/SEKBF00012197064/Artery.nrrd  \
            --outputK1 path_to_K1.nrrd \
            --outputK2 path_to_K2.nrrd \
            --outputMTT path_to_MTT.nrrd \
            --outputCBF path_to_CBF.nrrd \
            /Users/kliron/Downloads/multivolume.nrrd

    :param uid: Unique exam ID string
    :return: None
    """
    print(f'Computing perfusion maps for {uid}')
    root_path = os.path.join(DERIVED_ROOT, uid)
    args = ['--aifMask', os.path.join(root_path, 'Artery.nrrd'),
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


def get_statistics(write: bool = False) -> pd.DataFrame:
    """
    Extract statistics from the perfusion maps.
    :param write: If True, will write an excel file of the results at PLEXUS_STATS_FILE
    :return: Returns the Pandas dataframe of results
    """
    stats = pd.DataFrame()
    for uid in paths.keys():
        print(f'Reading perfusion maps for {uid}')
        img_root = os.path.join(DERIVED_ROOT, uid)
        # By default GetArrayFromImage returns float32 values which may cause precision errors. We convert to float64.
        k1, k2, mtt, cbf = (np.float64(sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(img_root, f)))) for f in [
            'K1.nrrd', 'K2.nrrd', 'MTT.nrrd', 'CBF.nrrd'])
        plexus_segm = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(img_root, 'final_Plexus.nrrd')))

        # Get a flat array of pixel intensity values for the whole brain and the plexus
        # A few voxels may have NaN values after skull-stripping we set them to 0.0 by callin np.nan_to_num()
        brain_values = {name: arr for arr, name in zip([k1, k2, mtt, cbf], ['K1', 'K2', 'MTT', 'CBF'])}
        plexus_values = {name: arr[plexus_segm == 1] for arr, name in zip([k1, k2, mtt, cbf], ['K1', 'K2', 'MTT', 'CBF'])}

        # Descriptive statistics
        ustats = {'uid': uid, 'n_voxels_brain': brain_values['K1'].size, 'n_voxels_plexus': plexus_values['K1'].size}
        for k, v in plexus_values.items():
            # Get the mean brain value and normalize the mean plexus value
            brain_vals = brain_values[k]
            brain_mean = np.mean(brain_vals)
            ustats[f'{k}'] = np.mean(v)
            ustats[f'{k}_norm'] = ustats[f'{k}']/brain_mean
            ustats[f'{k}_std'] = np.std(v)
            ustats[f'{k}_min'] = np.min(v)
            ustats[f'{k}_max'] = np.max(v)

        stats = stats.append(ustats, ignore_index=True)

    # Read optic sheath and mars data
    data = pd.read_csv('/Users/kliron/Downloads/sheath.csv')
    df = stats.merge(data, how='left', on='uid')
    df['CBV_norm'] = df['CBF_norm'] * df['MTT_norm']
    df['CBV'] = df['CBF'] * df['MTT']
    df = df[['uid', 'optic_sheath', 'mars',
             'CBF_norm', 'CBV_norm', 'MTT_norm', 'K1_norm', 'K2_norm',
             'CBF', 'CBV', 'MTT', 'K1', 'K2',
             'CBF_std', 'MTT_std', 'K1_std', 'K2_std',
             'CBF_max', 'MTT_max', 'K1_max', 'K2_max',
             'CBF_min', 'MTT_min', 'K1_min', 'K2_min',
             'n_voxels_plexus', 'n_voxels_brain']]
    if write:
        df.to_excel(PLEXUS_STATS_FILE, index=False)

    high_icp_cbf = df.loc[df['optic_sheath'] == 1, 'CBF']
    norm_icp_cbf = df.loc[df['optic_sheath'] == 0, 'CBF']
    high_icp_mtt = df.loc[df['optic_sheath'] == 1, 'MTT']
    norm_icp_mtt = df.loc[df['optic_sheath'] == 0, 'MTT']
    high_icp_cbv = high_icp_cbf * high_icp_mtt
    norm_icp_cbv = norm_icp_cbf * norm_icp_mtt
    high_icp_k1 = df.loc[df['optic_sheath'] == 1, 'K1']
    norm_icp_k1 = df.loc[df['optic_sheath'] == 0, 'K1']
    high_icp_k2 = df.loc[df['optic_sheath'] == 1, 'K2']
    norm_icp_k2 = df.loc[df['optic_sheath'] == 0, 'K2']
    t1, p1, d1 = ttest_ind(high_icp_cbv, norm_icp_cbf)
    t2, p2, d2 = ttest_ind(high_icp_cbv, norm_icp_cbv)
    t3, p3, d3 = ttest_ind(high_icp_cbv, norm_icp_mtt)
    t4, p4, d4 = ttest_ind(high_icp_k1, norm_icp_k1)
    t5, p5, d5 = ttest_ind(high_icp_k2, norm_icp_k2)

    print(f'Independent Two-sample T-test. Null hypothesis is average CBF values are equal. P-value: {p1}')
    print(f'Independent Two-sample T-test. Null hypothesis is average CBV values are equal. P-value: {p2}')
    print(f'Independent Two-sample T-test. Null hypothesis is average MTT values are equal. P-value: {p3}')
    print(f'Independent Two-sample T-test. Null hypothesis is average K1 values are equal. P-value: {p4}')
    print(f'Independent Two-sample T-test. Null hypothesis is average K2 values are equal. P-value: {p5}')

    plt.hist(high_icp_cbf)
    plt.hist(norm_icp_cbf)
    plt.close()

    return df


def clean_derived_files(wipeout: bool = False):
    """
    Deletes all machine-derived files in the DERIVED_ROOT folder. If wipeout is True, even deletes the
    f'{uid}.nrrd perfusion file
    """
    for uid in paths.keys():
        remove_files = [os.path.join(DERIVED_ROOT, uid, f) for f in ['K1.nrrd', 'K2.nrrd', 'CBF.nrrd', 'MTT.nrrd',
                                                                     'final_Plexus.nrrd',
                                                                     FINAL_PERFUSION_MULTIVOLUME,
                                                                     NII_DIR]]
        if wipeout:
            remove_files.append(os.path.join(DERIVED_ROOT, uid, f'{uid}.nrrd'))

        for f in remove_files:
            print(f'Removing {f}')
            try:
                shutil.rmtree(f) if os.path.isdir(f) else os.unlink(f)
            except FileNotFoundError:
                continue


def preprocess(uids: [str]) -> None:
    for uid in uids:
        params = paths[uid]
        series = edit_segmentation(uid, params)
        extract_brain(uid, series)


def main():
    parser = argparse.ArgumentParser(description='DSC perfusion analysis.',
                                     usage=f'{sys.argv[0]} [--uids UID1 UID2 ...] --preprocess | --perfusion\n')
    parser.add_argument('-u', '--uids', help='Specify UIDs to process', nargs='*')
    parser.add_argument('-p', '--preprocess', action='store_true', help='Pre-process images and segmentations')
    parser.add_argument('-c', '--copytags', action='store_true', help='Copy DICOM tags from original to skull-stripped 4D image')
    parser.add_argument('-r', '--perfusion', action='store_true', help='Run perfusion analysis')
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)

    uids = args.uids if args.uids else [uid for uid in paths.keys()]

    if args.preprocess:
        preprocess(uids)

    if args.copytags:
        for uid in uids:
            correct_and_copy_dicom_tags_to_final_multivolume(uid)

    if args.perfusion:
        for uid in uids:
            run_perfusion(uid)


if __name__ == '__main__':
    main()


