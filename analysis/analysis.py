#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os.path
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
import shutil
import re
from statsmodels.stats.weightstats import ttest_ind

# This is to make relative imports work both when this file is run "on the fly" in an IDE
# and when it is run as a script (as __main__). For the first case to work, the console
# must of course cd in the directory this file is in.
parent_module = sys.modules['.'.join(__name__.split('.')[:-1]) or '__main__']
if __name__ == '__main__' or parent_module.__name__ == '__main__':
    from util.paths import *
    from util.util import *
else:
    # Create an empty CMakeLists.txt file and reload it.
    # In order for CLion to be able to import local python modules you have to mark the directory containing them as
    # "sources". There is a bug in Clion where the context menu is missing "mark directory as". To workaround create an
    # empty CMakeLists.txt file, reload the project and the option will appear.
    from .util.paths import *
    from .util.util import *

os.environ['SITK_SHOW_COMMAND'] = 'H:/Slicer 4.11.20210226/Slicer.exe' if os.name == 'nt' else '/Applications/Slicer.app/Contents/MacOS/Slicer'
os.environ['FSLOUTPUTTYPE'] = 'NIFTI'  # this is to make bet2 tool not bitch when we execute it with subprocess


# Step 1
def edit_segmentation(uid: str,
                      foreground_label: int = 1,
                      background_label: int = 0,
                      show: bool = False) -> [sitk.Image]:
    """
    We remove voxels that contain calcifications or don't take up gadolinium from the plexus segmentation.
    Segmentations are saved as labelmaps where each voxel belonging to the segmentation has the value 1 (default)
    and the rest the value 0. Thus, by using orig_img_arr[segm == 1] we can index into the original image and get the
    voxels corresponding to the volume of the segmentation.
    :param uid: Unique identification number for an exam.
    :param foreground_label: Foreground integer value of the segmentation LabelMap (Default 1).
    :param background_label: Background integer value of the segmentation LabelMap (Default 0).
    :param show: If true, show an overlay of the final segmentation and the first slice of the first frame of the series.
    :return: The 3D image list to process in other functions downstream
    """
    print(f'Editing segmentation for {uid}')
    analysis_root = os.path.join(ANALYSIS_ROOT, uid)
    params = paths[uid]
    dcm_path = os.path.join(DICOM_ROOT, uid, 'DICOM', params['dicom_dir'])
    series, dicom_tags = read_4d_series(dcm_path)
    segm = sitk.ReadImage(os.path.join(analysis_root, SEGMENTATION_FILENAME))

    # Get the average voxel values for the non-contrast series
    noncontrast_imgs = [series[frame] for frame in params['noncontrast_frames']]
    contrast_imgs = [series[frame] for frame in params['contrast_frames']]
    noncontrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in noncontrast_imgs]) / len(noncontrast_imgs)).astype(int)
    contrast_avg_array = np.round(np.add(*[sitk.GetArrayFromImage(img) for img in contrast_imgs]) / len(contrast_imgs)).astype(int)

    segm_arr = sitk.GetArrayFromImage(segm)
    # Get a list of all voxel indexes that belong to the segmentation
    mask_indexes = list(zip(*np.where(segm_arr == foreground_label)))

    # Remove voxels where the absolute value of the intensity difference between noncontrast minus contrast average
    # intensity is lower than one standard deviation of the
    # average noncontrast intensity.
    remove_condition = np.abs((noncontrast_avg_array[segm_arr == 1] - contrast_avg_array[segm_arr == 1])) < np.std(noncontrast_avg_array[segm_arr == 1])

    # Setting the voxel to 0 (the default label for background) removes it from the segmentation
    for idx, cond in zip(mask_indexes, remove_condition):
        if cond:
            segm_arr[idx] = background_label

    # The final segmentation contains voxels that take up contrast and don't contain calcifications
    final_segm = sitk.GetImageFromArray(segm_arr)
    final_segm.CopyInformation(segm)
    show_overlay(noncontrast_imgs[0], final_segm) if show else ()
    derived_dir = os.path.join(analysis_root, DERIVED_DIR)
    if not os.path.exists(derived_dir):
        os.mkdir(derived_dir)
    sitk.WriteImage(final_segm, os.path.join(derived_dir, f'final_{SEGMENTATION_FILENAME}'))
    return series


# Step 2
def extract_brain(uid: str) -> None:
    """
    For each of N 3D volumes in `series` a NIFTI file will be saved under ANALYSIS_ROOT/${uid}/tmpdir/ directory.
    BET will be run on N parallel subprocesses to extract the brain and save a new NIFTI file. Each of the new files
    will be joined in a 4D image with JoinSeriesImageFilter which will be saved as a .nrrd file.
    :param uid: Unique exam id
    :return: Paths to the final .nrrd file that was created
    """
    params = paths[uid]
    dcm_path = os.path.join(DICOM_ROOT, uid, 'DICOM', params['dicom_dir'])
    series, dicom_tags = read_4d_series(dcm_path)
    nii_path = os.path.join(ANALYSIS_ROOT, uid, DERIVED_DIR, NII_DIR)
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

#####################################################################################
#                                                                                   #
# STEP 3: In slicer, import nifti as MultiVolume and save the file as EXT_PERFUSION #
#                                                                                   #
#####################################################################################


# Step 4
def correct_and_copy_dicom_tags_to_final_multivolume(uid: str) -> None:
    """
    1. Copies the image metadata needed for DSC module to be able to understand the MultiVolume from the original perfusion
       file.
    2. Copies the Origin and Spacing information from the original image. That is because skull-stripping with BET
       causes a miniscule (in the order of 0.000001 cm) distortion of the image which is still above Tolerance level
       for ITK's ImageFilter class and causes DSC module to error out with a "Inputs do not occupy the same physical
       space!" message.
       As the AIF function is *crucial* for correct perfusion analysis, make sure to select the Artery segmentation on
       the FINAL image that DSC will process. Otherwise you will have to set the Origin and Spacing even in the
       segmentation!
    :param uid: Exam unique id
    :return: None
    """
    print(f'Copying image metadata from original MultiVolume for {uid}...')
    orig = sitk.ReadImage(os.path.join(ANALYSIS_ROOT, uid, ORIGINAL_PERFUSION_MULTIVOLUME))
    ext = sitk.ReadImage(os.path.join(ANALYSIS_ROOT, uid, DERIVED_DIR, EXTRACTED_PERFUSION))
    for k in orig.GetMetaDataKeys():
        ext.SetMetaData(k, orig.GetMetaData(k))

    print('Correcting origin and spacing...')
    ext.SetOrigin(orig.GetOrigin())
    ext.SetSpacing(orig.GetSpacing())
    sitk.WriteImage(ext, os.path.join(ANALYSIS_ROOT, uid, DERIVED_DIR, FINAL_PERFUSION_MULTIVOLUME))


def run_perfusion(uid: str, show=False) -> None:
    """
    Uses DSCMRIAnalysis Slicer CLI module (https://www.slicer.org/w/index.php/Documentation/Nightly/Modules/DSC_MRI_Analysis)
    to extract perfusion maps for all examinations using the final segmentations created by the step above as ROI.

    This is an example of the command that will run in a subprocess:
    /Applications/Slicer.app/Contents/Extensions-30329/DSCMRIAnalysis/lib/Slicer-4.13/cli-modules/DSCMRIAnalysis \
            --aifMask /Users/kliron/Data/dsc_covid19/0_SEKBF00012197064/Artery.nrrd  \
            --outputK1 /Users/kliron/Data/dsc_covid19/Derived/0_SEKBF00012197064/K1.nrrd \
            --outputK2 /Users/kliron/Data/dsc_covid19/Derived/0_SEKBF00012197064/K2.nrrd \
            --outputMTT /Users/kliron/Data/dsc_covid19/Derived/0_SEKBF00012197064/MTT.nrrd \
            --outputCBF /Users/kliron/Data/dsc_covid19/Derived/0_SEKBF00012197064/CBF.nrrd \
            /Users/kliron/Downloads/multivolume.nrrd

    :param uid: Unique exam ID string
    :param show: If True, show resulting CBF map
    :return: None
    """
    print(f'Computing perfusion maps for {uid}')
    fpath = os.path.join(ANALYSIS_ROOT, uid)
    args = ['--aifMask', os.path.join(fpath, 'Artery.nrrd'),
            '--outputK1', os.path.join(fpath, DERIVED_DIR, 'K1.nrrd'),
            '--outputK2', os.path.join(fpath, DERIVED_DIR, 'K2.nrrd'),
            '--outputMTT', os.path.join(fpath, DERIVED_DIR, 'MTT.nrrd'),
            '--outputCBF', os.path.join(fpath, DERIVED_DIR, 'CBF.nrrd'),
            os.path.join(fpath, DERIVED_DIR, FINAL_PERFUSION_MULTIVOLUME)]

    print('Executing DSC CLI module...')
    print(DSC_EXECUTABLE + ' \\ \n' +
          ' \\ \n    '.join([f'{args[i]} {args[j]}' for i, j in zip(range(0, len(args), 2), range(1, len(args), 2))]) +
          ' \\ \n    ' + args[len(args) - 1])
    exit_codes = run_n_subprocesses([DSC_EXECUTABLE + ' ' + ' '.join(args)])
    if any(exit_codes):
        print('Warning: some of the subprocesses returned nonzero exit status')
    if show:
        img = sitk.ReadImage(os.path.join(fpath, DERIVED_DIR, 'CBF.nrrd'))
        sitk.Show(img, title=f'{uid}, CBF map')


def get_plexus_measurements(use_isp_maps=True, uids: List[str] = None) -> pd.DataFrame:
    """
    Extract statistics from the perfusion maps.
    :param use_isp_maps: If true (default) will use the Philips ISP perfusion maps for all measurements
    :param uids: get statistics only for these UIDs
    :return: Returns the Pandas dataframe of results
    """
    stats = pd.DataFrame()
    for uid in paths.keys():
        if uids and uid not in uids:
            continue
        print(f'Reading perfusion maps for {uid}')
        img_root = os.path.join(ANALYSIS_ROOT, uid) if use_isp_maps else os.path.join(ANALYSIS_ROOT, uid, DERIVED_DIR)
        # By default GetArrayFromImage returns float32 values which may cause precision errors. We convert to float64.
        k1, k2, mtt, cbf, cbv = (np.float64(sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(img_root, f)))) for f in [
            'K1.nrrd', 'K2.nrrd', 'MTT.nrrd', 'CBF.nrrd', 'CBV.nrrd'])
        plexus_segm = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(ANALYSIS_ROOT, uid, DERIVED_DIR, 'final_Plexus.nrrd')))

        # Get a flat array of pixel intensity values for the whole brain and the plexus
        # A few voxels may have NaN values after skull-stripping we set them to 0.0 by calling np.nan_to_num()
        brain_values = {name: arr for name, arr in zip(['K1', 'K2', 'MTT', 'CBF', 'CBV'], [k1, k2, mtt, cbf, cbv])}
        plexus_values = {name: arr[plexus_segm == 1] for name, arr in zip(['K1', 'K2', 'MTT', 'CBF', 'CBV'], [k1, k2, mtt, cbf, cbv])}

        # Descriptive statistics
        ustats = {'uid': uid, 'n_voxels_brain': brain_values['K1'].size, 'n_voxels_plexus': plexus_values['K1'].size}
        for k, v in plexus_values.items():
            # Get the mean brain value and normalize the mean plexus value
            brain_vals = brain_values[k]
            # Some voxels may have NaN values after processing, that is why we call np.nanmean() instead of np.mean()
            brain_mean = np.nanmean(brain_vals)
            ustats[f'{k}'] = np.mean(v)
            ustats[f'{k}_std'] = np.std(v)
            ustats[f'{k}_min'] = np.min(v)
            ustats[f'{k}_max'] = np.max(v)
            # Normalized by the whole brain mean
            ustats[f'{k}_norm'] = ustats[f'{k}']/brain_mean
        stats = stats.append(ustats, ignore_index=True)

    return stats[['uid', 'CBF_norm', 'CBV_norm', 'MTT_norm', 'K1_norm', 'K2_norm', 'CBF', 'CBV', 'MTT', 'K1', 'K2',
                  'CBF_std', 'CBV_std', 'MTT_std', 'K1_std', 'K2_std', 'CBF_max', 'CBV_max', 'MTT_max', 'K1_max', 'K2_max',
                  'CBF_min', 'CBV_min', 'MTT_min', 'K1_min', 'K2_min', 'n_voxels_plexus', 'n_voxels_brain']]


def merge_and_save_excel(stats: pd.DataFrame) -> None:
    # Read optic sheath and mars data
    data = pd.read_excel(ADDITIONAL_DATA_FILE)
    d1 = pd.read_excel('/Users/kliron/Data/covid19/Radiology2020/clinical.xlsx')[['pid', 'days_on_ventilator', 'ventilator_off']]
    d2 = pd.read_excel('/Users/kliron/Data/covid19/Radiology2020/radiology_old.xlsx')[['pid', 'unr', 'scan_date']]
    d3 = d1.merge(d2, how='inner', on='pid')
    df = stats.merge(data, how='left', on='uid')
    df['unr'] = df['uid'].str.split(r'[0-9]{1,2}_').str[-1]
    df = df.merge(d3, how='left', on='unr')
    df['days_from_extubation_to_scan'] = df['scan_date'] - df['ventilator_off']
    df = df[['pid', 'uid', 'optic_sheath', 'mars', 'iva', 'days_on_ventilator', 'days_from_extubation_to_scan',
             'CBF_norm', 'CBV_norm', 'MTT_norm', 'K1_norm', 'K2_norm', 'CBF', 'CBV', 'MTT', 'K1', 'K2',
             'CBF_std', 'CBV_std', 'MTT_std', 'K1_std', 'K2_std', 'CBF_max', 'CBV_max', 'MTT_max', 'K1_max', 'K2_max',
             'CBF_min', 'CBV_min', 'MTT_min', 'K1_min', 'K2_min', 'n_voxels_plexus', 'n_voxels_brain']]
    df.to_excel(PLEXUS_STATS_FILE, index=False)


def clean_derived_files():
    """Deletes all machine-derived files in the ANALYSIS_ROOT."""
    if input(f'This will delete ALL "{DERIVED_DIR}" directories under {os.path.join(DATA_ROOT, ANALYSIS_ROOT) + "/{uid}/"}, are you sure? (Y/N)') != 'Y':
        exit(0)
    for uid in paths.keys():
        derived = os.path.join(ANALYSIS_ROOT, uid, DERIVED_DIR)
        if os.path.exists(derived):
            shutil.rmtree(derived)


def save_nrrd_from_isp_dicom() -> None:
    """Reads the DICOM files from Philips ISP perfusion results. Asks user for name of sequence to save."""
    imply_yes = False
    for uid, vals in paths.items():
        dcm_root = os.path.join(DICOM_ROOT, uid, 'DICOM', '/'.join(vals['dicom_dir'].split('/')[:-1]))
        dcm_dirs = vals.get('isp_dicom_dirs', None)
        dcm_dirs = os.listdir(dcm_root) if dcm_dirs is None else dcm_dirs
        for dcm_dir in dcm_dirs:
            dcm_path = os.path.join(dcm_root, dcm_dir)
            print(f'Trying to read 3D series at {dcm_path}')
            try:
                img, tags = read_3d_series(dcm_path)
            except RuntimeError as e:
                print(f'Error trying to read 3D DICOM series: {e}')
                continue

            # ISP saves the type of image under 'Image Type' tag
            img_type = tags['Image Type']
            m = re.match(r'.*?(MTT|CBF|CBV|K1|K2).*', img_type)
            if not m or not m.group(1):
                print(f'Could not find any matching perfusion type for {uid}.')
                print(f'Image type information in the DICOM was: {img_type}')
                print('Skipping...')
                continue
            else:
                ptype = m.group(1)
                out_path = os.path.join(DATA_ROOT, ANALYSIS_ROOT, uid, f'{ptype}.nrrd')
                print(f'Found the type {ptype}. Will save at {out_path}')
                if not imply_yes:
                    choice = input('Is this correct [y/n/A]?')
                    imply_yes = choice == 'A'
                    if choice == 'y':
                        sitk.WriteImage(img, out_path)
                    else:
                        print('Skipping...')
                else:
                    sitk.WriteImage(img, out_path)


def run_stats() -> None:
    df = pd.read_excel('/Users/kliron/Data/covid19/plexus_perfusion/Analysis/final_analysis.xlsx')
    print(f'Cases # = {len(df.loc[df["case"] == 1,])} ({len(df.loc[(df["case"] == 1) & (df["sex"] == "F"),])} female)')
    print(f'Controls # = {len(df.loc[df["case"] == 0,])} ({len(df.loc[(df["case"] == 0) & (df["sex"] == "F"),])} female)')
    print(f'Case age mean: {df.loc[df["case"] == 1, "age"].mean():.1f}, median: {df.loc[df["case"] == 1, "age"].median():.1f}, '
          f'standard deviation: {df.loc[df["case"] == 1, "age"].std():.1f}')
    print(f'Control age mean: {df.loc[df["case"] == 0, "age"].mean():.1f}, median: {df.loc[df["case"] == 0, "age"].median():.1f}, '
          f'standard deviation: {df.loc[df["case"] == 0, "age"].std():.1f}')

    selectors = {
        'high_icp': (df['case'] == 1) & (df['optic_sheath'] == 1),
        'norm_icp': (df['case'] == 1) & (df['optic_sheath'] == 0),
        'ctrl': df['case'] == 0,
        'ctrl_old': (df['case'] == 0) & (df['age'] > 40),
        'ctrl_young': (df['case'] == 0) & (df['age'] < 40)
    }
    comparison_pairs = [('high_icp', 'norm_icp'),
                        ('high_icp', 'ctrl'), ('norm_icp', 'ctrl'),
                        ('high_icp', 'ctrl_young'), ('norm_icp', 'ctrl_young'),
                        ('high_icp', 'ctrl_old'), ('norm_icp', 'ctrl_old'),
                        ('ctrl_young', 'ctrl_old')]
    comparison_vars = ['CBF_norm', 'CBV_norm', 'MTT_norm', 'K1_norm', 'K2_norm']

    for left, right in comparison_pairs:
        sel_left = selectors[left]
        sel_right = selectors[right]
        print('')
        for var in comparison_vars:
            rows_left = df.loc[sel_left, var]
            rows_right = df.loc[sel_right, var]
            t, p, d = ttest_ind(rows_left, rows_right)
            print(f'{var.replace("_norm", "")} {left.replace("_", " ")} vs {right.replace("_", " ")}, p = {p:.3f} {"**" if p < 0.01 else "*" if p < 0.05 else ""}')


def main():
    parser = argparse.ArgumentParser(description='DSC perfusion analysis.', usage=f'{sys.argv[0]} [--uids UID1 UID2 ...] --preprocess | --perfusion\n')
    parser.add_argument('-u', '--uids', help='Specify UIDs to process', nargs='*')
    parser.add_argument('-e', '--edit_segmentations', action='store_true', help='Remove voxels with calcifications and voxels that do not take up contrast from plexus segmentations')
    parser.add_argument('-b', '--extract_brain', action='store_true', help='Call BET tool to extract brain from perfusion series')
    parser.add_argument('-c', '--copytags', action='store_true', help='Copy DICOM tags from original to skull-stripped 4D image')
    parser.add_argument('-r', '--perfusion', action='store_true', help='Run perfusion analysis')
    parser.add_argument('-x', '--clean', action='store_true', help='Clean all derived data')
    parser.add_argument('-i', '--isp', action='store_true', help='Save DICOM ISP perfusion maps as .nrrd')
    parser.add_argument('-a', '--analysis', action='store_true', help='Get plexus measurements')
    parser.add_argument('-s', '--stats', action='store_true', help='Run statistics tests')
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)

    uids = args.uids if args.uids else [uid for uid in paths.keys()]

    if args.clean:
        clean_derived_files()

    if args.isp:
        save_nrrd_from_isp_dicom()

    if args.edit_segmentations:
        for uid in uids:
            edit_segmentation(uid)

    if args.extract_brain:
        for uid in uids:
            extract_brain(uid)

    if args.copytags:
        for uid in uids:
            correct_and_copy_dicom_tags_to_final_multivolume(uid)

    if args.perfusion:
        for uid in uids:
            run_perfusion(uid, show=len(args.uids) != 0)

    if args.analysis:
        get_plexus_measurements()

    if args.stats:
        run_stats()


if __name__ == '__main__':
    main()
