# -*- coding: utf-8 -*-
import os

DATA_ROOT = 'H:/Dokument' if os.name == 'nt' else '/Users/kliron'
DICOM_ROOT = os.path.join(DATA_ROOT, 'Data/dsc_covid19/Examinations')
ANALYSIS_ROOT = os.path.join(DATA_ROOT, 'Data/dsc_covid19/Analysis')
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
DERIVED_DIR = 'derived'
EXTRACTED_PERFUSION = 'XPerf.nrrd'
# This is the file produced by reading the skull stripped nii files
# This is the file DSCMRIAnalysis module will work on
FINAL_PERFUSION_MULTIVOLUME = 'final_Perfusion.nrrd'
DSC_EXECUTABLE = '/Applications/Slicer.app/Contents/Extensions-30329/DSCMRIAnalysis/lib/Slicer-4.13/cli-modules/DSCMRIAnalysis'
PLEXUS_STATS_FILE = os.path.join(ANALYSIS_ROOT, 'plexus_stats.xlsx')
PERFUSION_DATA_FILE = os.path.join(ANALYSIS_ROOT, 'perfusion_data.xlsx')
ADDITIONAL_DATA_FILE = os.path.join(ANALYSIS_ROOT, 'data.xlsx')
ISP_ROOT = '/Users/kliron/Data/dsc_covid19/ISP'

paths = {
    '0_SEKBF00012197064': {
        'dicom_dir': '0000294F/AACBE554/AA0C1C1E/0000376F',
        'isp_dicom_root': '0000EAE1/AAAE18D4/AACD1A05',
        'isp_dicom_dirs': ['0000F728', '00003B27', '00003C27', '0000823D', '00003858'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [4, 5, 6]
    },
    '1_SEKBF00012199907': {
        'dicom_dir': '00008DF4/AA870D89/AACE36DB/0000EA45',
        'isp_dicom_root': '00001CF3/AAFE59BD/AA936F39',
        'isp_dicom_dirs': ['0000A0F7', '0000A1BD', '0000C7D2', '000029C7', '000044F7'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [7, 8, 9]
    },
    '2_SEKBF00012200759': {
        'dicom_dir': '0000124A/AA99CE0B/AA0D55EA/0000A176',
        'isp_dicom_root': '0000BD95/AA1E088C/AA41F74F',
        'isp_dicom_dirs': ['0000AE06', '0000C9CD', '0000F851', '00009BCA', '00008579'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [10, 11, 12]
    },
    '3_SEKBF00012202695': {
        'dicom_dir': '0000732E/AA82036F/AA293FC2/0000854C',
        'isp_dicom_root': '0000D276/AA4A2F2E/AABB770E',
        'isp_dicom_dirs': ['0000A2D3', '0000AB98', '00002ED9', '0000856A', '00009337'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [6, 7, 8]
    },
    '4_SEKBF00012204240': {
        'dicom_dir': '00009BC9/AA2D8A4A/AAB031AE/0000C8EB',
        'isp_dicom_root': '00009F05/AA9AB432/AAC193B9',
        'isp_dicom_dirs': ['0000BEDE', '00003A6E', '00007AA9', '000075C5', '0000667D'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [9, 10, 11]
    },
    '5_SEKBF00012207046': {
        'dicom_dir': '0000FBF7/AA51C680/AA035F69/0000CDDA',
        'isp_dicom_root': '0000A7AE/AA4DDB2C/AACC88D0',
        'isp_dicom_dirs': ['0000DFFF', '00008ABA', '00009DF2', '000027E7', '0000461E'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [8, 9, 10]
    },
    # Large ICH
    '6_SEKBF00012208882': {
        'dicom_dir': '0000CD6C/AA40BA99/AAC0EC4F/0000FC82',
        'isp_dicom_root': '0000CDD6/AAE3FAAC/AAFEFFF2',
        'isp_dicom_dirs': ['0000D65D', '0000FFE6', '0000421E', '00002510', '00006498'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [4, 5, 6]
    },
    '7_SEKBF00012212884': {
        'dicom_dir': '0000E9BE/AA112ECA/AAC9E5B4/00005935',
        'isp_dicom_root': '0000422F/AA2065DB/AAD8E667',
        'isp_dicom_dirs': ['0000A2DC', '0000F987', '00002A9D', '000054A4', '00008438'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [14, 15, 16]
    },
    '8_SEKBF00012214682': {
        'dicom_dir': '0000180D/AAC8BD9E/AAB7B29B/00006420',
        'isp_dicom_root': '0000B653/AA11024E/AA639A18',
        'isp_dicom_dirs': ['0000A374', '0000DBC4', '000080CE', '00002273', '00008808'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [16, 17, 18]
    },
    '9_SEKBF00012216532': {
        'dicom_dir': '0000AE50/AA544B16/AA57B434/0000B05B',
        'isp_dicom_root': '0000D0DD/AA879784/AA076E53',
        'isp_dicom_dirs': ['0000A713', '0000EDDB', '000040C7', '0000213C', '0000496B'],
        'noncontrast_frames': [3, 4, 5],
        'contrast_frames': [14, 15, 16]
    },
    '10_SEKBF00012217381': {
        'dicom_dir': '0000D98A/AA5060C5/AA08AD0F/00001946',
        'isp_dicom_root': '0000451D/AA7FC8DE/AA85A279',
        'isp_dicom_dirs': ['0000AF40', '00003E68', '00009CDC', '00000283', '00004488'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [16, 17, 18]
    },
    '11_SEKBF00012221435': {
        'dicom_dir': '0000229C/AA605811/AAA98CA4/00000848',
        'isp_dicom_root': '00001999/AA735A61/AA641891',
        'isp_dicom_dirs': ['00004A8A', '0000115A', '0000160B', '00003537', '00007589'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [16, 17, 18]
    },
    '12_SEKBF00012223948': {
        'dicom_dir': '00000272/AA759D6D/AA444757/00007313',
        'isp_dicom_root': '00002CCC/AAEF3724/AA3FE5FB',
        'isp_dicom_dirs': ['0000BA18', '0000D2EE', '000049F7', '0000812B', '00002509'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [8, 9, 10]
    },
    '13_SEKBF00012226261': {
        'dicom_dir': '000080AC/AA650728/AA8B09A5/000098AC',
        'isp_dicom_root': '0000C94F/AA334A42/AAD0E8D5',
        'isp_dicom_dirs': ['0000C36C', '00003F56', '00004BF5', '0000025D', '0000090E'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [8, 9, 10]
    },
    '14_SEKBF00012226431': {
        'dicom_dir': '00006E1C/AA2DF666/AA1A00EA/00006444',
        'isp_dicom_root': '00008E88/AA9FAC02/AA1580F6',
        'isp_dicom_dirs': ['0000EFF7', '0000F5FD', '0000F646', '00009A2A', '000055F5'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [9, 10, 11]
    },
    '15_SEKBF00012227656': {
        'dicom_dir': '0000AF1E/AAC44726/AAA982BD/0000CD29',
        'isp_dicom_root': '00009D54/AA38771F/AA2CA2C0',
        'isp_dicom_dirs': ['0000ADB9', '0000EF47', '0000FBA5', '000037AF', '0000314C'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [10, 11, 12]
    },
    '16_SEKBF00012228858': {
        'dicom_dir': '0000D7FA/AA4FE859/AA71394E/00004471',
        'isp_dicom_root': '000045E2/AA012ACA/AA1E9229',
        'isp_dicom_dirs': ['0000AD94', '0000C5AA', '00000C14', '00003465', '00006070'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [8, 9, 10]
    }
}