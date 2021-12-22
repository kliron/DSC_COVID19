# -*- coding: utf-8 -*-
import os

DATA_ROOT = 'H:/Dokument' if os.name == 'nt' else '/Users/kliron/Data/covid19/perfusion/data'
ANALYSIS_ROOT = os.path.join(DATA_ROOT, '/Users/kliron/Data/covid19/perfusion/analysis')
PLEXUS_SEGMENTATION_FILENAME = 'Plexus.nrrd'
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

paths = {
    'SEKBF00012197064': {
        'case': 1,
        'perf_dicom_dir': '0000376F',
        'isp_dicom_dirs': ['0000F728', '00003B27', '00003C27', '0000823D', '00003858'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [4, 5, 6]
    },
    'SEKBF00012199907': {
        'case': 1,
        'perf_dicom_dir': '0000EA45',
        'isp_dicom_dirs': ['0000A0F7', '0000A1BD', '0000C7D2', '000029C7', '000044F7'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [7, 8, 9]
    },
    'SEKBF00012200759': {
        'case': 1,
        'perf_dicom_dir': '0000A176',
        'isp_dicom_dirs': ['0000AE06', '0000C9CD', '0000F851', '00009BCA', '00008579'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [10, 11, 12]
    },
    # Severe motion artifacts
    # 'SEKBF00012202695': {
    #     'perf_dicom_dir': '0000854C',
    #     'isp_dicom_dirs': ['0000A2D3', '0000AB98', '00002ED9', '0000856A', '00009337'],
    #     'noncontrast_frames': [0, 1, 2],
    #     'contrast_frames': [6, 7, 8]
    # },
    'SEKBF00012204240': {
        'case': 1,
        'perf_dicom_dir': '0000C8EB',
        'isp_dicom_dirs': ['0000BEDE', '00003A6E', '00007AA9', '000075C5', '0000667D'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [9, 10, 11]
    },
    'SEKBF00012207046': {
        'case': 1,
        'perf_dicom_dir': '0000CDDA',
        'isp_dicom_dirs': ['0000DFFF', '00008ABA', '00009DF2', '000027E7', '0000461E'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [8, 9, 10]
    },
    # Large ICH
    # 'SEKBF00012208882': {
    #     'perf_dicom_dir': '0000FC82',
    #     'isp_dicom_dirs': ['0000D65D', '0000FFE6', '0000421E', '00002510', '00006498'],
    #     'noncontrast_frames': [0, 1, 2],
    #     'contrast_frames': [4, 5, 6]
    # },
    'SEKBF00012212884': {
        'case': 1,
        'perf_dicom_dir': '00005935',
        'isp_dicom_dirs': ['0000A2DC', '0000F987', '00002A9D', '000054A4', '00008438'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [14, 15, 16]
    },
    'SEKBF00012214682': {
        'case': 1,
        'perf_dicom_dir': '00006420',
        'isp_dicom_dirs': ['0000A374', '0000DBC4', '000080CE', '00002273', '00008808'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [16, 17, 18]
    },
    'SEKBF00012216532': {
        'case': 1,
        'perf_dicom_dir': '0000B05B',
        'isp_dicom_dirs': ['0000A713', '0000EDDB', '000040C7', '0000213C', '0000496B'],
        'noncontrast_frames': [3, 4, 5],
        'contrast_frames': [14, 15, 16]
    },
    'SEKBF00012217381': {
        'case': 1,
        'perf_dicom_dir': '00001946',
        'isp_dicom_dirs': ['0000AF40', '00003E68', '00009CDC', '00000283', '00004488'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [16, 17, 18]
    },
    'SEKBF00012221435': {
        'case': 1,
        'perf_dicom_dir': '00000848',
        'isp_dicom_dirs': ['00004A8A', '0000115A', '0000160B', '00003537', '00007589'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [16, 17, 18]
    },
    'SEKBF00012223948': {
        'case': 1,
        'perf_dicom_dir': '00007313',
        'isp_dicom_dirs': ['0000BA18', '0000D2EE', '000049F7', '0000812B', '00002509'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [8, 9, 10]
    },
    'SEKBF00012226261': {
        'case': 1,
        'perf_dicom_dir': '000098AC',
        'isp_dicom_dirs': ['0000C36C', '00003F56', '00004BF5', '0000025D', '0000090E'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [8, 9, 10]
    },
    'SEKBF00012226431': {
        'case': 1,
        'perf_dicom_dir': '00006444',
        'isp_dicom_dirs': ['0000EFF7', '0000F5FD', '0000F646', '00009A2A', '000055F5'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [9, 10, 11]
    },
    'SEKBF00012227656': {
        'case': 1,
        'perf_dicom_dir': '0000CD29',
        'isp_dicom_dirs': ['0000ADB9', '0000EF47', '0000FBA5', '000037AF', '0000314C'],
        'noncontrast_frames': [2, 3, 4],
        'contrast_frames': [10, 11, 12]
    },
    'SEKBF00012228858': {
        'case': 1,
        'perf_dicom_dir': '00004471',
        'isp_dicom_dirs': ['0000AD94', '0000C5AA', '00000C14', '00003465', '00006070'],
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [8, 9, 10]
    },
    # Controls
    '1MNE23': {
        'case': 0,
        'perf_dicom_dir': '0000E8AC',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [29, 30, 31],
        'optic_sheaths': {'r': 6.6, 'l': 4.7}
    },
    '2FJA20': {
        'case': 0,
        'perf_dicom_dir': '0000A465',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [14, 15, 16],
        'optic_sheaths': {'r': 7, 'l': 8}
    },
    '3FAA26': {
        'case': 0,
        'perf_dicom_dir': '00001B15',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [16, 17, 18],
        'optic_sheaths': {'r': 7.3, 'l': 7.8}
    },
    '4FSM57': {
        'case': 0,
        'perf_dicom_dir': '00007468',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [14, 15, 16],
        'optic_sheaths': {'r': 7, 'l': 5.9}
    },
    '5MSA56': {
        'case': 0,
        'perf_dicom_dir': '0000B401',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [31, 32, 33],
        'optic_sheaths': {'r': 6.7, 'l': 5.2}
    },
    '6MHJ24': {
        'case': 0,
        'perf_dicom_dir': '00000D08',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [28, 29, 30],
        'optic_sheaths': {'r': 7.2, 'l': 6.3}
    },
    '7MGP54': {
        'case': 0,
        'perf_dicom_dir': '0000B5C7',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [28, 29, 30],
        'optic_sheaths': {'r': 6, 'l': 4.9}
    },
    '8FSM28': {
        'case': 0,
        'perf_dicom_dir': '0000C729',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [25, 26, 27],
        'optic_sheaths': {'r': 5.9, 'l': 5}
    },
    '9FLL23': {
        'case': 0,
        'perf_dicom_dir': '0000FF5A',
        'noncontrast_frames': [0, 1, 2],
        'contrast_frames': [27, 28, 29],
        'optic_sheaths': {'r': 4.7, 'l': 4.1}
    }
}
