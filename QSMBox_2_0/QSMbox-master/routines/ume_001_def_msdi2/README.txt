# ume_001_def_msdi2
#
#  Default QSMbox pipeline for coil-uncombined multi-GRE data
# 
#
# It consists of three steps:
#
#  STEP #1: Brain mask calculation
#
#  STEP #2: Coil-wise multi-echo combination
#
#  STEP #3: Coil combination, background removal and QSM inversion (MSDI2)
#
#
# Instructions
# ============
#
# 1. All sM*.nii and sP*.nii files must be copied/moved to this directory.
#
#  Nb. This assumes you converted nifti files from dicom using SPM's
#  'Dicom Import' utility
#
#
# 2. Import dicom file (this will be used to read in scan parameters).
#
#  mkdir dicom
#
#  cp sample dicom file to the 'dicom' dir
#
#  If all your dicom files are in one directory, you can search for it
#  printing the dicom 'Series Description' field for all files as follows:
#
#  for i in *dcm; do echo $i; mri_probedicom --i ${i} --t 0008 103e; done
#
#  Nb. mri_probedicom comes with FreeSurfer
#
#
# 3. Edit (input arguments) and run EDIT_n_RUN_ume_001_def_msdi2.sh
#
#
# Created by Julio Acosta-Cabronero

