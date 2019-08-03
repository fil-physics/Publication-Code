#!/bin/bash

# QSMBOX'S LOCAL SETTINGS FILE FOR LINUX/OSX.
# READ BELOW AND EDIT ACCORDINGLY.


####################
# OPERATING SYSTEM #
####################

# The 3DSRNCP best-path unwrapper used in 'prep.init_offset_corr'
# 'prep.unwr4d.srncp', 'prep.bipolar_corr' and 'unwr.srncp' has 
# been compiled from C++ code - source code and binaries here:
#
#  QSMbox/master/ptb/_3DSRNCP/
# 
# The following README file contains more information: 
#
#  QSMbox/master/ptb/_3DSRNCP/3DSRNCP_jac_README.txt 
#
# QSMbox will read in the line below to use the correct binary:

OS="OS=2.11;" ### EDIT HERE ###

# Supported operating systems:
#
# [1.00] Windows         placeholder - currently unavailable (!)
# [2.10] Debian Linux    compiled and tested on Ubuntu 16.04
# [2.11] Debian Linux    compiled and tested on Ubuntu 14.04
# [2.20] Red Hat Linux   compiled and tested on CentOS 6.9
# [2.30] OSX             compiled and tested on 10.13.6
#
# If you are on a OS different from those above, but based 
# on Debian, Red Hat or OSX, you could try the most relevant option(s).
#
# If that did not work, please read the compilation 
# instructions at:
#
#  QSMbox/master/ptb/_3DSRNCP/3DSRNCP_jac_README.txt 
#
# Windows support needs more than this compilation - this is work 
# in progress


#######
# SPM # 
#######

# QSMbox's default brain mask calculation procedure uses SPM under the bonnet.
# SPM (Statistical Parametric Mapping, http://www.fil.ion.ucl.ac.uk/spm).
#
# Provide SPM path below:

SPMPATH="addpath /root/matlab/spm12.3" ### EDIT HERE ###

# If you do not have SPM installed uncomment below:

#NOSPM=true

# N.B. QSMbox will also let you choose BET for brain mask calculation.
# You will need Linux OS/OSX and an FSL installation.
# FSL BET (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET/UserGuide).


#####################################
# Bash actions below - DO NOT EDIT: #
#####################################

LBOXDIR=${HOME}/.QSMbox

if [ ! -d "$LBOXDIR" ]; then 
	mkdir ${LBOXDIR}
fi

rm -f ~/.QSMbox/ptb_OS.m
echo "function OS = ptb_OS" > ~/.QSMbox/ptb_OS.m
echo ${OS} >> ~/.QSMbox/ptb_OS.m

rm -f ~/.QSMbox/ptb_spmpath.m
echo ${SPMPATH} > ~/.QSMbox/ptb_spmpath.m

if [ "$NOSPM" == "true" ]; then
	rm -f ~/.QSMbox/ptb_spmpath.m
fi

