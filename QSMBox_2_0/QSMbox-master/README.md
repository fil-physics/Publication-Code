QSMbox v2.0 (release candidate)
===============================

MATLAB software for quantitative susceptibility mapping (QSM) reconstruction.

QSMbox contains two implementations of the multi-scale dipole inversion (MSDI) 
algorithm described in Acosta-Cabronero et al. Neuroimage (2018), 
https://doi.org/10.1016/j.neuroimage.2018.07.065

N.b. MSDI2 is the latest implementation with preservation of phase noise texture 
and low susceptibility contrast features.


INSTALL
=======

1. Copy the QSMbox directory to e.g. '/home/[USERNAME]/matlab'.

2. Edit '/home/[USERNAME]/matlab/startup.m' to include: 

	    addpath /home/[USERNAME]/matlab/QSMbox
	    addpath /home/[USERNAME]/matlab/QSMbox/master

   Or add them with the Set Path utility from the Matlab GUI.

3. Select your operating system and set the SPM path in: 
   'QSMbox/bashutils/QSMbox\_local\_settings.sh', then run:
	    
	    ./QSMbox_local_settings.sh


USAGE
=====

1. Make working directory, cd to it and open MATLAB.

2. In MATLAB, run: qsmbox, or help qsmbox

    a) Select working directory

    b) Import data from a previous run?

    c) Select datatype

	    ->[1] Coil-combined single-GRE   | full pipeline | 'cse' |
	    ->[2] Coil-combined multi-GRE    | full pipeline | 'cme' |
	    ->[3] Coil-uncombined single-GRE | full pipeline | 'use' |
	    ->[4] Add-ons                    | partial proc. | 'ppd' |
	    
	    e.g. Select: 1

    d) Select preset option

            Default MSDI pipeline for combined single-echo data
             ->[1] cse_001_def_msdi
 
            Default MSDI v2 pipeline for combined single-echo data
             ->[2] cse_002_def_msdi2

	    e.g. Select: 1

    e) Select operation
    
        ->[1] Edit settings, or ->[2] Run? 

	Both operations will import the settings file, 'ptbs_XXX.m', to pwd.

	If you selected '1' (i.e. Edit settings), the file will open automatically. 
	Read through it, prepare input data (optional, you will be prompted for 
	manual input otherwise), edit 'CUSTOM SETTINGS' if necessary, then run (F5).
	
	If you uncommented ->[flag = 'preset'], your existing settings file will be
	shown as a 'PRESET OPTION' under 'ptb.pipename' next time you run 'qsmbox'.
	Do not forget to also update the function name (line 1) and 'ptb.pipename'
	in the 'DEFINE PIPELINE' section, as well as pipeline descriptions both in the 
	help dialogue (line2) and in 'DEFINE PIPELINE' (disp message).
    
	If you have several datasets acquired with the same protocol you can propagate
	(and run) the 'ptbs_XXX.m' file across several working directories. To run 
	them from the command line you can use the QSMbox_matcml.sh program in 
	'bashutils'. N.b. The 'bashutils' directory can be added to the Bash path for 
	convenience.
	

ADDITIONAL NOTES
================

The 'routines' directory contains additional (non-standard) pipelines. 
E.g. ume\_001\_def\_msdi2 is the default MSDI2 routine for coil-uncombined, 
multi-echo data.


OUTPUTS
=======

A brief description of MSDI outputs can be found below:

* wC.nii.gz: Phase unreliability map
* ptb\_ppm2rad.txt: ppm to radians conversion factor
* ptb\_rad2ppm.txt: Radians to ppm conversion factor
* ptb\_ppm2rad\_ref.txt: ppm to radians conversion factor (TE\*B0=60 ms\*T)
* ptb\_rad2ppm\_ref.txt: Radians to ppm conversion factor (TE\*B0=60 ms\*T)
* nfm\_0.nii.gz: Normalised (TE\*B0=60 ms\*T) local field map in radians
* Mask.nii.gz: ROI mask
* W\_0.nii.gz: Initial consistency weighting (whitening) map (magn/mean\_magn) 
* N\_std\_0.nii.gz: Phase noise normalised standard deviation map (1/W\_0)
* qsm\_msk\_wC\_[N].nii.gz: Scale-wise, reliable phase mask
* SMVII[SMV\_radius\_N]\_nfm.nii.gz: Scale-wise, SMV-complement filtered, local field map in radians
* W\_fid.nii.gz: SMV-adjusted, consistency whitening map (if ptb.qsm.MSDI.smv\_N\_std\_0 enabled this won't be used)
* edgeima.nii.gz: Image sharing edges with susceptibility map (default: magnitude)
* wG[x,y,z,sum].nii.gz: Edge masks
* qsm\_wGsum\_[N].nii.gz: Scale-wise, wGsum.nii.gz
* x.nii.gz: Temporary susceptibility map
* W\_merit.nii.gz: MERIT-adjusted, consistency whitening map
* qsmII\_[N]\_MSDI\_l[lambda].nii.gz: Scale-wise susceptibility map
* qsm\_INTEGRAL\_[N]\_MSDI\_l[lambda].nii.gz: Scale-wise susceptibility-sum map
* qsm\_UNMASKED\_[N]\_MSDI.nii.gz: Scale-wise, unmasked susceptibility-sum map
* qsm\_nfmdiff\_[N]\_MSDI\_l[lambda].nii.gz: Scale-wise, local field differential

[N] denotes a scale number. Each scale is defined by the spherical mean value (SMV) kernel radius vector in ’ptb.qsm.MSDI.R\_smv’ (see ‘ptbs’ file).

qsm\_INTEGRAL\_[N < N\_max]\_MSDI\_l[lambda].nii.gz are high-pass susceptibility maps (HPSM) 

qsm\_INTEGRAL\_[N\_max]\_MSDI\_l[lambda].nii.gz is the full-scale MSDI solution


DOCUMENTATION
=============

See QSMbox/doc.


SUPPORT
=======

This software was developed by Julio Acosta-Cabronero at the Wellcome Centre
for Human Neuroimaging (University College London), supported by core funding 
from the Wellcome (203147/Z/16/Z).


DISCLAIMER
==========

QSMbox is experimental software conceived and developed for academic use by Julio Acosta-Cabronero. 

QSMbox includes third-party software - for more details see CONTRIBUTIONS.md.

QSMbox is licensed under the BSD 3-clause "New" or "Revised" License, https://choosealicense.com/licenses/bsd-3-clause. 

Julio Acosta-Cabronero and all other QSMbox contributors assume no responsibility for its use by other parties, and make no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.
