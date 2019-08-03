MSDI stands for “multi-scale dipole inversion” - a new QSM approach developed by Julio Acosta-Cabronero. 

MSDI is described in detail in Acosta-Cabronero et al. Neuroimage (2018).

A brief description of MSDI outputs can be found below:

* wC.nii.gz: Phase unreliability map
* ptb_ppm2rad.txt: ppm to radians conversion factor
* ptb_rad2ppm.txt: Padians to ppm conversion factor
* ptb_ppm2rad_ref.txt: ppm to radians conversion factor (TE*B0=60 ms*T)
* ptb_rad2ppm_ref.txt: Radians to ppm conversion factor (TE*B0=60 ms*T)
* nfm_0.nii.gz: Normalised (TE*B0=60 ms*T) local field map in radians
* Mask.nii.gz: ROI mask
* W_0.nii.gz: Initial consistency weighting (whitening) map (magn/mean_magn) 
* N_std_0.nii.gz: Phase noise normalised standard deviation map (1/W_0)
* qsm_msk_wC_[N].nii.gz: Scale-wise, reliable phase mask
* SMVII[SMV_radius_N]_nfm.nii.gz: Scale-wise, SMV-complement filtered, local field map in radians
* W_fid.nii.gz: SMV-adjusted, consistency whitening map (if ptb.qsm.MSDI.smv_N_std_0 enabled this won't be used)
* edgeima.nii.gz: Image sharing edges with susceptibility map (default: magnitude)
* wG[x,y,z,sum].nii.gz: Edge masks
* qsm_wGsum_[N].nii.gz: Scale-wise, wGsum.nii.gz
* x.nii.gz: Temporary susceptibility map
* W_merit.nii.gz: MERIT-adjusted, consistency whitening map
* qsmII_[N]_MSDI_l[lambda].nii.gz: Scale-wise susceptibility map
* qsm_INTEGRAL_[N]_MSDI_l[lambda].nii.gz: Scale-wise susceptibility-sum map
* qsm_UNMASKED_[N]_MSDI_l[lambda].nii.gz: Scale-wise, unmasked susceptibility-sum map
* qsm_nfmdiff_[N]_MSDI_l[lambda].nii.gz: Scale-wise, local field differential

[N] denotes a scale number. Each scale is defined by the spherical mean value (SMV) kernel radius vector in ’ptb.qsm.MSDI.R_smv’ (see ‘ptbs’ file).

qsm_INTEGRAL_[N<N_max]_MSDI_l[lambda].nii.gz are high-pass susceptibility maps (HPSM) 

qsm_INTEGRAL_[N_max]_MSDI_l[lambda].nii.gz is the full-scale MSDI solution