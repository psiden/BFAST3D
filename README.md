# BFAST3D
This document describes the Bayesian Fast Accurate Spatial Tricks in 3D (BFAST3D) code
that can be used to run the Spatial variational Bayes (SVB) and Markov chain Monte Carlo
(MCMC) methods for fMRI analysis in Sidén et al. (2016). The code is an add-on to the
SPM12 software (Ashburner et al., 2013) and to its Bayesian single subject method (Penny
et al., 2003, 2005b,a, 2007; Penny and Flandin, 2005). To use the code, follow these steps

(1) Download/duplicate your spm12-directory (downloadable at
http://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

(2) Delete the spm_spm_vb.m-file in the new spm12-directory.

(3) Copy the svb-folder from this package into the new spm12-directory.

runExample.m gives an example on how to run the code and calls the functions runSVB.m
and runMCMC.m which are adapted for the OpenfMRI data ds105, available at
https://openfmri.org/dataset/ds000105/. Use at own risk.

In the current version, various settings are defined in different files and some of these are
described below:
* To not run SVB in parallel, in spm12/svb/2D/spm_spm_vb.m or
spm12/svb/3D/spm_spm_vb.m, change the variable SPM.ParallelGMRFSampling
to 0.
* Change the number of SVB iterations in spm12/svb/2D/spm_spm_vb.m or
spm12/svb/3D/spm_spm_vb.m, variable SPM.PPM.maxits (default 50).
* In the beginning of spm12b/svb/spm_svb_init.m one can change some SVB settings,
e.g. the number of expectation approximation samples Ns (default 100) and
the PCG tolerance d (default 10􀀀8).
* For the MCMC-method, settings are changed in runMCMC.m, e.g. the number of
MCMC iterations, warmup iterations, thinning factor and PCG tolerance d can be
changed by changing the variables niter, warmup, thinningFactor and PCGTol.

Post-processing computePPMs.m depends on the Tools for NIfTI and ANALYZE imagepackage
(http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).
To compute joint PPMs, R is required and so is the excursions-package (Bolin and Lindgren,
2015), development version (https://bitbucket.org/davidbolin/excursions) and also
the R-package R.matlab. Be aware that for large data sets the time and memory requirements
can be quite demanding, especially for MCMC, which can be quick-fixed by limiting
the number of iterations.
