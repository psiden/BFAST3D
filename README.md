# BFAST3D
This document describes the **B**ayesian **F**ast **A**ccurate **S**patial **T**ricks in **3D** (BFAST3D) code
that can be used to run the Spatial variational Bayes (SVB) and Markov chain Monte Carlo
(MCMC) methods for fMRI analysis in Sidén et al. (2017a). The code is an add-on to the
SPM12 software (Ashburner et al., 2013) and to its Bayesian single subject method (Penny
et al., 2003, 2005b,a, 2007; Penny and Flandin, 2005). Posterior standard deviations and PPMs
are computed using the simple RBMC method described in Sidén et al. (2017b). To use the
code, follow these steps

(1) Download/duplicate your spm12-directory (downloadable at
http://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

(2) Delete the spm_spm_vb.m-file in the new spm12-directory.

(3) Copy the svb-folder from this package into the new spm12-directory.

runExample.m gives an example on how to run the code and calls the functions runSVB.m
and runMCMC.m which are adapted for the OpenfMRI data ds105, available at
https://openfmri.org/dataset/ds000105/. The current version supports BIDS-format (ds105
v.2.0.2), but requires the condition and contrast files (condition_key.txt and task_contrasts.txt)
from an earlier version (ds105 v.1.0.1) Use at own risk.

In the current version, various settings are defined in different files and some of these are
described below:
* To not run SVB in parallel, in spm12/svb/2D/spm_spm_vb.m or
spm12/svb/3D/spm_spm_vb.m, change the variable SPM.ParallelGMRFSampling
to 0.
* Change the number of SVB iterations in spm12/svb/2D/spm_spm_vb.m or
spm12/svb/3D/spm_spm_vb.m, variable SPM.PPM.maxits (default 50).
* In the beginning of spm12b/svb/spm_svb_init.m one can change some SVB settings,
e.g. the number of expectation approximation samples Ns (default 100) and
the PCG tolerance d (default 10^-8).
* For the MCMC-method, settings are changed in runMCMC.m, e.g. the number of
MCMC iterations, warmup iterations, thinning factor and PCG tolerance d can be
changed by changing the variables niter, warmup, thinningFactor and PCGTol.

The code requires MATLAB_R2016a or later. Post-processing computePPMs.m depends on the Tools for NIfTI and ANALYZE imagepackage
(http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).
To compute joint PPMs, R is required and so is the excursions-package (Bolin and Lindgren,
2015), and also the R-package R.matlab. Be aware that for large data sets the time and memory requirements
can be quite demanding, especially for MCMC, which can be quick-fixed by limiting
the number of iterations.

REFERENCES
* Ashburner, J., Barnes, G., Chen, C.-c., Daunizeau, J., Moran, R., Henson, R., Glauche, V., and Phillips, C. (2013). SPM12 Manual The FIL Methods Group ( and honorary members ). Functional Imaging Laboratory, pages 475–1.
* Bolin, D. and Lindgren, F. (2015). Excursion and contour uncertainty regions for latent Gaussian models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(1):85–106.
* Penny, W. and Flandin, G. (2005). Bayesian analysis of fMRI data with spatial priors. In Proceedings of the Joint Statistical Meeting (JSM). American Statistical Association. 
* Penny, W., Flandin, G., and Trujillo-Barreto, N. (2007). Bayesian comparison of spatially regularised general linear models. Human Brain Mapping, 28(4):275–293.
* Penny, W., Kiebel, S., and Friston, K. (2003). Variational Bayesian inference for fMRI time series. NeuroImage, 19(3):727–741.
* Penny,W. D., Trujillo-Bareto, N., and Flandin, G. (2005a). Bayesian analysis of single-subject fMRI data: SPM implementation. Technical report, Wellcome Department of Imaging Neuroscience. 
* Penny, W. D., Trujillo-Barreto, N. J., and Friston, K. J. (2005b). Bayesian fMRI time series analysis with spatial priors. NeuroImage, 24(2):350–362. 
* Sidén, P., Eklund, A., Bolin, D., and Villani, M. (2017a). Fast Bayesian whole-brain fMRI analysis with spatial 3D priors. NeuroImage, 146:211–225.
* Sidén, P., Lindgren, F., Bolin, D., and Villani, M. (2017b). Efficient Covariance Approximations for Large Sparse Precision Matrices. arXiv preprint 1705.08656v1.
