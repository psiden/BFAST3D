%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run excursions in R (just contrasts)
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
% REVISED:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function computeExcursions(outputPath,subject,method,contrastNbr,sliceNbr)
%% Setup

% Needed for MCMC for some additional input, like contrasts
if findstr(method,'MCMC')
    VBMethod = strcat('SVB',method(5:6));
end

PPMpThresh = 0.90;
PPMwPercThresh = .005;

%% Run
if findstr(method,'MCMC3D'); sliceNbr = 1; end;

subjStr = ['-',num2str(subject)];
resultsPath = strcat(outputPath,'sub',subjStr,'/',method);
num = strcat('000',num2str(contrastNbr));
num = num(end-3:end);
load(strcat(outputPath,'sub',subjStr,'/',VBMethod,'/SPM.mat'));
load(strcat(resultsPath,'/MCMC.mat'));
PPMwThresh = 100*PPMwPercThresh/(max(SPM.xBF.bf(:,1))/SPM.xBF.dt);

K = size(SPM.xX.X,2);
contrast = SPM.xCon(contrastNbr).c;

wVeck = zeros(size(squeeze(MCMC.b(sliceNbr).wVec2(1,:,:))));
for kk = 1:K
    wVeck = wVeck + contrast(kk) * squeeze(MCMC.b(sliceNbr).wVec2(kk,:,:));
end

% Save
save(strcat(resultsPath,'/excursionsInput',method,'.mat'),'wVeck');
        
% Run R and load results
system(sprintf('%s%s %s %s %f %f','unset DYLD_LIBRARY_PATH; Rscript runExcursions.R ',...
        outputPath,method,subjStr,PPMwThresh,PPMpThresh));
load(strcat(resultsPath,'/excursionsResults',method,'.mat'),'excuFuncMCMC');

% Save nii
vol = spm_vol(strcat(resultsPath,'/con_PPM_0001','.nii'));
[PPMMCMC,XYZ] = spm_read_vols(vol);
PPMMCMC(~isnan(PPMMCMC)) = 0;
bmask = ~isnan(PPMMCMC);

if findstr(method,'MCMC2D');
    
    mask2d = ~isnan(PPMMCMC(:,:,MCMC.a.sliceNbrs(sliceNbr)));
    tempSlice = nan(size(mask2d));
    tempSlice(mask2d) = excuFuncMCMC;
    PPMMCMC(:,:,MCMC.a.sliceNbrs(sliceNbr)) = tempSlice;
    
elseif findstr(method,'MCMC3');
        
    PPMMCMC(bmask) = excuFuncMCMC;
    
end

PPMMCMCThresh = PPMMCMC;
PPMMCMCThresh(bmask) = PPMMCMC(bmask) .* (PPMMCMC(bmask)>PPMpThresh);

nii = load_untouch_nii(strcat(resultsPath,'/con_PPM_',num,'.nii'));
nii.img = PPMMCMC;
save_untouch_nii(nii,strcat(resultsPath,'/con_jointPPM_',num,'.nii'));
nii.img = PPMMCMCThresh;
save_untouch_nii(nii,strcat(resultsPath,'/con_jointPPMThresh_',num,'.nii'));
