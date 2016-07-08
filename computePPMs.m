%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Compute Posterior mean and PPMs for all contrasts and
%               regressors for a certain method. Save result as nifti.
%
%               Assuming using design matrix with canonical HRF and
%               temporal derivative + 6 motion regressors and intercept.
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
function computePPMs(outputPath,subject,method)
%% Setup

% Needed for MCMC for some additional input, like contrasts
if findstr(method,'MCMC')
    VBMethod = strcat('SVB',method(5:6));
end

PPMpThresh = 0.90;
PPMwPercThresh = .005;

%% Run

subjStr = strcat('00',num2str(subject));
subjStr = subjStr(end-2:end);
resultsPath = strcat(outputPath,'sub',subjStr,'/',method);

if findstr(method,'VB')
    
    load(strcat(resultsPath,'/SPM.mat'));
    PPMwThresh = 100*PPMwPercThresh/(max(SPM.xBF.bf(:,1))/SPM.xBF.dt);
    K = size(SPM.xX.X,2);
    nContr = length(SPM.xCon);
    
    % Compute PPMs for relevant regressors (exclude temporal derivative
    % hrf, hmp and intercept regressors)
    for k = 1:2:(K-7)
        num = strcat('000',num2str(k));
        num = num(end-3:end);
        
        % Load
        vol = spm_vol(strcat(resultsPath,'/Cbeta_',num,'.nii'));
        [wVB,XYZ] = spm_read_vols(vol);
        
        % Load stds
        vol = spm_vol(strcat(resultsPath,'/SDbeta_',num,'.nii'));
        [stdVB,XYZ] = spm_read_vols(vol);
        
        % Compute marginal PPM
        PPMVB = nan(size(wVB));
        PPMVBThresh = nan(size(wVB));
        bmask = ~isnan(wVB);
        PPMVB(bmask) = normcdf((wVB(bmask)-PPMwThresh) ./ stdVB(bmask));
        PPMVBThresh(bmask) = PPMVB(bmask) .* (PPMVB(bmask)>PPMpThresh);
        
        % Save nii
        nii = load_untouch_nii(strcat(resultsPath,'/Cbeta_',num,'.nii'));
        nii.img = PPMVB;
        save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_PPM_',num,'.nii'));
        nii.img = PPMVBThresh;
        save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_PPMThresh_',num,'.nii'));

    end
    
    % Compute PPMs for contrasts
    for k = 1:nContr
        num = strcat('000',num2str(k));
        num = num(end-3:end);
        
        % Load
        vol = spm_vol(strcat(resultsPath,'/con_',num,'.nii'));
        [wVB,XYZ] = spm_read_vols(vol);
        
        % Load stds
        vol = spm_vol(strcat(resultsPath,'/con_sd_',num,'.nii'));
        [stdVB,XYZ] = spm_read_vols(vol);
        
        % Compute marginal PPM
        PPMVB = nan(size(wVB));
        PPMVBThresh = nan(size(wVB));
        bmask = ~isnan(wVB);
        PPMVB(bmask) = normcdf((wVB(bmask)-PPMwThresh) ./ stdVB(bmask));
        PPMVBThresh(bmask) = PPMVB(bmask) .* (PPMVB(bmask)>PPMpThresh);
        
        % Save nii
        nii = load_untouch_nii(strcat(resultsPath,'/con_',num,'.nii'));
        nii.img = PPMVB;
        save_untouch_nii(nii,strcat(resultsPath,'/con_PPM_',num,'.nii'));
        nii.img = PPMVBThresh;
        save_untouch_nii(nii,strcat(resultsPath,'/con_PPMThresh_',num,'.nii'));
        
    end

elseif findstr(method,'MCMC')
    
    load(strcat(resultsPath,'/MCMC.mat'));
    VBResultsPath = strcat(outputPath,'sub',subjStr,'/',VBMethod);
    load(strcat(VBResultsPath,'/SPM.mat'));
    PPMwThresh = 100*PPMwPercThresh/(max(SPM.xBF.bf(:,1))/SPM.xBF.dt);
    K = size(SPM.xX.X,2);
    nContr = length(SPM.xCon);
    nLb = length(MCMC.a.sliceNbrs);
    
    % Load template from VB results and later overwrite this to save
    vol = spm_vol(strcat(VBResultsPath,'/Cbeta_0001.nii'));
    [wMCMC,XYZ] = spm_read_vols(vol);
    [stdMCMC,XYZ] = spm_read_vols(vol);
    [PPMMCMC,XYZ] = spm_read_vols(vol);
    [PPMMCMCThresh,XYZ] = spm_read_vols(vol);
    
    % Save means and stds as nifti
    for k = 1:K
        
        num = strcat('000',num2str(k));
        num = num(end-3:end);
        
        if strcmp(method,'MCMC2D')
            
            for j = 1:nLb

                % mean
                mask2d = ~isnan(wMCMC(:,:,MCMC.a.sliceNbrs(j)));
                tempSlice = nan(size(mask2d));
                tempSlice(mask2d) = MCMC.b(j).wPostMean(k,:);
                wMCMC(:,:,MCMC.a.sliceNbrs(j)) = tempSlice;
                
                % std
                stdSlice = nan(size(mask2d));
                tempSlice(mask2d) = MCMC.b(j).wPostStd(k,:);
                stdMCMC(:,:,MCMC.a.sliceNbrs(j)) = tempSlice;
                
            end
        
        elseif strcmp(method,'MCMC3D')
            
            % mean
            wMCMC(~isnan(wMCMC)) = MCMC.b.wPostMean(k,:);

            % std
            stdMCMC(~isnan(stdMCMC)) = MCMC.b.wPostStd(k,:);
            
        end
        
        nii = load_untouch_nii(strcat(VBResultsPath,'/Cbeta_',num,'.nii'));
        nii.img = wMCMC;
        save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_',num,'.nii'));
        
        nii = load_untouch_nii(strcat(VBResultsPath,'/SDbeta_',num,'.nii'));
        nii.img = stdMCMC;
        save_untouch_nii(nii,strcat(resultsPath,'/SDbeta_',num,'.nii'));
        
    end
    
    % Compute PPMs for relevant regressors (exclude temporal derivative
    % hrf, hmp and intercept regressors)
    for k = 1:2:(K-7)
        num = strcat('000',num2str(k));
        num = num(end-3:end);
        
        if strcmp(method,'MCMC2D')
            
            for j = 1:nLb
                
                % PPM
                mask2d = ~isnan(PPMMCMC(:,:,MCMC.a.sliceNbrs(j)));
                wVec2Temp = squeeze(MCMC.b(j).wVec2(k,:,:));
                tempSlice = nan(size(mask2d));
                tempSlice(mask2d) = mean(wVec2Temp > PPMwThresh,2);
                PPMMCMC(:,:,MCMC.a.sliceNbrs(j)) = tempSlice;
                
            end
        
        elseif strcmp(method,'MCMC3D')
            
            % PPM
            wVec2Temp = squeeze(MCMC.b.wVec2(k,:,:));
            PPMMCMC(~isnan(PPMMCMC)) = mean(wVec2Temp > PPMwThresh,2);
            
        end
        
        % Thresholded PPM
        bmask = ~isnan(PPMMCMC);
        PPMMCMCThresh(bmask) = PPMMCMC(bmask) .* (PPMMCMC(bmask)>PPMpThresh);
           
        % Save nii
        nii = load_untouch_nii(strcat(resultsPath,'/Cbeta_',num,'.nii'));
        nii.img = PPMMCMC;
        save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_PPM_',num,'.nii'));
        nii.img = PPMMCMCThresh;
        save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_PPMThresh_',num,'.nii'));

    end
    
    % Compute mean, stds and PPMs for contrasts
    for k = 1:nContr
        num = strcat('000',num2str(k));
        num = num(end-3:end);
        
        contrast = SPM.xCon(k).c;
        
        if strcmp(method,'MCMC2D')
            
            for j = 1:nLb
                
                mask2d = ~isnan(PPMMCMC(:,:,MCMC.a.sliceNbrs(j)));
                wVec2Temp = zeros(size(squeeze(MCMC.b(j).wVec2(1,:,:))));
                for kk = 1:K
                	wVec2Temp = wVec2Temp + contrast(kk) * squeeze(MCMC.b(j).wVec2(kk,:,:));
                end
                
                % Contrast mean
                tempSlice = nan(size(mask2d));
                tempSlice(mask2d) = mean(wVec2Temp,2);
                wMCMC(:,:,MCMC.a.sliceNbrs(j)) = tempSlice;
                
                % Contrast std
                stdTemp = zeros(size(wVec2Temp,1),1);
                for i = 1:length(stdTemp);
                    stdTemp(i) = std(wVec2Temp(i,:));
                end
                tempSlice = nan(size(mask2d));
                tempSlice(mask2d) = stdTemp;
                stdMCMC(:,:,MCMC.a.sliceNbrs(j)) = tempSlice;
                
                % Constrast PPM
                tempSlice = nan(size(mask2d));
                tempSlice(mask2d) = mean(wVec2Temp > PPMwThresh,2);
                PPMMCMC(:,:,MCMC.a.sliceNbrs(j)) = tempSlice;
                
            end
        
        elseif strcmp(method,'MCMC3D')
            
            wVec2Temp = zeros(size(squeeze(MCMC.b.wVec2(1,:,:))));
            for kk = 1:K
                wVec2Temp = wVec2Temp + contrast(kk) * squeeze(MCMC.b.wVec2(kk,:,:));
            end
            
            % Contrast mean
            wMCMC(~isnan(wMCMC)) = mean(wVec2Temp,2);

            % Contrast std
            stdTemp = zeros(size(wVec2Temp,1),1);
            for i = 1:length(stdTemp);
                stdTemp(i) = std(wVec2Temp(i,:));
            end
            stdMCMC(~isnan(wMCMC)) = stdTemp;
            
            % Contrast PPM
            PPMMCMC(~isnan(PPMMCMC)) = mean(wVec2Temp > PPMwThresh,2);
            
        end
        
        % Thresholded PPM
        bmask = ~isnan(PPMMCMC);
        PPMMCMCThresh(bmask) = PPMMCMC(bmask) .* (PPMMCMC(bmask)>PPMpThresh);
           
        % Save nii
        nii = load_untouch_nii(strcat(resultsPath,'/Cbeta_',num,'.nii'));
        nii.img = wMCMC;
        save_untouch_nii(nii,strcat(resultsPath,'/con_',num,'.nii'));
        nii.img = stdMCMC;
        save_untouch_nii(nii,strcat(resultsPath,'/con_sd_',num,'.nii'));
        nii.img = PPMMCMC;
        save_untouch_nii(nii,strcat(resultsPath,'/con_PPM_',num,'.nii'));
        nii.img = PPMMCMCThresh;
        save_untouch_nii(nii,strcat(resultsPath,'/con_PPMThresh_',num,'.nii'));
        
    end    
    
end