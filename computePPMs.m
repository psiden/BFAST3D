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
function computePPMs(outputPath,subject,method,fS)
%% Setup

% Needed for MCMC for some additional input, like contrasts
if findstr(method,'MCMC')
  VBMethod = strcat('SVB',method(5:6));
end

PPMpThresh = 0.90;
PPMwPercThresh = .005;

%% Run

% subjStr = strcat('00',num2str(subject));
% subjStr = subjStr(end-2:end);
subjStr = ['-',num2str(subject)];

if findstr(method,'EBMatern')
  resultsPath = strcat(outputPath,'sub',subjStr,'/',fS.ResultsFolder);
else
  resultsPath = strcat(outputPath,'sub',subjStr,'/',method);
end

% Main if-statement
if findstr(method,'EBMatern')
  
  load(strcat(resultsPath,'/Output.mat'));
  sliceNbrs = Output.sliceNbrs;
  SPMMatPath = [fS.dataPath,fS.SPMResultsFolder,'/'];
  load([SPMMatPath,'SPM.mat']);
  model = fS.ResultsFolder;
  doSlices = 1:500;
%   
%   nu = 2 - 3/2; ndim = 3;
%   if strfind(model,'M2Iso')
%     range = sqrt(8*nu) ./ sqrt(Output.b.kappa2)
%     margstd = sqrt(gamma(nu) ./ (gamma(nu+ndim/2)*((4*pi)^(ndim/2)).*Output.b.kappa2.^nu.*...
%       Output.b.tau2))
%   end
  
  PPMwThresh = 100*PPMwPercThresh/(max(SPM.xBF.bf(:,1))/SPM.xBF.dt);
  K = Output.fS.K; % size(SPM.xX.X,2); %
  P = SPM.PPM.AR_P;
  if strfind(resultsPath,'Simu')
    nContr = 0;
  else
    nContr = length(SPM.xCon);
  end
  nLb = length(Output.b);
  if strfind(model,'2d'); nLb = max(length(Output.b),length(sliceNbrs));end
  sz = SPM.xY.VY(1).dim;
  
  % Load template from SPM results and later overwrite this to save
  vol = spm_vol(strcat(SPMMatPath,'/Cbeta_0001.nii'));
  [w,XYZ] = spm_read_vols(vol);
  [stdw,XYZ] = spm_read_vols(vol);
  [PPM,XYZ] = spm_read_vols(vol);
  [PPMThresh,XYZ] = spm_read_vols(vol);
  [SDerror,XYZ] = spm_read_vols(vol);
  [AR,XYZ] = spm_read_vols(vol);
  
  % Save means and stds as nifti
  for k = 1:K
    
    num = strcat('000',num2str(k));
    num = num(end-3:end);
    
    if strfind(model,'2d')
      
      for j = 1:nLb
        
        % mean
        mask2d = ~isnan(w(:,:,sliceNbrs(j)));
        tempSlice = nan(size(mask2d));
        
        if find(doSlices==j)
          wTemp = Output.b(j).w(k,:)';
          tempSlice(mask2d) = wTemp;
          N = size(Output.b(j).w,2);
        end
        w(:,:,sliceNbrs(j)) = tempSlice;
        
        
        % std
        if find(doSlices==j)
          stdSlice = nan(size(mask2d));
          stdtmp = diag(Output.b(j).wCovMat);
          stdTemp = sqrt(stdtmp(k:K:N*K));
          tempSlice(mask2d) = stdTemp;
        end
        stdw(:,:,sliceNbrs(j)) = tempSlice;
        
        
        % PPM, PPMThresh
        tempSlice = nan(size(mask2d));
        if find(doSlices==j)
          tempSlice(mask2d) = normcdf((wTemp-PPMwThresh)./stdTemp);
        end
        PPM(:,:,sliceNbrs(j)) = tempSlice;
        PPMThresh(:,:,sliceNbrs(j)) = tempSlice.*(tempSlice>PPMpThresh);
        
      end
      
    elseif strfind(model,'3d')
      
      for j = 1:nLb
        
        % mean
        if strfind(resultsPath,'Simu')
          maskInd = Output.b(j).maskInd;
        else
          maskParc = SPM.xVol.XYZ(:,find(SPM.xVol.labels==j));
          maskInd = sub2ind(sz,maskParc(1,:)',maskParc(2,:)',maskParc(3,:)');
        end
        
        
        if find(doSlices==j)
          wTemp = Output.b(j).w(k,:)';
          N = size(Output.b(j).w,2);
        else
          wTemp = nan(size(maskInd));
        end
        w(:,:,:) = nan;
        w(maskInd) = wTemp;
        
        % std
        if find(doSlices==j)
          stdtmp = diag(Output.b(j).wCovMat);
          stdTemp = sqrt(stdtmp(k:K:N*K));
        else
          stdTemp = nan(size(maskInd));
        end
        stdw(:,:,:) = nan;
        stdw(maskInd) = stdTemp;
        
        % PPM, PPMThresh
        if find(doSlices==j)
          PPMTemp = normcdf((wTemp-PPMwThresh)./stdTemp);
        else
          PPMTemp = nan(size(maskInd));
        end
        PPM(:,:,:) = nan;
        PPM(maskInd) = PPMTemp;
        PPMThresh(:,:,:) = nan;
        PPMThresh(maskInd) = PPMTemp.*(PPMTemp>PPMpThresh);
        
      end
    end
    
    nii = load_untouch_nii(strcat(SPMMatPath,'/Cbeta_',num,'.nii'));
    nii.img = w;
    save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_',num,'.nii'));
    
    nii = load_untouch_nii(strcat(SPMMatPath,'/SDbeta_',num,'.nii'));
    nii.img = stdw;
    save_untouch_nii(nii,strcat(resultsPath,'/SDbeta_',num,'.nii'));
    
    % Save nii
    nii = load_untouch_nii(strcat(resultsPath,'/Cbeta_',num,'.nii'));
    nii.img = PPM;
    save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_PPM_',num,'.nii'));
    nii.img = PPMThresh;
    save_untouch_nii(nii,strcat(resultsPath,'/Cbeta_PPMThresh_',num,'.nii'));
    
  end
  
  % Compute mean, stds and PPMs for contrasts
  for k = 1:nContr
    num = strcat('000',num2str(k));
    num = num(end-3:end);
    
    contrast = SPM.xCon(k).c;
    
    if strfind(model,'2d')
      
      for j = 1:nLb
        
        % mean
        mask2d = ~isnan(PPM(:,:,sliceNbrs(j)));
        tempSlice = nan(size(mask2d));
                
        if find(doSlices==j)
          N = size(Output.b(j).w,2);
          iQData = reshape(repmat(1:K,N*K,1)',N*K*K,1) + ...
            reshape(repmat(K*(0:N-1),K*K,1),N*K*K,1);
          jQData = reshape(repmat(1:N*K,K,1),N*K*K,1);
          indQData = sub2ind([N*K,N*K],iQData,jQData);
          wMean = Output.b(j).w' * contrast;
        end
        
        % Contrast mean
        tempSlice = nan(size(mask2d));
        if find(doSlices==j)
          tempSlice(mask2d) = wMean;
        end
        w(:,:,sliceNbrs(j)) = tempSlice;
        
        % Contrast std
        tempSlice = nan(size(mask2d));
        if find(doSlices==j)
          wCov = reshape(full(Output.b(j).wCovMat(indQData)),[K,K,N]);
          covTemp = zeros(K,N);
          for kk = 1:K
            covTemp(kk,:) = dot(squeeze(wCov(kk,:,:)),repmat(contrast,[1,N]));
          end
          stdTemp = sqrt(contrast' * covTemp)';
          tempSlice(mask2d) = stdTemp;
        end
        stdw(:,:,sliceNbrs(j)) = tempSlice;
        
        % Constrast PPM, PPMThresh
        tempSlice = nan(size(mask2d));
        if find(doSlices==j)
          tempSlice(mask2d) = normcdf((wMean-PPMwThresh)./stdTemp);
        end
        PPM(:,:,sliceNbrs(j)) = tempSlice;
        PPMThresh(:,:,sliceNbrs(j)) = tempSlice.*(tempSlice>PPMpThresh);
      end
      
    elseif strfind(model,'3d')
      
      for j = 1:nLb
        
        % mean
        if strfind(resultsPath,'Simu')
          maskInd = Output.b(j).maskInd;
        else
          maskParc = SPM.xVol.XYZ(:,find(SPM.xVol.labels==j));
          maskInd = sub2ind(sz,maskParc(1,:)',maskParc(2,:)',maskParc(3,:)');
        end
        
        if find(doSlices==j)
          N = size(Output.b(j).w,2);
          iQData = reshape(repmat(1:K,N*K,1)',N*K*K,1) + ...
            reshape(repmat(K*(0:N-1),K*K,1),N*K*K,1);
          jQData = reshape(repmat(1:N*K,K,1),N*K*K,1);
          indQData = sub2ind([N*K,N*K],iQData,jQData);
          wMean = Output.b(j).w' * contrast;
        else
          wMean = nan(size(maskInd));
        end
        w(:,:,:) = nan;
        w(maskInd) = wMean;
        
        % Contrast std
        if find(doSlices==j)
          wCov = reshape(full(Output.b(j).wCovMat(indQData)),[K,K,N]);
          covTemp = zeros(K,N);
          for kk = 1:K
            covTemp(kk,:) = dot(squeeze(wCov(kk,:,:)),repmat(contrast,[1,N]));
          end
          stdTemp = sqrt(contrast' * covTemp)';
        else
          stdTemp = nan(size(maskInd));
        end
        stdw(:,:,:) = nan;
        stdw(maskInd) = stdTemp;
        
        % Constrast PPM, PPMThresh
        if find(doSlices==j)
          PPMTemp = normcdf((wMean-PPMwThresh)./stdTemp);
        else
          PPMTemp = nan(size(maskInd));
        end
        PPM(:,:,:) = nan;
        PPM(maskInd) = PPMTemp;
        PPMThresh(:,:,:) = nan;
        PPMThresh(maskInd) = PPMTemp.*(PPMTemp>PPMpThresh);
      end
      
    end
    
    % Save nii
    nii = load_untouch_nii(strcat(resultsPath,'/Cbeta_',num,'.nii'));
    nii.img = w;
    save_untouch_nii(nii,strcat(resultsPath,'/con_',num,'.nii'));
    nii.img = stdw;
    save_untouch_nii(nii,strcat(resultsPath,'/con_sd_',num,'.nii'));
    nii.img = PPM;
    save_untouch_nii(nii,strcat(resultsPath,'/con_PPM_',num,'.nii'));
    nii.img = PPMThresh;
    save_untouch_nii(nii,strcat(resultsPath,'/con_PPMThresh_',num,'.nii'));
    
  end
  
  % Save measurement noise variances
  if strfind(model,'2d')
    
    for j = 1:nLb
      mask2d = ~isnan(SDerror(:,:,sliceNbrs(j)));
      tempSlice = nan(size(mask2d));
      if find(doSlices==j)
        sdTemp = 1./sqrt(Output.b(j).lambda);
        tempSlice(mask2d) = sdTemp;
        %       N = size(Output.b(j).w,2);
      end
      SDerror(:,:,sliceNbrs(j)) = tempSlice;
    end
    
  elseif strfind(model,'3d')
    
    %   sdTemp = 1./sqrt(Output.b.lambda);
    %   SDerror(~isnan(SDerror)) = sdTemp;
    
    for j = 1:nLb
      if strfind(resultsPath,'Simu')
        maskInd = Output.b(j).maskInd;
      else
        maskParc = SPM.xVol.XYZ(:,find(SPM.xVol.labels==j));
        maskInd = sub2ind(sz,maskParc(1,:)',maskParc(2,:)',maskParc(3,:)');
      end
      if find(doSlices==j)
        sdTemp = 1./sqrt(Output.b(j).lambda);
      else
        sdTemp = nan(size(maskInd));
      end
      SDerror(:,:,:) = nan;
      SDerror(maskInd) = sdTemp;
    end
    
  end
  
  nii = load_untouch_nii(strcat(SPMMatPath,'/Sess1_SDerror.nii'));
  nii.img = SDerror;
  save_untouch_nii(nii,strcat(resultsPath,'/Sess1_SDerror.nii'));
  
  % AR parameters
  for p = 1:P
    
    num = strcat('000',num2str(p));
    num = num(end-3:end);
    
    if strfind(model,'2d')
      
      for j = 1:nLb
        
        % AR
        mask2d = ~isnan(AR(:,:,sliceNbrs(j)));
        tempSlice = nan(size(mask2d));
        
        if find(doSlices==j)
          ARTemp = Output.b(j).a(p,:)';
          tempSlice(mask2d) = ARTemp;
        end
        AR(:,:,sliceNbrs(j)) = tempSlice;
      end
      
    elseif strfind(model,'3d')
      
      %     % AR
      %     ARTemp = Output.b.a(p,:)';
      %     AR(~isnan(AR)) = ARTemp;
      
      for j = 1:nLb
        if strfind(resultsPath,'Simu')
          maskInd = Output.b(j).maskInd;
        else
          maskParc = SPM.xVol.XYZ(:,find(SPM.xVol.labels==j));
          maskInd = sub2ind(sz,maskParc(1,:)',maskParc(2,:)',maskParc(3,:)');
        end
        if find(doSlices==j)
          ARTemp = Output.b(j).a(p,:)';
        else
          ARTemp = nan(size(maskInd));
        end
        AR(:,:,:) = nan;
        AR(maskInd) = ARTemp;
      end
    end
    
    nii = load_untouch_nii(strcat(SPMMatPath,'/Sess1_AR_',num,'.nii'));
    nii.img = AR;
    save_untouch_nii(nii,strcat(resultsPath,'/Sess1_AR_',num,'.nii'));
  end
  
  
  
elseif findstr(method,'VB')
  
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