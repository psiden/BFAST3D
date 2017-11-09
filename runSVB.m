%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run SVB estimation for BIDS format fMRI data
%
% INPUT:        dS - data settings struct
%               VBMethod - string containing estimation method,
%                                  possible options: SPMsVB2d, SPMsVB3d,
%                                  ImprovedVB2d, ImprovedVB3d.
%               
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
% REVISED:      2017-10-26
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runSVB(dS,VBMethod)
%% Run

% Read settings
dataPath = dS.dataPath; outputPath = dS.outputPath; SPMPath = dS.SPMPath; 
subjStr = dS.subjStr; fileStr = dS.fileStr; TR = dS.TR;
condStr = dS.condStr; contrastStr = dS.contrastStr;

% SPM Paths
SPMPath2 = strcat(SPMPath,'/svb');
SPMPath3 = strcat(SPMPath2,'/',VBMethod);
addpath(SPMPath);
addpath(SPMPath2);
addpath(SPMPath3);

% BOLD path
BOLDDataPath = [dataPath,'sub-',subjStr,'/func/'];

% Results path
resultsPath = [outputPath,'sub-',subjStr,'/',VBMethod];
mkdir(resultsPath);

% Delete old SPM file
SPMFile = strcat(resultsPath,'/SPM.mat'); 
if exist(SPMFile,'file')==2 
  %system(['rm' SPMFile]); % Linux 
  delete(SPMFile) 
end

% Initialise SPM defaults
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Number of fMRI volumes
load([BOLDDataPath,fileStr,'bold.mat'],'mat');
TRs = size(mat,3); 
clear mat;

% BOLD data
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(resultsPath);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
scans = cell(TRs,1);
for t = 1:TRs
  scans{t} = [BOLDDataPath,'r',fileStr,'bold.nii,',num2str(t)];
end
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;

% Conditions
condCell = readTaskFile([dataPath,condStr]); 
eventCell = readTaskFile([BOLDDataPath,fileStr,'events.tsv']); 
nCond = size(condCell,1);
nEvents = size(eventCell,1);
for k = 1:nCond
  onsets = []; durations = [];
  for j = 1:nEvents
    if strcmp(strcat(condCell{k,3}),strcat(eventCell{j,3}))
      onsets = [onsets;str2num(eventCell{j,1})]; 
      durations = [durations;str2num(eventCell{j,2})]; 
    end
  end  
  matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).name = condCell{k,3};
  matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).onset = onsets;
  matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).duration = durations;
end

% Contrasts
fid = fopen([dataPath,contrastStr]);
commandString = ['%s%s']; % Always start by reading two strings
for k = 1:nCond % Read one float for each regressor
    commandString = [commandString '%f'];
end 
text = textscan(fid,commandString);  
fclose(fid);
nContrasts = size(text{3},1);
contrastNames = text{2};
contrastVectors = zeros(nContrasts,2*nCond);
for k = 1:nContrasts
    for l = 1:nCond
        scalars = text{2+l};
        contrastVectors(k,2*l-1) = scalars(k);
    end
    contrastVectors(k,:) = contrastVectors(k,:) / sum(abs(contrastVectors(k,:)));
end

% Head motion parameters
hmpFile = [BOLDDataPath,'rp_',fileStr,'bold.txt'];
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(hmpFile);

% High-pass filter
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

% HRF with one temporal derivative
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1,0];

% Estimate model
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(SPMFile);    
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.LogEv = 'No';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.anova.first = 'No';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.anova.second = 'Yes';
for k = 1:nContrasts
    contrastName = contrastNames{k};
    contrastVector = contrastVectors(k,:);
    matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.gcon(k).name = contrastName; 
    matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.gcon(k).convec = contrastVector;
end
spm_jobman('run',matlabbatch);

% Remove SPM path
rmpath(SPMPath);
