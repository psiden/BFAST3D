%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run SVB estimation
%
% INPUT:        dataPath - string containing openfMRI data folder
%               outputPath - string containing output folder
%               SPMPath - string containing SPM path
%               subject - subject number (int)
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
% REVISED:      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runSVB(dataPath,outputPath,SPMPath,subject,VBMethod)
%% Run

% SPM Paths
SPMPath2 = strcat(SPMPath,'/svb');
SPMPath3 = strcat(SPMPath2,'/',VBMethod);
addpath(SPMPath);
addpath(SPMPath2);
addpath(SPMPath3);

% BOLD path
subjStr = strcat('00',num2str(subject));
subjStr = subjStr(end-2:end);
BOLDDataPath = strcat(dataPath,'sub',subjStr,'/BOLD/task001_run001/');

% Results path
resultsPath = strcat(outputPath,'sub',subjStr,'/',VBMethod);
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

% Repetition time
TRCell = textscan(fopen(strcat(dataPath,'scan_key.txt')),'%s%f'); 
TR = TRCell{2};

% Number of fMRI volumes
load(strcat(BOLDDataPath,'bold.mat'),'mat');
TRs = size(mat,3); 
clear mat;

% BOLD data
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(resultsPath);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
scans = cell(TRs,1);
for t = 1:TRs
    scans{t} = strcat(BOLDDataPath,'rbold.nii,',num2str(t));
end
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;

% Conditions
condPath = strcat(dataPath,'sub',subjStr,'/model/model001/onsets/task001_run001');
condCell = readTaskFile(strcat(strcat(dataPath,'models/model001/condition_key.txt'))); 
nCond = size(condCell,1);
for k = 1:nCond
    kStr = strcat('00',num2str(k));
    kStr = kStr(end-2:end);
    fid = fopen(strcat(condPath,'/cond',kStr,'.txt'));
    text = textscan(fid,'%f%f%f'); % Read 3 floats
    fclose(fid);

    onsets = text{1}; 
    durations = text{2};
    values = text{3};

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).name = condCell{k,3};
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).onset = onsets;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).duration = durations;
end

% Contrasts
fid = fopen(strcat(dataPath,'/models/model001/task_contrasts.txt'));
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
hmpFile = strcat(BOLDDataPath,'rp_bold.txt');
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
