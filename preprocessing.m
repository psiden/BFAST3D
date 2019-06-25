%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Motion correction for openneuro datasets
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
% REVISED:      2019-02-18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preprocessing(dS)
%% Run

% Read settings
dataPath = dS.dataPath; SPMPath = dS.SPMPath; subjStr = dS.subjStr;
fileStr = dS.fileStr;

% Add SPM path
addpath(SPMPath);

% BOLD path
BOLDDataPath = [dataPath,'sub-',subjStr,'/func/'];

% Unpack nifti file
gunzip([BOLDDataPath,fileStr,'bold.nii.gz']);

% Initialise SPM defaults
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Realign
clear matlabbatch
f = [BOLDDataPath,fileStr,'bold.nii'];
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cellstr(f);
spm_jobman('run',matlabbatch);

% Remove SPM path
rmpath(SPMPath);
