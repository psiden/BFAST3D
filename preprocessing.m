%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Unpack and realignment for OpenfMRI datasets
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
function preprocessing(dataPath,SPMPath,subject)
%% Run

% Add SPM path
addpath(SPMPath);

% BOLD path
subjStr = strcat('00',num2str(subject));
subjStr = subjStr(end-2:end);
BOLDDataPath = strcat(dataPath,'sub',subjStr,'/BOLD/task001_run001/');

% Unpack nifti file
gunzip(strcat(BOLDDataPath,'bold.nii.gz'));

% Initialise SPM defaults
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Realign
clear matlabbatch
f = strcat(BOLDDataPath,'bold.nii');
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cellstr(f);
spm_jobman('run',matlabbatch);

% Remove SPM path
rmpath(SPMPath);
