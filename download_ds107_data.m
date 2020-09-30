%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Download ds107 data from openneuro.org. Task and condition
%               data are pre-stored in the BFAST3D repo since openneuro
%               seems to have removed them.
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2019-02-18
% REVISED:      2019-06-25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function download_ds107_data(dS)

% Read settings
dataPath = dS.dataPath;
subjStr = dS.subjStr;
fileStr = dS.fileStr;

% BOLD path
BOLDDataPath = [dataPath,'sub-',subjStr,'/func/'];

% Download data
disp('Downloading word and object dataset.');
mkdir(BOLDDataPath);
urlwrite(['https://openneuro.org/crn/datasets/ds000107/snapshots/58054c33cce88d000ca33620',...
          '/files/sub-10:func:sub-10_task-onebacktask_run-01_bold.nii.gz'],...
          [BOLDDataPath,dS.fileStr,'bold.nii.gz']);
urlwrite(['https://openneuro.org/crn/datasets/ds000107/snapshots/58054c33cce88d000ca33620',...
          '/files/sub-10:func:sub-10_task-onebacktask_run-01_events.tsv'],...
          [BOLDDataPath,dS.fileStr,'events.tsv']);
urlwrite(['https://github.com/psiden/BFAST3D/blob/master/data/ds000107_R2.0.2/condition_key.txt'],...
          [dataPath,'condition_key.txt']);
urlwrite(['https://github.com/psiden/BFAST3D/blob/master/data/ds000107_R2.0.2/task_contrasts.txt'],...
          [dataPath,'task_contrasts.txt']);

% Unpack nifti file
gunzip([BOLDDataPath,fileStr,'bold.nii.gz']);