%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Settings for OpenfMRI ds105
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2017-10-26
% REVISED:      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

dS.dataPath = '.../ds000105_R2.0.2/';
dS.outputPath = '.../Results/';
dS.SPMPath = '.../spm12/';

dS.subjStr = '1';
dS.taskStr = 'objectviewing';
dS.runStr = '_run-01_';
dS.fileStr = ['sub-',dS.subjStr,'_task-',dS.taskStr,dS.runStr];
dS.TR = 2.5;
dS.condStr = 'condition_key.txt';
dS.contrastStr = 'task_contrasts.txt';