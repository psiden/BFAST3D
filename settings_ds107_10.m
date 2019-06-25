%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Settings for OpenfMRI ds107
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2017-11-24
% REVISED:      2019-06-25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

dS.dataPath = '.../ds000107_R2.0.2/';
dS.outputPath = '.../Results/';
dS.SPMPath = '.../spm12/';

dS.subjStr = '10';
dS.taskStr = 'onebacktask';
dS.runStr = '_run-01_';
dS.fileStr = ['sub-',dS.subjStr,'_task-',dS.taskStr,dS.runStr];
dS.TR = 3;
dS.condStr = 'condition_key.txt';
dS.contrastStr = 'task_contrasts.txt';