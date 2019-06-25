%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run estimation files for SVB/MCMC for BIDS format fMRI data
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
% REVISED:      2019-06-25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify
clear all
close all
clc

%% Setup data settings (dS)
settings_ds107_10;
subject = str2num(dS.subjStr);

%% Download data
download_ds107_data(dS);

%% Preprocessing

% Head motion correction
preprocessing(dS);

%% Run 3D

% Run Dummy VB (Dummy method to convert data into SPM.mat format)
VBMethod = 'DVB3D';
runSVB(dS,VBMethod);

% Run EB Matern
fS.boldPath = dS.boldPath; 
fS.dataPath = [dS.outputPath,'sub-',dS.subjStr,'/'];
fS.outputPath = fS.dataPath;
fS.SPMResultsFolder = VBMethod;
fS.ResultsFolder = 'M2Iso3dEyeFix';
fS.doPlot = 0;
fS.doParallel = 0;
fS.maxiter = 200;
addpath(genpath('eb'));
runEBMatern(fS);

% Run SVB
VBMethod = 'SVB3D';
runSVB(dS,VBMethod);

% Run MCMC
MCMCMethod = 'MCMC3D';
samplingMethod = 'PCG';
runMCMC(dS,MCMCMethod,VBMethod,samplingMethod);

% Compute marginal and joint PPMs
addpath(dS.SPMPath);

% Requires Tools for NIfTI and ANALYZE image Matlab-package
computePPMs(dS.outputPath,subject,'EBMatern',fS);
computePPMs(dS.outputPath,subject,VBMethod);
computePPMs(dS.outputPath,subject,MCMCMethod);

% Requires excursions R-package
contrastNbr = 1;
computeExcursions(dS.outputPath,subject,MCMCMethod,contrastNbr,1)

% %% Run 2D
% 
% % Run SVB
% VBMethod = 'SVB2D';
% runSVB(dS,VBMethod);
% 
% % Run MCMC
% MCMCMethod = 'MCMC2D';
% samplingMethod = 'PCG';
% runMCMC(dS,MCMCMethod,VBMethod,samplingMethod);
% 
% % Compute marginal and joint PPMs
% addpath(dS.SPMPath);
% 
% % Requires Tools for NIfTI and ANALYZE image Matlab-package
% computePPMs(dS.outputPath,subject,VBMethod);
% computePPMs(dS.outputPath,subject,MCMCMethod);
% 
% % Requires excursions R-package
% contrastNbr = 1;
% sliceNbr = 7;
% computeExcursions(dS.outputPath,subject,MCMCMethod,contrastNbr,sliceNbr)

rmpath(dS.SPMPath);
