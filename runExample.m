%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run estimation files for SVB/MCMC for BIDS format fMRI data
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
% REVISED:      2017-10-27
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify
clear all
close all
clc

%% Setup data settings (dS)
settings_ds105_01;
subject = str2num(dS.subjStr);

%% Preprocessing

% Head motion correction
preprocessing(dS);

%% Run 3D

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
computePPMs(dS.outputPath,subject,VBMethod);
computePPMs(dS.outputPath,subject,MCMCMethod);

% Requires excursions R-package
contrastNbr = 5;
computeExcursions(dS.outputPath,subject,MCMCMethod,contrastNbr,1)

rmpath(dS.SPMPath);

%% Run 2D

% Run SVB
VBMethod = 'SVB2D';
runSVB(dS,VBMethod);

% Run MCMC
MCMCMethod = 'MCMC2D';
samplingMethod = 'PCG';
runMCMC(dS,MCMCMethod,VBMethod,samplingMethod);

% Compute marginal and joint PPMs
addpath(dS.SPMPath);

% Requires Tools for NIfTI and ANALYZE image Matlab-package
computePPMs(dS.outputPath,subject,VBMethod);
computePPMs(dS.outputPath,subject,MCMCMethod);

% Requires excursions R-package
contrastNbr = 5;
sliceNbr = 30;
computeExcursions(dS.outputPath,subject,MCMCMethod,contrastNbr,sliceNbr)

rmpath(dS.SPMPath);