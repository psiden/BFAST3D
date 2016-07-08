%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run estimation files for SVB/MCMC
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
%% Specify
clear all
close all
clc

%% Setup

dataPath = '.../ds105/';
outputPath = '.../Results/';
SPMPath = '.../spm12/';
subject = 1;

%% Preprocessing

% Head motion correction
preprocessing(dataPath,SPMPath,subject)

%% Run 3D

% Run SVB
VBMethod = 'SVB3D';
runSVB(dataPath,outputPath,SPMPath,subject,VBMethod);

% Run MCMC
MCMCMethod = 'MCMC3D';
samplingMethod = 'PCG';
runMCMC(outputPath,subject,MCMCMethod,VBMethod,samplingMethod);

% Compute marginal and joint PPMs
addpath(SPMPath);

% Requires Tools for NIfTI and ANALYZE image Matlab-package
computePPMs(outputPath,subject,VBMethod);
computePPMs(outputPath,subject,MCMCMethod);

% Requires excursions R-package
contrastNbr = 5;
computeExcursions(outputPath,subject,MCMCMethod,contrastNbr,1)

rmpath(SPMPath);

%% Run 2D

% Run SVB
VBMethod = 'SVB2D';
runSVB(dataPath,outputPath,SPMPath,subject,VBMethod);

% Run MCMC
MCMCMethod = 'MCMC2D';
samplingMethod = 'PCG';
runMCMC(outputPath,subject,MCMCMethod,VBMethod,samplingMethod);

% Compute marginal and joint PPMs
addpath(SPMPath);

% Requires Tools for NIfTI and ANALYZE image Matlab-package
computePPMs(outputPath,subject,VBMethod);
computePPMs(outputPath,subject,MCMCMethod);

% Requires excursions R-package
contrastNbr = 5;
sliceNbr = 30;
computeExcursions(outputPath,subject,MCMCMethod,contrastNbr,sliceNbr)

rmpath(SPMPath);
