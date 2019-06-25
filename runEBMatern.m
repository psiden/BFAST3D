%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run Empirical Bayes estimation of the hyperparameters in
%               the model with Matern spatial priors
%
%               Based on the Paper "Matern 3D fMRI" (2019).
%
% INPUT:        fS - fMRI analysis settings struct
%               
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2018-06-19
% REVISED:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runEBMatern(fS)

% Check if the Results folder is defined and sets it otherwise
if ~isfield(fS,'ResultsFolder'); fS.ResultsFolder = 'M2Iso3dEyeFix'; end
% Whether to run multiple threads and to show convergence plots
if ~isfield(fS,'doParallel'); fS.doParallel = 1; end 
if ~isfield(fS,'doPlot'); fS.doPlot = 0; end
% Similar checks could be applied to the other settings, to make them
% chooseable outside this function, but this is left to the user

% Spatial prior settings: {Prior,startvalues,fixed?}
fS.hrfPrior = {'Matern2Iso',[.01,.1],[0,0]}; % HRF and HRF derivative
fS.hmpPrior = {'Eye',1e-12,1}; % Head motion parameters
fS.interceptPrior = {'Eye',1e-12,1}; % Intercept

% Random seed to use
fS.rngseed = 102; if fS.rngseed > 0; rng(fS.rngseed); end

% Plot and print settings
fS.plotUpdateFreq = 10;fS.plotLength = 1000;
fS.doPrint = 1;

% Optimization settings
if ~isfield(fS,'maxiter'); fS.maxiter = 200; end % Nbr of iterations for 
                                        % the optimization algorithm (opt)
fS.nWarmupIterations = 5; % Nbr of warmup iterations in the opt
fS.Ns = 50; fS.NsRBMC = 50; % Number of samples used in opt and 
                            % for the simple RBMC variance estimation
fS.solveMethod = 'PCG'; % 'Chol'; % % Method for solving matrix equations
fS.PCGtol = 1e-9; fS.PCGmaxiter = 500; fS.maxIcholEps = 100; % PCG settings
fS.updatePrecond = 1; % If 1, the pre-conditioner is frequently updated
fS.precondUpdateFreq = 10; % How often update preconditioner?
fS.gradGain = 1e-4; % Gradient gain in the SGD, only used during warmup
fS.doMargLikComp = 0; % Option 1 is not yet supported
fS.useHutchRandomPool = 0; % 1 means that a pool of random samples is created
                           % and samples will be reused for efficient PCG
fS.SHutchRandomPool = 100; % Number of samples in the pool
fS.useHess = 1; % If 1 use second order "Hessian" estimates in the opt
fS.hessAlphaStart = .9; % Starting value for opt parameter alpha^(j)
fS.maxDelta = 3; % Maximum change for any log-parameter during opt iteration
fS.gam1 = 0.2;  % Opt parameter gamma_1
fS.gam2Start = 0.9; % Opt param gamma_2 (Only used until startIterRobbinsMonro)
fS.momPrior = .5; % Opt parameter alpha_mom
fS.useRobbinsMonro = 1; % If 1 turns off the use of the Hessian in late iterations
fS.startIterRobbinsMonro = 100; % After which iteration to turn off the Hessian
fS.nPolyakVal = 4; % Nbr of the last iterations to base the final estimate on
fS.normalizeX = 1; % 0; % % If 1, the design matrix X is normalized
% Settings for non-spatial preoptimizing the noise parameters
fS.preOptNoiseParams = 1; fS.preOptNoiseMaxiter = 100;
fS.gradGainNoise = 1e-3; fS.useHessNoise = 0; fS.momNoise = 0.0;
% Whether to check AR process stability and tolerance for that
fS.checkARStability = 1; fS.ARCheckUnitCircleOffset = 1e-4;


%% Run
% SPMMat to be used as input
SPMMatPath = [fS.dataPath,fS.SPMResultsFolder,'/'];
load([SPMMatPath,'SPM.mat']);
  
% Create output folder
resultsPath = [fS.outputPath,fS.ResultsFolder];
mkdir(resultsPath);

% Slice independent parameters
X = SPM.xX.X; 

gSF = SPM.xGX.gSF;
fS.K = size(X,2); % number of regressors, code not safe for K = 1
fS.KSpat = fS.K - 7; % number of regressors with spatial priors
fS.T = size(X,1);
fS.sz = SPM.xY.VY(1).dim;
fS.P = SPM.PPM.AR_P; % 0; % 
fS.isGLMNoiseAR = (fS.P > 0); % 0; %
fS.sf = max(SPM.xBF.bf(:,1))/SPM.xBF.dt;

% 2d or 3d
if contains(fS.SPMResultsFolder,'2')
  fS.ndim = 2; fS.sz = fS.sz(1:2);
  sliceNbrs = unique(SPM.xVol.XYZ(3,:));
elseif contains(fS.SPMResultsFolder,'3')
  fS.ndim = 3; sliceNbrs = 0;
end

% Number of slices (or just one whole volume)
fS.nParcels = length(SPM.SVB);

% Ouputmat
Output.fS = fS; Output.X = X; Output.gSF = gSF; Output.sliceNbrs = sliceNbrs;
K = fS.K; T = fS.T; P = fS.P; sz = fS.sz; ndim = fS.ndim;

% Normalize X
if fS.normalizeX
  XOld = X;
  XNorm = sqrt(1/fS.T*dot(XOld,XOld))'; % sqrt(diag(Xold'*Xold))';
  meanHrfFirstRegNorm = mean(XNorm(1:2:fS.KSpat-1));
  XNorm([1:2:fS.KSpat-1,fS.K])=1;
  nuisRegInd = [2:2:fS.KSpat,fS.KSpat+1:fS.K];
  XNorm(nuisRegInd) = XNorm(nuisRegInd) ./ meanHrfFirstRegNorm;
  X = XOld ./ XNorm';
  Output.XNorm = XNorm;
  Output.XOld = XOld;
  Output.X = X;
end

% List of models
fMRIModels = cell(fS.nParcels,1);
for i = 1:fS.nParcels; fMRIModels{i} = FMRIModel(fS);end

% Start parpool
if fS.doParallel; pool = parpool; end

%% Start loop over slices/parcels
for iParcel = 1:fS.nParcels % end around line 202
% iParcel = 1; % 11; %
% iParcel = 37;
disp(['Running on slice/parcel nbr ',num2str(iParcel),...
                           ', Time: ',datestr(now,'dd-mmm-yyyy HH:MM:SS')]);

tic

% Get data and mask, check for voxels with zero signal
fM = fMRIModels{iParcel}; 
fM.X = X;
fM.Y = SPM.SVB(iParcel).Y;
if fS.ndim == 2
  mask = SPM.xVol.XYZ(1:2,SPM.xVol.XYZ(3,:) == sliceNbrs(iParcel));
  maskInd = sub2ind(fS.sz,mask(1,:)',mask(2,:)');
elseif fS.ndim == 3
  mask = SPM.xVol.XYZ(:,find(SPM.xVol.labels==iParcel));
  maskInd = sub2ind(fS.sz,mask(1,:)',mask(2,:)',mask(3,:)');
end
zeroSignalMask = (std(fM.Y) < 1e-14);
if any(zeroSignalMask)
  disp(['Warning: Voxels with zero signal inside mask. Removing ',...
    num2str(sum(zeroSignalMask)),' voxels.']);
  maskInd = maskInd(~zeroSignalMask);
  fM.Y = fM.Y(~zeroSignalMask);
end
mask_mat = false(fS.sz);
mask_mat(maskInd) = 1;
fM.maskInd = maskInd;
fM.N = size(fM.Y,2); N = fM.N;

% Parameter settings
fM.priorList = cell(fM.fS.K,1);

for k = 1:fS.KSpat
  fM.priorList{k} = feval(fS.hrfPrior{1},fM,k,fS.hrfPrior{2},fS.hrfPrior{3});
end
for k = fS.KSpat+1:fS.K-1 % head motion parameters
  fM.priorList{k} = feval(fS.hmpPrior{1},fM,k,fS.hmpPrior{2},fS.hmpPrior{3});
end
fM.priorList{fS.K} = feval(fS.interceptPrior{1},fM,fS.K,fS.interceptPrior{2},...
                           fS.interceptPrior{3});

% Initializations
fM.doPrecalculations;
fM.nM = NoiseModel(fM);
fM.initialize;

%% GD loop
while(fM.iter <= fS.maxiter)
    fM.gradientStep;
end

%% Compute final posterior
fM.doPolyakAveraging;
fM.computePostMeanAndStd;

%% Store results
if fS.normalizeX
  Output.b(iParcel).w = fM.M ./ XNorm;
  XNormLongInv = 1 ./ repmat(XNorm,N,1);
  XNormInvMat = spdiags(XNormLongInv,0,N*K,N*K);
  Output.b(iParcel).wCovMat = (XNormInvMat * fM.wCovMat) * XNormInvMat;
else
  Output.b(iParcel).w = fM.M;
  Output.b(iParcel).wCovMat = fM.wCovMat;
end

tau2 = fM.priorList{1}.tau2; kappa2 = fM.priorList{1}.kappa2;
tau2Vec = fM.priorList{1}.tau2Vec'; kappa2Vec = fM.priorList{1}.kappa2Vec';
for k = 2:fM.fS.K
    tau2 = [tau2;fM.priorList{k}.tau2];
    kappa2 = [kappa2;fM.priorList{k}.kappa2];
    tau2Vec = [tau2Vec;fM.priorList{k}.tau2Vec'];
    kappa2Vec = [kappa2Vec;fM.priorList{k}.kappa2Vec'];
end
Output.b(iParcel).fMThin = fM.saveModelThin;
Output.b(iParcel).tau2 = tau2; Output.b(iParcel).tau2Vec = tau2Vec;
Output.b(iParcel).kappa2 = kappa2; Output.b(iParcel).kappa2Vec = kappa2Vec;
Output.b(iParcel).lambda = fM.nM.lambda; Output.b(iParcel).lambdaVec = fM.nM.lambdaVec;
Output.b(iParcel).marglikVec = fM.marglikVec;
Output.b(iParcel).maskInd = fM.maskInd;

if strcmp(fM.priorList{1}.priorName,'Matern2Aniso')  
  hx = fM.priorList{1}.hx; hy = fM.priorList{1}.hy;
  hxVec = fM.priorList{1}.hxVec'; hyVec = fM.priorList{1}.hyVec';
  for k = 2:fM.fS.K
    if strcmp(fM.priorList{k}.priorName,'Matern2Aniso')
      hx = [hx;fM.priorList{k}.hx];
      hy = [hy;fM.priorList{k}.hy];
      hxVec = [hxVec;fM.priorList{k}.hxVec'];
      hyVec = [hyVec;fM.priorList{k}.hyVec'];
    end
  end
  Output.b(iParcel).hx = hx; Output.b(iParcel).hxVec = hxVec;
  Output.b(iParcel).hy = hy; Output.b(iParcel).hyVec = hyVec;
end

if fM.fS.isGLMNoiseAR
  Output.b(iParcel).aVec = fM.nM.aVec; Output.b(iParcel).a = fM.nM.a;
end

Output.b(iParcel).maxiter = fM.fS.maxiter;
Output.b(iParcel).timeTaken = toc;
Output.b(iParcel).timeStamp = datestr(now,'dd-mmm-yyyy HH:MM:SS');
save([fS.outputPath,fS.ResultsFolder,'/Output.mat'],'Output','-v7.3');
disp(['Results written to ',fS.outputPath,fS.ResultsFolder,'/Output.mat']);

end
if fS.doParallel; delete(pool); end
disp('done.');
