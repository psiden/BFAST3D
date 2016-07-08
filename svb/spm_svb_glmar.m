function [block] = spm_svb_glmar(Y,block)
% Modified version for SVB, Per Siden 2016-06-09
%
% Variational Bayes for GLM-AR modelling in a block of fMRI
% FORMAT [block] = spm_vb_glmar (Y,block)
%
% Y     -  [T x N] time series with T time points, N voxels
%
% block -  data structure containing the following fields:
%
%          .X              [T x k] the design matrix
%          .p              order of AR model
%          .D              [N x N] spatial precision matrix
%                          (see spm_vb_set_priors.m)
%
%          The above fields are mandatory. The fields below are
%          optional or are filled in by this function.
%
%          OPTIMISIATION PARAMETERS:
%
%          .tol            termination tolerance (default = 0.01% increase in F)
%          .maxits         maximum number of iterations (default=4)
%          .verbose        '1' for description of actions (default=1)
%          .update_???     set to 1 to update parameter ??? (set to 0 to fix)
%                          eg. update_alpha=1; % update prior precision on W
%
%          ESTIMATED REGRESSION COEFFICIENTS:
%
%          .wk_mean        [k x N] VB regression coefficients
%          .wk_ols         [k x N] OLS "  "
%          .w_cov          N-element cell array with entries [k x k]
%          .w_dev          [k x N] standard deviation of regression coeffs
%
%          ESTIMATED AR COEFFICIENTS:
%
%          .ap_mean        [p x N] VB AR coefficients
%          .ap_ols         [p x N] OLS AR coefficients
%          .a_cov          N-element cell array with entries [p x p]
%
%          ESTIMATED NOISE PRECISION:
%
%          .b_lambda       [N x 1] temporal noise precisions
%          .c_lambda
%          .mean_lambda
%
%          MODEL COMPARISON AND COEFFICIENT RESELS:
%
%          .gamma_tot      [k x 1] Coefficient RESELS
%          .F              Negative free energy (used for model selection)
%          .F_record       [its x 1] record of F at each iteration
%          .elapsed_seconds  estimation time
%          PRIORS:
%
%          .b_alpha        [k x 1] spatial prior precisions for W
%          .c_alpha
%          .mean_alpha
%
%          .b_beta         [p x 1] spatial prior precisions for AR
%          .c_beta
%          .mean_beta
%
%          .b              [k x N] prior precision matrix
%
%          HYPERPRIORS:
%
%          .b_alpha_prior   priors on alpha
%          .c_alpha_prior
%
%          .b_beta_prior    priors on beta
%          .c_beta_prior
%
%          .b_lambda_prior  priors on temporal noise precisions
%          .c_lambda_prior
%
%          There are other fields that are used internally
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_glmar.m 5219 2013-01-29 17:07:07Z spm $


t0 = clock;

%-Set defaults
%--------------------------------------------------------------------------
[T,N]   = size(Y);
block.T = T;
block.N = N;

try
    X = block.X;
catch
    error('Mandatory field X missing.');
end

[tmp,k] = size(X);
if tmp ~= T
    error('X is not of compatible size to Y.');
end

block.k = k;
p = block.p;

if ~isfield(block,'Dw')
    disp('Error in spm_vb_glmar: mandatory field Dw missing');
    return
end

%-Initialise block
%--------------------------------------------------------------------------
F      = 0;
last_F = 0;
block = spm_svb_init_block(Y,block);

%-Run Variational Bayes
%--------------------------------------------------------------------------
if block.verbose
    disp(' ');
    disp('Starting VB-GLM-AR-BLOCK');
end

block.SVB.timeTakenVec = zeros(block.maxits,1);
block = spm_svb_init(Y,block);
if block.SVB.ParallelGMRFSampling
    pool = parpool;
end

do_plot = 0;

block.SVB.alphaNonConvexJumpSize = 20;
block.SVB.alphaMaxChangeFactor = 5;
block.SVB.alphaSave = zeros(k,block.maxits);

if p > 0;
    block.SVB.betaNonConvexJumpSize = 20;
    block.SVB.betaMaxChangeFactor = 5;
    block.SVB.betaSave = zeros(p,block.maxits);
end;

for it = 1:block.maxits, % Loop over iterations
    
    block.SVB.it = it;
    if it == (block.SVB.warmupIterations + 1)
        block.SVB.inWarmup = 0;

        block.SVB.wSampSim(:,(block.SVB.NsWarmup+1):block.SVB.Ns) = ...
              repmat(block.w_mean,1,block.SVB.Ns-block.SVB.NsWarmup);
        block.SVB.wSampSimP = block.SVB.wSampSim(block.SVB.HwInd,:);
        block.SVB.wSampSimP = block.SVB.wSampSimP(block.SVB.BI,:);

        block.SVB.aSampSim(:,(block.SVB.NsWarmup+1):block.SVB.Ns) = ...
              repmat(block.a_mean,1,block.SVB.Ns-block.SVB.NsWarmup);
        block.SVB.aSampSimP = block.SVB.aSampSim(block.SVB.HaInd,:);
        block.SVB.aSampSimP = block.SVB.aSampSimP(block.SVB.JI,:);
    end

    if block.update_w
        block = spm_svb_w(Y,block);
    end
    if (block.p>0) && (block.update_a)
        block = spm_svb_a(Y,block);
    end
    if block.update_lambda
        block = spm_svb_lambda(Y,block);
    end
    if block.update_alpha
        block = spm_svb_alpha(Y,block);
    end
    if (block.p>0) && (block.update_beta)
        block = spm_svb_beta(Y,block);
    end
    if block.update_F
        [F, Lav, KL] = spm_vb_F (Y,block);
    end
    if block.verbose
        disp(sprintf('Iteration %d, F=%1.2f',it,F));
    end

    if block.update_F
        block.F_record(it)=F;
        delta_F=F-last_F;
        if it > 2
            if delta_F < 0
                disp(sprintf('********** Warning: decrease in F of %1.4f per cent *************',100*(delta_F/F)));
            end;     
        end
        last_F = F;
    end
    
    block.SVB.alphaSave(:,it) = block.mean_alpha;
    block.SVB.FSave(it) = F;
    if p > 0
        block.SVB.betaSave(:,it) = block.mean_beta;
    end
    if do_plot
        afig = figure(98);
        set(afig,'Position',[100 1000 1000 1000]);
        scK = ceil(sqrt(k));
        for kk = 1:k
            subplot(scK,scK,kk);
            plot(block.SVB.alphaSave(kk,1:it)); title(strcat('alpha',num2str(kk)));
        end
        drawnow
        if p > 0
            bfig = figure(99);
            set(bfig,'Position',[1200 1000 600 600]);
            scP = ceil(sqrt(p));
            for pp = 1:p
                subplot(scP,scP,pp);
                plot(block.SVB.betaSave(pp,1:it)); title(strcat('beta',num2str(pp)));
            end
            drawnow
        end
    end
    if block.verbose
        disp(datestr(now,'HH:MM:SS'));
    end
    block.SVB.timeTakenVec(it) = etime(clock,t0);
end

if block.update_F
    block.F   = F;
    block.Lav = Lav;
    block.KL  = KL;
end

for n = 1:N,
    block.w_cov{n} = cov(block.SVB.wSampSim(((n-1)*k+1):n*k,:)');
end
    
block = spm_svb_gamma(Y,block);

if block.SVB.ParallelGMRFSampling
    delete(pool);
end

block.elapsed_seconds = etime(clock,t0);
