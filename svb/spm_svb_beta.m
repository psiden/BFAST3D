function [block] = spm_svb_beta(Y,block)
% Modified version for SVB, Per Siden 2016-06-09
%
% Variational Bayes for GLM-AR modelling in a block - Update beta 
% FORMAT [block] = spm_vb_beta (Y,block)
%
% Y             [T x N] time series 
% block         data structure 
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_beta.m 2451 2008-11-10 16:20:32Z lee $

if block.verbose
    disp('Updating beta');
end

N=block.N;
p=block.p;

if block.SVB.inWarmup
    Ns = block.SVB.NsWarmup;
else
    Ns = block.SVB.Ns;
end

for j = 1:p,
    subblock_p = j:p:N*p;

    H = mean(dot(block.SVB.aSampSim(subblock_p,1:Ns),(block.Da*block.SVB.aSampSim(subblock_p,1:Ns))));

    % Equation 16 in paper VB4
    block.b_beta(j)    = 1./(H./2 + 1./block.b_beta_prior(j));
    block.mean_beta(j) = block.c_beta(j)*block.b_beta(j);

    % Try to make beta converge faster by quadratic approximation 
    % y = ax^2 + bx + c
    if block.SVB.it > 3 && mod(block.SVB.it,2)==0
        a0 = block.SVB.betaSave(j,block.SVB.it-2);
        a1 = block.SVB.betaSave(j,block.SVB.it-1);
        a2 = block.mean_beta(j);
        if abs(a2 - a1) < abs(a1 - a0)
            a = .5*a2 - a1 + .5*a0;
            b = .5*a2 - .5*a0;
            c = a1;
            new_beta = -(b^2 - 4*a*c)/(4*a);
        else
            new_beta = a1 + block.SVB.betaNonConvexJumpSize*(a2-a1);
        end
        new_beta = min(new_beta,a2*block.SVB.betaMaxChangeFactor);
        new_beta = max(new_beta,a2/block.SVB.betaMaxChangeFactor);
        block.mean_beta(j) = new_beta;
    end
end