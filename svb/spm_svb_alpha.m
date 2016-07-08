function [block] = spm_svb_alpha(Y,block)
% Modified version for SVB, Per Siden 2016-06-09
%
% Variational Bayes for GLM-AR models - Update alpha 
% FORMAT [block] = spm_vb_alpha (Y,block)
%
% Y             [T x N] time series 
% block         data structure 
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_alpha.m 2451 2008-11-10 16:20:32Z lee $

if block.verbose
    disp('Updating alpha');
end

N=block.N;
k=block.k;
  
if block.SVB.inWarmup
    Ns = block.SVB.NsWarmup;
else
    Ns = block.SVB.Ns;
end

for j = 1:k,
    subblock_k = j:k:N*k;

    H = mean(dot(block.SVB.wSampSim(subblock_k,1:Ns),(block.Dw*block.SVB.wSampSim(subblock_k,1:Ns))));
    
    % Equation 15 in paper VB4
    block.b_alpha(j)    = 1./(H./2 + 1./block.b_alpha_prior(j));
    block.mean_alpha(j) = block.c_alpha(j)*block.b_alpha(j);
    
    
    % Try to make alpha converge faster by quadratic approximation 
    % y = ax^2 + bx + c
    if block.SVB.it > 5 && mod(block.SVB.it,2)==0

        % regular quadratic approximation
        a0 = block.SVB.alphaSave(j,block.SVB.it-2);
        a1 = block.SVB.alphaSave(j,block.SVB.it-1);
        a2 = block.mean_alpha(j);
        if abs(a2 - a1) < abs(a1 - a0)
            a = .5*a2 - a1 + .5*a0;
            b = .5*a2 - .5*a0;
            c = a1;
            new_alpha = -(b^2 - 4*a*c)/(4*a);
        else
            new_alpha = a1 + block.SVB.alphaNonConvexJumpSize*(a2-a1);
        end
        new_alpha = min(new_alpha,a2*block.SVB.alphaMaxChangeFactor);
        new_alpha = max(new_alpha,a2/block.SVB.alphaMaxChangeFactor);
        block.mean_alpha(j) = new_alpha;
    end
end
    
