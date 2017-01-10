function [block] = spm_svb_a(Y,block)
% Modified version for SVB, Per Siden 2016-06-09, 2017-01-10
%
% Update AR coefficients in VB GLM-AR model 
% FORMAT [block] = spm_vb_a (Y,block)
%
% Y             [T x N] time series 
% block         data structure 
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_a.m 2451 2008-11-10 16:20:32Z lee $

if block.verbose
    disp('Updating a');
end

P=block.p;
K=block.k;
N=block.N;

if block.SVB.inWarmup
    Ns = block.SVB.NsWarmup;
else
    Ns = block.SVB.Ns;
end

% MC estimations
JDataVecTempEst = zeros(P,P,N,Ns);
JTildeNuMatEst = zeros(N,P,Ns);
for i = 1:Ns  

    wTemp = reshape(block.SVB.wSampSim(:,i)',K,N);

    % for precision
    for p = 1:P
        for q = 1:P
            block.SVB.WSW(p,q,:) = dot(wTemp,squeeze(block.SVB.S(p,:,:,q)) * ...
                                        wTemp)';
            block.SVB.DW(p,q,:) = dot(squeeze(block.SVB.D(p,:,q,:)),wTemp);
        end
    end
    JDataVecTempEst(:,:,:,i) = block.SVB.ddT - block.SVB.DW - ...
                                permute(block.SVB.DW,[2,1,3]) + block.SVB.WSW;
    % for mean
    for p = 1:P
        block.SVB.RW(:,p,:) = permute(squeeze(block.SVB.R(:,:,p)) * wTemp,[1,3,2]);
    end
    for p = 1:P
        JTildeNuMatEst(:,p,i) = block.mean_lambda .* (block.SVB.YTdT(:,p) - ...
             dot(wTemp,squeeze(block.SVB.B(p,:,:)) - squeeze(block.SVB.RW(:,p,:)))');
    end
end
block.SVB.JDataVecTemp = mean(JDataVecTempEst,4);
block.SVB.JTildeNuMat = mean(JTildeNuMatEst,3);

% a: compute precision
if P <= 1
    block.SVB.JDataVec(:,1,:) = block.mean_lambda .* squeeze(block.SVB.JDataVecTemp(:,p,:));
else
    for p = 1:P
        block.SVB.JDataVec(:,p,:) = bsxfun(@times,block.mean_lambda,...
                                squeeze(block.SVB.JDataVecTemp(:,p,:))')';
    end
end
block.SVB.JData = sparse(block.SVB.iJData,block.SVB.jJData,block.SVB.JDataVec(:));
for p = 1:P
    block.SVB.JBetaList{p} = block.mean_beta(p) * block.SVB.JList{p};
    block.SVB.GaBetaList{p} = sqrt(block.mean_beta(p)) * block.SVB.GaList{p};
end
block.SVB.JTilde = block.SVB.JData(block.SVB.HaInd,block.SVB.HaInd) + ...
                        blkdiag(block.SVB.JBetaList{:});

% Reorder and update Preconditioner
block.SVB.JTildeP = block.SVB.JTilde(block.SVB.JI,block.SVB.JI);
try
    block.SVB.icholJTildeP = ichol(block.SVB.JTildeP);
catch
    if block.SVB.it == 1
        disp('Warning: ichol failed at first iteration');
        block.SVB.icholJTildeP = spdiags(sqrt(diag(block.SVB.JTildeP),0,N,N));
    else
        disp('Warning: ichol failed');
    end
end  

% a: compute mean
ba = block.SVB.JTildeNuMat(:);
baP = ba(block.SVB.JI);     
[block.SVB.a_meanP,flag,relres,iter] = pcg(block.SVB.JTildeP,baP,block.SVB.MeanPCGTol,...
                    500,block.SVB.icholJTildeP,block.SVB.icholJTildeP',block.SVB.a_meanP);
block.SVB.PCGaIterSave = [block.SVB.PCGaIterSave,iter];
block.a_mean(block.SVB.JI) = block.SVB.a_meanP;
block.a_mean(block.SVB.HaInd) = block.a_mean;

block.SVB.PCGaIterSave = [block.SVB.PCGaIterSave,iter];

% samples for expectation approximation
block.SVB.cholJData = spm_svb_cholSafe(block.SVB.JData,1e-8);
block.SVB.baSamp = blkdiag(block.SVB.GaBetaList{:})' * block.SVB.aSampPriorRand + ...
      block.SVB.cholJData(block.SVB.HaInd,block.SVB.HaInd)' * block.SVB.aSampDataRand + ...
      repmat(block.SVB.JTildeNuMat(:),1,block.SVB.Ns);
block.SVB.baSampP = block.SVB.baSamp(block.SVB.JI,:);
if block.SVB.ParallelGMRFSampling
    [samples,iterSave] = spm_svb_parallel_GMRF_sampling(block.SVB.JTildeP,block.SVB.baSampP,...
                                block.SVB.SimPCGTol,block.SVB.icholJTildeP,block.SVB.aSampSimP(:,1:Ns));
    block.SVB.aSampSimP(:,1:Ns) = samples;
    block.SVB.PCGaIterSave = [block.SVB.PCGaIterSave,iterSave];
else
    for i = 1:Ns
        [aP,flag,relres,iter] = pcg(block.SVB.JTildeP,block.SVB.baSampP(:,i),...
                block.SVB.SimPCGTol,500,block.SVB.icholJTildeP,block.SVB.icholJTildeP',block.SVB.aSampSimP(:,i));
        block.SVB.PCGaIterSave = [block.SVB.PCGaIterSave,iter];
        block.SVB.aSampSimP(:,i) = aP;
    end
end
block.SVB.aSampSim(block.SVB.JI,:) = block.SVB.aSampSimP;      
block.SVB.aSampSim(block.SVB.HaInd,:) = block.SVB.aSampSim;

block.SVB.a = reshape(block.a_mean,P,N);
