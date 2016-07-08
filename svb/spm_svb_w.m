function [block] = spm_svb_w(Y,block)
% Modified version for SVB, Per Siden 2016-06-09
%
% Variational Bayes for GLM-AR modelling in a block - update w
% FORMAT [block] = spm_vb_w (Y,block)
%
% Y             [T x N] time series 
% block         data structure 
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_w.m 2451 2008-11-10 16:20:32Z lee $


if block.verbose
    disp('Updating w');
end

P = block.p;
N = block.N;
K = block.k;
if block.SVB.inWarmup
    Ns = block.SVB.NsWarmup;
else
    Ns = block.SVB.Ns;
end


if P <= 0 
    
    % w: compute precision
    block.SVB.BData = kron(spdiags(block.mean_lambda,0,N,N),block.SVB.XTX);
    for k = 1:K
        block.SVB.BAlphaList{k} = block.mean_alpha(k) * block.SVB.BList{k};
        block.SVB.GwAlphaList{k} = sqrt(block.mean_alpha(k)) * block.SVB.GwList{k};
    end
    block.SVB.BTilde = block.SVB.BData(block.SVB.HwInd,block.SVB.HwInd) + ...
                                         blkdiag(block.SVB.BAlphaList{:});
                
    % w: compute BTilde*mean
    block.SVB.BTildeMuMat = bsxfun(@times,block.mean_lambda,block.SVB.YTX);
    
else

    % MC estimations
    BDataVecTempEst = zeros(K,K,N,Ns);
    BTildeMuMatEst = zeros(N,K,Ns);
    for i = 1:Ns  

        aTemp = reshape(block.SVB.aSampSim(:,i)',P,N);

        % w: precision pre-steps
        for j = 1:K
            for k = 1:K
                block.SVB.ASA(j,k,:) = dot(aTemp,squeeze(block.SVB.S(:,j,k,:)) * ...
                                            aTemp,1)';
            end
            block.SVB.RA(:,j,:) = permute(squeeze(block.SVB.R(:,j,:)) * ...
                                        aTemp,[1,3,2]);
        end
        BDataVecTempEst(:,:,:,i) = repmat(block.SVB.XTX,[1,1,N]) - block.SVB.RA - ...
                                    permute(block.SVB.RA,[2,1,3]) + block.SVB.ASA;

        % w: compute BTilde*mean
        if P <= 1
            for j = 1:K
                block.SVB.DA(1,j,:) = squeeze(block.SVB.D(1,j,:,:)) .* aTemp';
            end
            for j = 1:K
                BTildeMuMatEst(:,j,i) = block.mean_lambda .* (block.SVB.YTX(:,j) - ...
                       aTemp' .* squeeze(block.SVB.B(:,j,:) - block.SVB.DA(:,j,:)));
            end
        else
            for p = 1:P
                for j = 1:K
                    block.SVB.DA(p,j,:) = dot(squeeze(block.SVB.D(p,j,:,:)),aTemp);
                end
            end
            for j = 1:K
                BTildeMuMatEst(:,j,i) = block.mean_lambda .* (block.SVB.YTX(:,j) - ...
                       dot(aTemp,squeeze(block.SVB.B(:,j,:) - block.SVB.DA(:,j,:)))');
            end
        end
    end
    block.SVB.BDataVecTemp = mean(BDataVecTempEst,4);
    block.SVB.BTildeMuMat = mean(BTildeMuMatEst,3);
    
    % w: compute precision
    block.SVB.BDataVec = zeros(K,K,N);
    for j = 1:K
        block.SVB.BDataVec(:,j,:) = bsxfun(@times,block.mean_lambda,...
                                        squeeze(block.SVB.BDataVecTemp(:,j,:))')';
    end
    block.SVB.BData = sparse(block.SVB.iBData,block.SVB.jBData,block.SVB.BDataVec(:));
    for k = 1:K
        block.SVB.BAlphaList{k} = block.mean_alpha(k) * block.SVB.BList{k};
        block.SVB.GwAlphaList{k} = sqrt(block.mean_alpha(k)) * block.SVB.GwList{k};
    end
    block.SVB.BTilde = block.SVB.BData(block.SVB.HwInd,block.SVB.HwInd) + ...
                                        blkdiag(block.SVB.BAlphaList{:});

end

% Reorder and update Preconditioner
block.SVB.BTildeP = block.SVB.BTilde(block.SVB.BI,block.SVB.BI);
try
    block.SVB.icholBTildeP = ichol(block.SVB.BTildeP);
catch
    if block.SVB.it == 1
        disp('Warning: ichol failed at first iteration');
        block.SVB.icholBTildeP = spdiags(sqrt(diag(block.SVB.BTildeP),0,N,N));
    else
        disp('Warning: ichol failed');
    end
end  
    
% w: compute mean  
b = block.SVB.BTildeMuMat(:);
bP = b(block.SVB.BI);
[block.SVB.w_meanP,flag,relres,iter] = pcg(block.SVB.BTildeP,bP,block.SVB.MeanPCGTol,...
                    500,block.SVB.icholBTildeP,block.SVB.icholBTildeP',block.SVB.w_meanP);
block.SVB.PCGwIterSave = [block.SVB.PCGwIterSave,iter];
                
block.w_mean(block.SVB.BI) = block.SVB.w_meanP;
block.w_mean(block.SVB.HwInd) = block.w_mean;
block.SVB.w = reshape(block.w_mean,K,N);

% samples for expectation approximation
block.SVB.cholBData = chol(block.SVB.BData);
block.SVB.bSamp = blkdiag(block.SVB.GwAlphaList{:})' * block.SVB.wSampPriorRand + ...
      block.SVB.cholBData(block.SVB.HwInd,block.SVB.HwInd)' * block.SVB.wSampDataRand + ...
      repmat(block.SVB.BTildeMuMat(:),1,block.SVB.Ns);
block.SVB.bSampP = block.SVB.bSamp(block.SVB.BI,:);
if block.SVB.ParallelGMRFSampling
    [samples,iterSave] = spm_svb_parallel_GMRF_sampling(block.SVB.BTildeP,block.SVB.bSampP,...
                                block.SVB.SimPCGTol,block.SVB.icholBTildeP,block.SVB.wSampSimP(:,1:Ns));
    block.SVB.wSampSimP(:,1:Ns) = samples;
    block.SVB.PCGwIterSave = [block.SVB.PCGwIterSave,iterSave];
else
    for i = 1:Ns
        [wP,flag,relres,iter] = pcg(block.SVB.BTildeP,block.SVB.bSampP(:,i),...
                block.SVB.SimPCGTol,500,block.SVB.icholBTildeP,block.SVB.icholBTildeP',block.SVB.wSampSimP(:,i));
        block.SVB.PCGwIterSave = [block.SVB.PCGwIterSave,iter];
        block.SVB.wSampSimP(:,i) = wP;
    end
end
block.SVB.wSampSim(block.SVB.BI,:) = block.SVB.wSampSimP;      
block.SVB.wSampSim(block.SVB.HwInd,:) = block.SVB.wSampSim;
