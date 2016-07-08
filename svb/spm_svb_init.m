function [block] = spm_svb_init(Y,block)
% Initialize Spatial Variational Bayes for GLM-AR model
%
% Author: Per Siden
% Created: 2016-06-09

K=block.k;
P=block.p;
N=block.N;
T=block.T;
X=block.X;

block.SVB.Ns = 100;
block.SVB.NsWarmup = 5;
block.SVB.MeanPCGTol = 1e-8;
block.SVB.SimPCGTol = 1e-8;
block.SVB.PCGwIterSave = [];
block.SVB.PCGaIterSave = [];
block.SVB.warmupIterations = 10;
block.SVB.inWarmup = 1;

rngseed = 102;
if rngseed > 0
    rng('default');
    rng(rngseed);
end

if P <= 0

    block.SVB.XTX = block.XTX;
    block.SVB.YTY = dot(Y,Y);
    block.SVB.YTX = block.XTY';
    
    % Find reordering
    block.SVB.HwInd = reshape(1:(K*N),K,N)';
    block.SVB.HwInd = block.SVB.HwInd(:);
    block.SVB.iHwInd = zeros(N*K,1);
    block.SVB.iHwInd(block.SVB.HwInd) = 1:(N*K);
    block.SVB.BAlphaList = cell(K,1);
    block.SVB.GwAlphaList = cell(K,1);
    block.SVB.BData = kron(spdiags(block.mean_lambda,0,N,N),block.SVB.XTX);
    BTypes = {''};
    for k = 1:K
        if strcmp(block.priors.W,'Spatial - UGL');
            BTypes{k} = 'LI';
        end
    end
    [block.SVB.BList,block.SVB.GwList] = spm_svb_setupPrecMats(BTypes,...
                         N,block.SVB.sz,block.SVB.bmask,block.SVB.ndim);
    for k = 1:K
        block.SVB.BAlphaList{k} = block.mean_alpha(k) * block.SVB.BList{k};
        block.SVB.GwAlphaList{k} = sqrt(block.mean_alpha(k)) * block.SVB.GwList{k};
    end
    block.SVB.BTilde = block.SVB.BData(block.SVB.HwInd,block.SVB.HwInd) + ...
                                         blkdiag(block.SVB.BAlphaList{:});
    block.SVB.BI = amd(block.SVB.BTilde);
    block.SVB.iBI = zeros(N*K,1);
    block.SVB.iBI(block.SVB.BI) = 1:(N*K);
    block.SVB.BTildeP = block.SVB.BTilde(block.SVB.BI,block.SVB.BI);
    block.SVB.icholBTildeP = ichol(block.SVB.BTildeP);
    block.SVB.NGw = size(blkdiag(block.SVB.GwAlphaList{:}),1);
    block.SVB.BTildeMuMat = zeros(N,K); % BTilde * mu on (N x K)-form
    block.SVB.w_meanP = block.w_mean(block.SVB.HwInd);
    block.SVB.w_meanP = block.SVB.w_meanP(block.SVB.BI);
    block.SVB.wSampPriorRand = randn(block.SVB.NGw,block.SVB.Ns);
    block.SVB.wSampDataRand = randn(K*N,block.SVB.Ns);
    block.SVB.wSampSim = repmat(block.w_mean,1,block.SVB.Ns);
    block.SVB.wSampSimP = block.SVB.wSampSim(block.SVB.HwInd,:);
    block.SVB.wSampSimP = block.SVB.wSampSimP(block.SVB.BI,:);
    
else
    
    block.SVB.w = reshape(block.w_mean,K,N);
    block.SVB.a = reshape(block.a_mean,P,N);
    
    % lagged X
    block.SVB.XTilde = zeros(P,T-P,K);
    for t = 1:(T-P)
        block.SVB.XTilde(:,t,:) = permute(block.X((t+P-1):-1:t,:),[1,3,2]);
    end
    block.SVB.X = block.X((P+1):T,:);
    
    % lagged Y
    block.SVB.d = zeros(P,T-P,N);
    for t = 1:(T-P)
        block.SVB.d(:,t,:) = permute(Y((t+P-1):-1:t,:),[1,3,2]); 
    end
    block.SVB.Y = Y((P+1):T,:);
    
    % voxel independent
    block.SVB.XTX = block.SVB.X'*block.SVB.X;
    block.SVB.R = zeros(K,K,P);
    block.SVB.S = zeros(P,K,K,P);
    for j = 1:K
        block.SVB.R(:,j,:) = permute(block.SVB.X' * ...
                                squeeze(block.SVB.XTilde(:,:,j))',[1,3,2]);
        for k = 1:K
            block.SVB.S(:,j,k,:) = permute(squeeze(block.SVB.XTilde(:,:,j)) * ...
                                        squeeze(block.SVB.XTilde(:,:,k))',[3,1,2,4]);
        end
    end
    
    % voxel dependent
    block.SVB.YTX = block.SVB.Y'*block.SVB.X;
    block.SVB.YTY = dot(block.SVB.Y,block.SVB.Y);
    block.SVB.B = zeros(P,K,N); % B = Gxy + Rxy in Penny (2005) (vb3)
    block.SVB.D = zeros(P,K,P,N);
    block.SVB.ddT = zeros(P,P,N);
    block.SVB.YTdT = zeros(N,P);
    
    for n = 1:N
        block.SVB.YTXTilde = zeros(P,K);
        block.SVB.dX = zeros(P,K);
        block.SVB.dXTilde = zeros(P,K,P);
        for p = 1:P
            block.SVB.YTXTilde(p,:) = block.SVB.Y(:,n)' * ...
                                        squeeze(block.SVB.XTilde(p,:,:));
            block.SVB.dX(p,:) = squeeze(block.SVB.d(p,:,n)) * block.SVB.X;
            block.SVB.dXTilde(:,:,p) = squeeze(block.SVB.d(:,:,n)) * ...
                                        squeeze(block.SVB.XTilde(p,:,:));
        end
        block.SVB.B(:,:,n) = block.SVB.YTXTilde + block.SVB.dX;
        block.SVB.D(:,:,:,n) = block.SVB.dXTilde;
        block.SVB.ddT(:,:,n) = squeeze(block.SVB.d(:,:,n)) * ...
                                    squeeze(block.SVB.d(:,:,n))';
        block.SVB.YTdT(n,:) = block.SVB.Y(:,n)' * ...
                                squeeze(block.SVB.d(:,:,n))';
    end

    % Find reorderings
    
    % w
    block.SVB.HwInd = reshape(1:(K*N),K,N)';
    block.SVB.HwInd = block.SVB.HwInd(:);
    block.SVB.iHwInd = zeros(N*K,1);
    block.SVB.iHwInd(block.SVB.HwInd) = 1:(N*K);
    block.SVB.BAlphaList = cell(K,1);
    block.SVB.GwAlphaList = cell(K,1);
    
    block.SVB.ASA = zeros(K,K,N);
    block.SVB.RA = zeros(K,K,N);
    for j = 1:K
        for k = 1:K
            block.SVB.ASA(j,k,:) = dot(block.SVB.a,squeeze(block.SVB.S(:,j,k,:)) * ...
                                        block.SVB.a,1)';
        end
        block.SVB.RA(:,j,:) = permute(squeeze(block.SVB.R(:,j,:)) * ...
                                    block.SVB.a,[1,3,2]);
    end
    block.SVB.BDataVecTemp = repmat(block.SVB.XTX,[1,1,N]) - block.SVB.RA - ...
                                permute(block.SVB.RA,[2,1,3]) + block.SVB.ASA;
    block.SVB.BDataVec = zeros(K,K,N);
    for j = 1:K
        block.SVB.BDataVec(:,j,:) = bsxfun(@times,block.mean_lambda,...
                                        squeeze(block.SVB.BDataVecTemp(:,j,:))')';
    end
    block.SVB.iBData = reshape(repmat(1:K,N*K,1)',N*K*K,1) + ...
                            reshape(repmat(K*(0:N-1),K*K,1),N*K*K,1);
    block.SVB.jBData = reshape(repmat(1:N*K,K,1),N*K*K,1);
    block.SVB.BData = sparse(block.SVB.iBData,block.SVB.jBData,block.SVB.BDataVec(:));
    BTypes = {''};
    for k = 1:K
        if strcmp(block.priors.W,'Spatial - UGL');
            BTypes{k} = 'LI';
        end
    end
    [block.SVB.BList,block.SVB.GwList] = spm_svb_setupPrecMats(BTypes,...
                         N,block.SVB.sz,block.SVB.bmask,block.SVB.ndim);
    for k = 1:K
        block.SVB.BAlphaList{k} = block.mean_alpha(k) * block.SVB.BList{k};
        block.SVB.GwAlphaList{k} = sqrt(block.mean_alpha(k)) * block.SVB.GwList{k};
    end
    block.SVB.BTilde = block.SVB.BData(block.SVB.HwInd,block.SVB.HwInd) + ...
                                        blkdiag(block.SVB.BAlphaList{:});
        
    block.SVB.BI = amd(block.SVB.BTilde);
    block.SVB.iBI = zeros(N*K,1);
    block.SVB.iBI(block.SVB.BI) = 1:(N*K);
    block.SVB.BTildeP = block.SVB.BTilde(block.SVB.BI,block.SVB.BI);
    block.SVB.icholBTildeP = ichol(block.SVB.BTildeP);
    block.SVB.NGw = size(blkdiag(block.SVB.GwAlphaList{:}),1);
    block.SVB.DA = zeros(P,K,N);
    block.SVB.BTildeMuMat = zeros(N,K); % BTilde * mu on (N x K)-form
    block.SVB.w_meanP = block.w_mean(block.SVB.HwInd);
    block.SVB.w_meanP = block.SVB.w_meanP(block.SVB.BI);
    block.SVB.wSampPriorRand = randn(block.SVB.NGw,block.SVB.Ns);
    block.SVB.wSampDataRand = randn(K*N,block.SVB.Ns);
    block.SVB.wSampSim = repmat(block.w_mean,1,block.SVB.Ns);
    block.SVB.wSampSimP = block.SVB.wSampSim(block.SVB.HwInd,:);
    block.SVB.wSampSimP = block.SVB.wSampSimP(block.SVB.BI,:);
    
    % a
    block.SVB.HaInd = reshape(1:(P*N),P,N)';
    block.SVB.HaInd = block.SVB.HaInd(:);
    block.SVB.iHaInd = zeros(N*P,1);
    block.SVB.iHaInd(block.SVB.HaInd) = 1:(N*P);
    block.SVB.JBetaList = cell(P,1);
    block.SVB.GaBetaList = cell(P,1);
    
    block.SVB.WSW = zeros(P,P,N);
    block.SVB.DW = zeros(P,P,N);
    for p = 1:P
        for q = 1:P
            block.SVB.WSW(p,q,:) = dot(block.SVB.w,squeeze(block.SVB.S(p,:,:,q)) * block.SVB.w)';
            block.SVB.DW(p,q,:) = dot(squeeze(block.SVB.D(p,:,q,:)),block.SVB.w);
        end
    end
    block.SVB.JDataVecTemp = block.SVB.ddT - block.SVB.DW - ...
                                permute(block.SVB.DW,[2,1,3]) + block.SVB.WSW;
    block.SVB.JDataVec = zeros(P,P,N);
    if P <= 1
        block.SVB.JDataVec(:,1,:) = block.mean_lambda .* squeeze(block.SVB.JDataVecTemp(:,p,:));
    else
        for p = 1:P
            block.SVB.JDataVec(:,p,:) = bsxfun(@times,block.mean_lambda,...
                                    squeeze(block.SVB.JDataVecTemp(:,p,:))')';
        end
    end
    block.SVB.iJData = reshape(repmat(1:P,N*P,1)',N*P*P,1) + ...
                            reshape(repmat(P*(0:N-1),P*P,1),N*P*P,1);
    block.SVB.jJData = reshape(repmat(1:N*P,P,1),N*P*P,1);
    block.SVB.JData = sparse(block.SVB.iJData,block.SVB.jJData,block.SVB.JDataVec(:));
    JTypes = {''};
    for p = 1:P
        if strcmp(block.priors.A,'Spatial - UGL');
            JTypes{p} = 'LI';
        end
    end
    [block.SVB.JList,block.SVB.GaList] = spm_svb_setupPrecMats(JTypes,...
                         N,block.SVB.sz,block.SVB.bmask,block.SVB.ndim);
    for p = 1:P
        block.SVB.JBetaList{p} = block.mean_beta(p) * block.SVB.JList{p};
        block.SVB.GaBetaList{p} = sqrt(block.mean_beta(p)) * block.SVB.GaList{p};
    end
    block.SVB.JTilde = block.SVB.JData(block.SVB.HaInd,block.SVB.HaInd) + ...
                            blkdiag(block.SVB.JBetaList{:});
        
    block.SVB.JI = amd(block.SVB.JTilde);
    block.SVB.iJI = zeros(N*P,1);
    block.SVB.iJI(block.SVB.JI) = 1:(N*P);
    block.SVB.JTildeP = block.SVB.JTilde(block.SVB.JI,block.SVB.JI);
    block.SVB.icholJTildeP = ichol(block.SVB.JTildeP);
    block.SVB.NGa = size(blkdiag(block.SVB.GaBetaList{:}),1);
    block.SVB.RW = zeros(K,P,N);
    block.SVB.JTildeNuMat = zeros(N,P); % JTilde * nu on (N x P)-form
    block.SVB.a_meanP = block.a_mean(block.SVB.HaInd);
    block.SVB.a_meanP = block.SVB.a_meanP(block.SVB.JI);
    block.SVB.aSampPriorRand = randn(block.SVB.NGa,block.SVB.Ns);
    block.SVB.aSampDataRand = randn(P*N,block.SVB.Ns);
    block.SVB.aSampSim = repmat(block.a_mean,1,block.SVB.Ns);
    block.SVB.aSampSimP = block.SVB.aSampSim(block.SVB.HaInd,:);
    block.SVB.aSampSimP = block.SVB.aSampSimP(block.SVB.JI,:); 
    
    % lambda
    block.SVB.BW  = zeros(p,N);
    block.SVB.WRW = zeros(p,N);
    block.SVB.ddTA = zeros(p,N);
    block.SVB.DWA = zeros(p,N);
    block.SVB.WSWA = zeros(p,N);

end
