%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run MCMC spatial model estimation
%
% INPUT:        dS - data settings struct
%               estimationMethod - string containing estimation method,
%                                  possible options: MCMC2D, MCMC3D.
%               VBMethod - string containing VB estimation method, possible 
%                          options: SVB2D, SVB3D.
%               samplingMethod - string containing method of sampling from
%               the GMRF, possible options: Cholesky, PCG.
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
function runMCMC(dS,estimationMethod,VBMethod,samplingMethod)
%% Run

% Read settings
outputPath = dS.outputPath; subjStr = dS.subjStr;

do_plot = 0;

% SPMMat to be used as input
SPMMatPath = strcat(outputPath,'sub-',subjStr,'/',VBMethod,'/');
load(strcat(SPMMatPath,'SPM.mat'));

% Create output folder
resultsPath = strcat(outputPath,'sub-',subjStr,'/',estimationMethod);
mkdir(resultsPath);

% Slice independent parameters
X = SPM.xX.X;
sz = SPM.xY.VY(1).dim;
gSF = SPM.xGX.gSF;
P = SPM.PPM.AR_P;
isGLMNoiseAR = (P > 0);
K = size(X,2); % number of regressors, code not safe for K = 1
T = size(X,1);

% 2D or 3D
if strcmp(estimationMethod,'MCMC2D')
    ndim = 2;
    sz = sz(1:2);
    sliceNbrs = unique(SPM.xVol.XYZ(3,:));
elseif strcmp(estimationMethod,'MCMC3D')
    ndim = 3;
    sliceNbrs = 0;
end

% Ouputmat
MCMC.a.estimationMethod = estimationMethod;
MCMC.a.samplingMethod = samplingMethod;
MCMC.a.X = X;
MCMC.a.sz = sz;
MCMC.a.gSF = gSF;
MCMC.a.P = P;
MCMC.a.isGLMNoiseAR = isGLMNoiseAR;
MCMC.a.K = K;
MCMC.a.T = T;
MCMC.a.ndim = ndim;
MCMC.a.sliceNbrs = sliceNbrs;

% Number of slices (or just one whole volume)
nLb = length(SPM.SVB);

for iLb = 1:nLb
    
    tic
    % Parameter settings
    Y = SPM.SVB(iLb).Y;
    if ndim == 2
        mask2d = SPM.xVol.XYZ(1:2,SPM.xVol.XYZ(3,:) == sliceNbrs(iLb));
        bmask = sub2ind(sz,mask2d(1,:)',mask2d(2,:)');
    elseif ndim == 3 
        mask3d = SPM.xVol.XYZ;
        bmask = sub2ind(sz,mask3d(1,:)',mask3d(2,:)',mask3d(3,:)');
    end
    bmask_mat = false(sz);
    bmask_mat(bmask) = 1;
    N = size(Y,2);
    
    BTypes = {''};
    for k = 1:K
        BTypes{k} = 'LI';
    end
    [BList,GwList] = setupPrecMats(BTypes,N,sz,bmask,ndim);
    if isGLMNoiseAR
        JTypes = {''};
        for p = 1:P
            JTypes{p} = 'LI';
        end
        [JList,GaList] = setupPrecMats(JTypes,N,sz,bmask,ndim);
    end

    % Initialize Gibbs
    niter = 10000; % nbr of Gibbs iterations
    warmup = 1000;
    thinningFactor = 5;
    
    if strcmp(samplingMethod,'PCG')
        updatePrecond = 1;
        PCGTol = 1e-8;
    end
    
    rngseed = 102;
    if rngseed > 0
        rng(rngseed);
    end

    % Set prior parameters
    % lambda
    u1 = 10;
    u2 = 0.1;
    u1Tilde = zeros(N,1);
    u2Tilde = ((T-P)/2 + u2)*ones(N,1);

    % alpha
    q1 = 10*ones(K,1);
    q2 = 0.1*ones(K,1);
    q1Tilde = zeros(K,1);
    q2Tilde = N/2 + q2;

    % beta
    if isGLMNoiseAR
        r1 = 10000*ones(P,1);
        r2 = 0.1*ones(P,1);
        r1Tilde = zeros(P,1);
        r2Tilde = N/2 + r2;
    end

    % Set start values
    lambda = u1*u2*ones(N,1);
    w = zeros(K,N);
    alpha = q1.*q2;
    if isGLMNoiseAR
        a = zeros(P,N);
        beta = r1.*r2;
    end

    % Results storage matrices
    thinningInd = zeros(warmup+niter,1);
    thinningInd([1:warmup,warmup+(1+mod((niter)-1,thinningFactor):thinningFactor:niter)]) = ...
                1:warmup+ceil((niter)/thinningFactor);
    nSavedDraws = nnz(thinningInd);

    wVec = zeros(K,N,nSavedDraws);
    wVec(:,:,1) = w;
    lambdaVec = zeros(N,nSavedDraws);
    lambdaVec(:,1) = lambda;
    alphaVec = zeros(K,nSavedDraws);
    alphaVec(:,1) = alpha;
    if isGLMNoiseAR
        aVec = zeros(P,N,nSavedDraws);
        aVec(:,:,1) = a;
        betaVec = zeros(P,nSavedDraws);
        betaVec(:,1) = beta;
    end

    % Precalculations
    if ~isGLMNoiseAR

        XTX = sparse(X'*X);
        YTY = dot(Y,Y);
        YTX = sparse(Y'*X);

        % Find reordering
        HwInd = reshape(1:(K*N),K,N)';
        HwInd = HwInd(:);
        BAlphaList = cell(K,1);
        GwAlphaList = cell(K,1);
        BData = kron(spdiags(lambda,0,N,N),XTX);
        for k = 1:K
            BAlphaList{k} = alpha(k) * BList{k};
            GwAlphaList{k} = sqrt(alpha(k)) * GwList{k};
        end
        BTilde = BData(HwInd,HwInd) + blkdiag(BAlphaList{:});
        BI = amd(BTilde);

        BTildeP = BTilde(BI,BI);
        M1 = ichol(BTildeP);
        NGw = size(blkdiag(GwAlphaList{:}),1);
        BTildeMuMat = zeros(N,K); % BTilde * mu on (N x K)-form
        wr = reshape(w',N*K,1); wrP = wr;

    else
        % lagged X
        XTilde = zeros(P,T-P,K);
        for t = 1:(T-P)
            XTilde(:,t,:) = permute(X((t+P-1):-1:t,:),[1,3,2]);
        end
        Xcut = X((P+1):T,:);

        % lagged Y
        d = zeros(P,T-P,N);
        for t = 1:(T-P)
            d(:,t,:) = permute(Y((t+P-1):-1:t,:),[1,3,2]); 
        end
        Y = Y((P+1):T,:);

        % voxel independent
        XTX = Xcut'*Xcut;
        R = zeros(K,K,P);
        S = zeros(P,K,K,P);
        for j = 1:K
            R(:,j,:) = permute(Xcut' * squeeze(XTilde(:,:,j))',[1,3,2]);
            for k = 1:K
                S(:,j,k,:) = permute(squeeze(XTilde(:,:,j)) * squeeze(XTilde(:,:,k))',[3,1,2,4]);
            end
        end

        % voxel dependent
        YTX = Y'*Xcut;
        YTY = dot(Y,Y);
        B = zeros(P,K,N); % B = Gxy + Rxy in Penny (2005) (vb3)
        D = zeros(P,K,P,N);
        ddT = zeros(P,P,N);
        YTdT = zeros(N,P);

        for n = 1:N
            YTXTilde = zeros(P,K);
            dX = zeros(P,K);
            dXTilde = zeros(P,K,P);
            for p = 1:P
                YTXTilde(p,:) = Y(:,n)' * squeeze(XTilde(p,:,:));
                dX(p,:) = squeeze(d(p,:,n)) * Xcut;
                dXTilde(:,:,p) = squeeze(d(:,:,n)) * squeeze(XTilde(p,:,:));
            end
            B(:,:,n) = YTXTilde + dX;
            D(:,:,:,n) = dXTilde;
            ddT(:,:,n) = squeeze(d(:,:,n)) * squeeze(d(:,:,n))';
            YTdT(n,:) = Y(:,n)' * squeeze(d(:,:,n))';
        end
        
        % Find reorderings
        % w
        HwInd = reshape(1:(K*N),K,N)';
        HwInd = HwInd(:);
        BAlphaList = cell(K,1);
        GwAlphaList = cell(K,1);

        ASA = zeros(K,K,N);
        RA = zeros(K,K,N);
        for j = 1:K
            for k = 1:K
                ASA(j,k,:) = dot(a,squeeze(S(:,j,k,:)) * a,1)';
            end
            RA(:,j,:) = permute(squeeze(R(:,j,:)) * a,[1,3,2]);
        end
        BDataVecTemp = repmat(XTX,[1,1,N]) -RA - permute(RA,[2,1,3]) + ASA;
        BDataVec = zeros(K,K,N);
        for j = 1:K
            BDataVec(:,j,:) = bsxfun(@times,lambda,squeeze(BDataVecTemp(:,j,:))')';
        end
        iBData = reshape(repmat(1:K,N*K,1)',N*K*K,1) + ...
                                reshape(repmat(K*(0:N-1),K*K,1),N*K*K,1);
        jBData = reshape(repmat(1:N*K,K,1),N*K*K,1);
        BData = sparse(iBData,jBData,BDataVec(:));
        for k = 1:K
            BAlphaList{k} = alpha(k) * BList{k};
            GwAlphaList{k} = sqrt(alpha(k)) * GwList{k};
        end
        BTilde = BData(HwInd,HwInd) + blkdiag(BAlphaList{:});

        BI = amd(BTilde);
        BTildeP = BTilde(BI,BI);
        M1 = ichol(BTildeP);
        NGw = size(blkdiag(GwAlphaList{:}),1);
        DA = zeros(P,K,N);
        BTildeMuMat = zeros(N,K); % BTilde * mu on (N x K)-form
        wr = reshape(w',N*K,1); wrP = wr;

        % a
        HaInd = reshape(1:(P*N),P,N)';
        HaInd = HaInd(:);
        JBetaList = cell(P,1);
        GaBetaList = cell(P,1);

        WSW = zeros(P,P,N);
        DW = zeros(P,P,N);
        for p = 1:P
            for q = 1:P
                WSW(p,q,:) = dot(w,squeeze(S(p,:,:,q)) * w)';
                DW(p,q,:) = dot(squeeze(D(p,:,q,:)),w);
            end
        end
        JDataVecTemp = ddT - DW - permute(DW,[2,1,3]) + WSW;
        JDataVec = zeros(P,P,N);
        if P <= 1
            JDataVec(:,1,:) = lambda .* squeeze(JDataVecTemp(:,p,:));
        else
            for p = 1:P
                JDataVec(:,p,:) = bsxfun(@times,lambda,squeeze(JDataVecTemp(:,p,:))')';
            end
        end
        iJData = reshape(repmat(1:P,N*P,1)',N*P*P,1) + ...
                                reshape(repmat(P*(0:N-1),P*P,1),N*P*P,1);
        jJData = reshape(repmat(1:N*P,P,1),N*P*P,1);
        JData = sparse(iJData,jJData,JDataVec(:));
        for p = 1:P
            JBetaList{p} = beta(p) * JList{p};
            GaBetaList{p} = sqrt(beta(p)) * GaList{p};
        end
        JTilde = JData(HaInd,HaInd) + blkdiag(JBetaList{:});

        JI = amd(JTilde);
        JTildeP = JTilde(JI,JI);
        M3 = ichol(JTildeP);
        NGa = size(blkdiag(GaBetaList{:}),1);
        RW = zeros(K,P,N);
        JTildeNuMat = zeros(N,P); % JTilde * nu on (N x P)-form
        ar = reshape(a',N*P,1); arP = ar;

        % lambda
        BW  = zeros(P,N);
        WRW = zeros(P,N);
        ddTA = zeros(P,N);
        DWA = zeros(P,N);
        WSWA = zeros(P,N);

    end

    % Gibbs loop
    tic
    for i = 2:(warmup+niter)

        if mod(i,100) == 0
            disp(strcat('Gibbs loop count: ',num2str(i),'/',num2str(warmup+niter),' Time: ',datestr(now,'HH:MM:SS')));
            plotEndi = thinningInd(i-thinningFactor);
            if plotEndi && do_plot
              	afig = figure(98);
                set(afig,'Position',[100 1000 1000 1000]);
                scK = ceil(sqrt(K));
                for kk = 1:K
                    subplot(scK,scK,kk);
                    plot(alphaVec(kk,1:plotEndi)); title(strcat('alpha',num2str(kk)));
                end
                drawnow
                if p > 0
                    bfig = figure(99);
                    set(bfig,'Position',[1200 1000 600 600]);
                    scP = ceil(sqrt(P));
                    for pp = 1:P
                        subplot(scP,scP,pp);
                        plot(betaVec(pp,1:plotEndi)); title(strcat('beta',num2str(pp)));
                    end
                    drawnow
                end
            end;

        end

        if isGLMNoiseAR

            % AR noise

            % w: compute precision
            for j = 1:K
                for k = 1:K
                    ASA(j,k,:) = dot(a,squeeze(S(:,j,k,:)) * a,1)';
                end
                RA(:,j,:) = permute(squeeze(R(:,j,:)) * a,[1,3,2]);
            end
            BDataVecTemp = repmat(XTX,[1,1,N]) -RA - permute(RA,[2,1,3]) + ASA;
            for j = 1:K
                BDataVec(:,j,:) = bsxfun(@times,lambda,squeeze(BDataVecTemp(:,j,:))')';
            end
            BData = sparse(iBData,jBData,BDataVec(:));
            for k = 1:K
                BAlphaList{k} = alpha(k) * BList{k};
                GwAlphaList{k} = sqrt(alpha(k)) * GwList{k};
            end
            BTilde = BData(HwInd,HwInd) + blkdiag(BAlphaList{:});

            % w: compute mean
            if P <= 1
                for j = 1:K
                    DA(1,j,:) = squeeze(D(1,j,:,:)) .* a';
                end
                for j = 1:K
                    BTildeMuMat(:,j) = lambda .* (YTX(:,j) - a' .* squeeze(B(:,j,:) - DA(:,j,:)));
                end
            else
                for p = 1:P
                    for j = 1:K
                        DA(p,j,:) = dot(squeeze(D(p,j,:,:)),a);
                    end
                end
                for j = 1:K
                    BTildeMuMat(:,j) = lambda .* (YTX(:,j) - dot(a,squeeze(B(:,j,:) - DA(:,j,:)))');
                end
            end

            % w: reorder and sample
            BTildeP = BTilde(BI,BI);
            
            if strcmp(samplingMethod,'Cholesky')
                cholBTildeP = chol(BTildeP);
                mu = BTildeMuMat(:);
                muP = mu(BI);
                muP = cholBTildeP' \ muP;
                muP = cholBTildeP \ muP;
                wr =  muP + (cholBTildeP \ randn(K*N,1));
                wr(BI) = wr;
            elseif strcmp(samplingMethod,'PCG') 
                if updatePrecond && (i < 5 || mod(i,50) == 0)
                    try 
                        M1 = ichol(BTildeP);
                    catch
                        disp('ichol error');
                    end
                end
                cholBData = chol(BData);
                bw = blkdiag(GwAlphaList{:})' * randn(NGw,1) + ...
                     cholBData(HwInd,HwInd)' * randn(K*N,1) + BTildeMuMat(:);
                bwP = bw(BI);
                [wrP,flag,relres] = pcg(BTildeP,bwP,PCGTol,500,M1,M1',wrP);
                wr(BI) = wrP;
            end
            w = reshape(wr,N,K)'; % MCMC posterior
            if thinningInd(i);wVec(:,:,thinningInd(i)) = w;end;

            % a: compute precision
            for p = 1:P
                for q = 1:P
                    WSW(p,q,:) = dot(w,squeeze(S(p,:,:,q)) * w)';
                    DW(p,q,:) = dot(squeeze(D(p,:,q,:)),w);
                end
            end
            JDataVecTemp = ddT - DW - permute(DW,[2,1,3]) + WSW;
            if P <= 1
                JDataVec(:,1,:) = lambda .* squeeze(JDataVecTemp(:,p,:));
            else
                for p = 1:P
                    JDataVec(:,p,:) = bsxfun(@times,lambda,squeeze(JDataVecTemp(:,p,:))')';
                end
            end
            JData = sparse(iJData,jJData,JDataVec(:));
            for p = 1:P
                JBetaList{p} = beta(p) * JList{p};
                GaBetaList{p} = sqrt(beta(p)) * GaList{p};
            end
            JTilde = JData(HaInd,HaInd) + blkdiag(JBetaList{:});

            % a: compute mean
            for p = 1:P
                RW(:,p,:) = permute(squeeze(R(:,:,p)) * w,[1,3,2]);
            end
            for p = 1:P
                JTildeNuMat(:,p) = lambda .* (YTdT(:,p) - dot(w,squeeze(B(p,:,:)) - squeeze(RW(:,p,:)))');
            end

            % a: reorder and sample
            JTildeP = JTilde(JI,JI);
            if strcmp(samplingMethod,'Cholesky')
                cholJTildeP = chol(JTildeP);
                nu = JTildeNuMat(:);
                nuP = nu(JI);
                nuP = cholJTildeP' \ nuP;
                nuP = cholJTildeP \ nuP;
                ar =  nuP + (cholJTildeP \ randn(P*N,1));
                ar(JI) = ar;
            elseif strcmp(samplingMethod,'PCG') 
                if updatePrecond && (i < 5 || mod(i,50) == 0)
                    try 
                        M3 = ichol(JTildeP);
                    catch
                        disp('ichol error');
                    end
                end
                cholJData = chol(JData);
                ba = blkdiag(GaBetaList{:})' * randn(NGa,1) + ...
                    cholJData(HaInd,HaInd)' * randn(P*N,1) + JTildeNuMat(:);
                baP = ba(JI);
                [arP,flag,relres] = pcg(JTildeP,baP,PCGTol,500,M3,M3',arP);
                ar(JI) = arP;
            end
            a = reshape(ar,N,P)'; % MCMC posterior
            if thinningInd(i);aVec(:,:,thinningInd(i)) = a;end;

            % lambda: compute u1Tilde
            for p = 1:P
                BW(p,:) = dot(squeeze(B(p,:,:)),w);
                WRW(p,:) = dot(w,squeeze(RW(:,p,:)));
            end
            if P <= 1
                ddTA = squeeze(ddT)' .* a;
                DWA = squeeze(DW)' .* a;
                WSWA = squeeze(WSW)' .* a;
            else
                for p = 1:P
                    ddTA(p,:) = dot(squeeze(ddT(p,:,:)),a);
                    DWA(p,:) = dot(squeeze(DW(p,:,:)),a);
                    WSWA(p,:) = dot(squeeze(WSW(p,:,:)),a);
                end
            end
            u1DataInv = YTY - 2*dot(YTX',w) + dot(w,XTX*w) ...
                        - 2*dot(YTdT',a,1) + 2*dot(BW,a,1) - 2*dot(WRW,a,1) ...
                        + dot(a,ddTA,1) - 2*dot(a,DWA,1) + dot(a,WSWA,1);
            u1Tilde = 1 ./ (1/u1 + .5*u1DataInv)';

            % lambda: sample
            lambda = gamrnd(u2Tilde,u1Tilde);
            if thinningInd(i);lambdaVec(:,thinningInd(i)) = lambda;end;

            % alpha: compute q1Tilde and sample
            for k = 1:K
                q1Tilde(k) = 1 / (1/q1(k) + .5*(w(k,:)*(BList{k}*w(k,:)')));
            end
            alpha = gamrnd(q2Tilde,q1Tilde); % MCMC posterior
            if thinningInd(i);alphaVec(:,thinningInd(i)) = alpha;end;

            % beta: compute r1Tilde and sample
            for p = 1:P
                r1Tilde(p) = 1 / (1/r1(p) + .5*(a(p,:)*(JList{p}*a(p,:)')));
            end
            beta = gamrnd(r2Tilde,r1Tilde); % MCMC posterior
            if thinningInd(i);betaVec(:,thinningInd(i)) = beta;end;

        else

            % White noise        

            % w: compute precision
            BData = kron(spdiags(lambda,0,N,N),XTX);
            for k = 1:K
                BAlphaList{k} = alpha(k) * BList{k};
                GwAlphaList{k} = sqrt(alpha(k)) * GwList{k};
            end
            BTilde = BData(HwInd,HwInd) + blkdiag(BAlphaList{:});

            % w: compute mean
            BTildeMuMat = bsxfun(@times,lambda,YTX);

            % w: reorder and sample
            BTildeP = BTilde(BI,BI);
            
            if strcmp(samplingMethod,'Cholesky')
                cholBTildeP = chol(BTildeP);
                mu = BTildeMuMat(:);
                muP = mu(BI);
                muP = cholBTildeP' \ muP;
                muP = cholBTildeP \ muP;
                wr =  muP + (cholBTildeP \ randn(K*N,1));
                wr(BI) = wr;
            elseif strcmp(samplingMethod,'PCG') 
                if updatePrecond && (i < 5 || mod(i,50) == 0)
                    try 
                        M1 = ichol(BTildeP);
                    catch
                        disp('ichol error');
                    end
                end
                cholBData = chol(BData);
                bw = blkdiag(GwAlphaList{:})' * randn(NGw,1) + ...
                     cholBData(HwInd,HwInd)' * randn(K*N,1) + BTildeMuMat(:);
                bwP = bw(BI);
                [wrP,flag,relres] = pcg(BTildeP,bwP,PCGTol,500,M1,M1',wrP);
                wr(BI) = wrP;
            end
            w = reshape(wr,N,K)'; % MCMC posterior
            if thinningInd(i);wVec(:,:,thinningInd(i)) = w;end;

            % lambda: compute u1Tilde
            u1DataInv = YTY - 2*dot(YTX',w) + dot(w,XTX*w);
            u1Tilde = 1 ./ (1/u1 + .5*u1DataInv)';

            % lambda: sample
            lambda = gamrnd(u2Tilde,u1Tilde);
            if thinningInd(i);lambdaVec(:,thinningInd(i)) = lambda;end;

            % sample alpha
            for k = 1:K
                q1Tilde(k) = 1 / (1/q1(k) + .5*(w(k,:)*(BList{k}*w(k,:)')));
            end
            alpha = gamrnd(q2Tilde,q1Tilde); % MCMC posterior
            if thinningInd(i);alphaVec(:,thinningInd(i)) = alpha;end;
        end

    end
    
    % Compute output summaries

    % Cut away warmup
    MCMC.b(iLb).wVec2 = wVec(:,:,warmup+1:end);
    MCMC.b(iLb).lambdaVec2 = lambdaVec(:,warmup+1:end);
    MCMC.b(iLb).alphaVec2 = alphaVec(:,warmup+1:end);
    nSavedDraws2 = nSavedDraws - warmup;

    % Initiate mean and std storage
    wPostMean = zeros(K,N);   % posterior mean for each voxel
    wPostStd = zeros(K,N);   % posterior std for each voxel
    lambdaPostMean = zeros(1,N);   % posterior mean for each voxel
    alphaPostMean = mean(MCMC.b(iLb).alphaVec2,2);

    % Mean computation w, lambda
    for i = 1:N
        for k = 1:K
            wPostMean(k,i) = mean(squeeze(MCMC.b(iLb).wVec2(k,i,:)));
            wPostStd(k,i) = std(squeeze(MCMC.b(iLb).wVec2(k,i,:)));
        end

        lambdaPostMean(i) = mean(MCMC.b(iLb).lambdaVec2(i,:));
    end

    if isGLMNoiseAR
        % Mean computation a
        MCMC.b(iLb).aVec2 = aVec(:,:,warmup+1:end);
        aPostMean = zeros(P,N);   % posterior mean for each voxel
        for i = 1:N
            for p = 1:P
                aPostMean(p,i) = mean(squeeze(MCMC.b(iLb).aVec2(p,i,:)));
            end
        end
        % beta
        MCMC.b(iLb).betaVec2 = betaVec(:,warmup+1:end);
        betaPostMean = mean(MCMC.b(iLb).betaVec2,2);
    end
    
    MCMC.b(iLb).nSavedDraws2 = nSavedDraws2;
    MCMC.b(iLb).wPostMean = wPostMean;
    MCMC.b(iLb).wPostStd = wPostStd;
    MCMC.b(iLb).lambdaPostMean = lambdaPostMean;
    MCMC.b(iLb).alphaPostMean = alphaPostMean;
    if isGLMNoiseAR
        MCMC.b(iLb).aPostMean = aPostMean;
        MCMC.b(iLb).betaPostMean = betaPostMean;
    end
    MCMC.b(iLb).niter = niter;
    MCMC.b(iLb).warmup = warmup;
    MCMC.b(iLb).thinningFactor = thinningFactor;
    MCMC.b(iLb).timeTaken = toc;
    
end

save(strcat(outputPath,'sub-',subjStr,'/',estimationMethod,'/MCMC.mat'),'MCMC','-v7.3');
disp('done.');

