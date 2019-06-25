classdef FMRIModel < handle
  %FMRIModel Model for fMRI analysis
  %
  %
  % AUTHOR:       Per Siden
  %               Division of Statistics and Machine Learning
  %               Department of Computer and Information Science
  %               Linkoping University
  %
  % FIRST VER.:   2017-09-05
  
  properties
    
    fS                  % fMRI analysis settings
    priorList           % List of priors
    nM                  % Noise model
    QkList              % List of prior precision matrices
    Q                   % Prior precision matrix
    X                   % TxK Design matrix
    Y                   % TxN Data matrix
    gSF                 % global scaling factor
    maskInd             % Brain mask indices
    N                   % Number of voxels
    priorNameList       % List of names of priors
    marglik             % Marginal likelihood
    marglikVec          % Results storage
    PNK, PKN            % Reorderings
    iQData, jQData      % For blockdiagonal matrices
    indQData
    XTX, YTX, YTY       % Precalculations
    XTilde, d, R, S
    B, D, ddT, YTdT
    YTXTilde, dX, dXTilde
    QnTilde             % Data part of the posterior precision matrix
    QnTildeLambda
    b, bReo
    QTilde, QTildeReo
    reoQ, ireoQ
    iQTilde, icholQTildeReo, RQTildeReo
    muTilde, muTildeReo, M
    wCovMat
    Vs, iQTildeVs
    randPool
    startIcholEps1, startIcholEps2
    iter
    hessAlpha, gam2
    
  end
  
  methods
    
    function fM = FMRIModel(fS)
      fM.fS = fS;
    end
    function doPrecalculations(fM)
      
      K = fM.fS.K; N = fM.N; T = fM.fS.T; P = fM.fS.P;
      
      % Dimension reordering
      fM.PNK = reshape(1:(K*N),K,N)';fM.PNK = fM.PNK(:);fM.PKN(fM.PNK) = 1:N*K;
      
      % if isGLMNoiseAR;PNP = reshape(1:(P*N),P,N)';PNP = PNP(:);PPN(PNP) = 1:N*P;end; % not used
      fM.iQData = reshape(repmat(1:K,N*K,1)',N*K*K,1) + ...
        reshape(repmat(K*(0:N-1),K*K,1),N*K*K,1);
      fM.jQData = reshape(repmat(1:N*K,K,1),N*K*K,1);
      fM.indQData = sub2ind([N*K,N*K],fM.iQData,fM.jQData);      
      
      if ~fM.fS.isGLMNoiseAR
        fM.XTX = fM.X'*fM.X; fM.YTY = dot(fM.Y,fM.Y); fM.YTX = fM.Y'*fM.X;
      else
        
        % lagged X
        fM.XTilde = zeros(P,T-P,K);
        for t = 1:(T-P)
          fM.XTilde(:,t,:) = permute(fM.X((t+P-1):-1:t,:),[1,3,2]);
        end
        Xcut = fM.X((P+1):T,:);
        
        % lagged Y
        fM.d = zeros(P,T-P,N);
        for t = 1:(T-P)
          fM.d(:,t,:) = permute(fM.Y((t+P-1):-1:t,:),[1,3,2]);
        end
        Ycut = fM.Y((P+1):T,:);
        
        % voxel independent
        fM.XTX = Xcut'*Xcut;
        fM.R = zeros(K,K,P);
        fM.S = zeros(P,K,K,P);
        for j = 1:K
          fM.R(:,j,:) = permute(Xcut' * squeeze(fM.XTilde(:,:,j))',[1,3,2]);
          for k = 1:K
            fM.S(:,j,k,:) = permute(squeeze(fM.XTilde(:,:,j)) * ...
              squeeze(fM.XTilde(:,:,k))',[3,1,2,4]);
          end
        end
        
        % voxel dependent
        fM.YTX = Ycut'*Xcut; fM.YTY = dot(Ycut,Ycut);
        fM.B = zeros(P,K,N); % B = Gxy + Rxy in Penny (2005) (vb3)
        fM.D = zeros(P,K,P,N); fM.ddT = zeros(P,P,N); fM.YTdT = zeros(N,P);
        
        for n = 1:N
          fM.YTXTilde = zeros(P,K);
          fM.dX = zeros(P,K);
          fM.dXTilde = zeros(P,K,P);
          for p = 1:P
            fM.YTXTilde(p,:) = Ycut(:,n)' * squeeze(fM.XTilde(p,:,:));
            fM.dX(p,:) = squeeze(fM.d(p,:,n)) * Xcut;
            fM.dXTilde(:,:,p) = squeeze(fM.d(:,:,n)) * squeeze(fM.XTilde(p,:,:));
          end
          fM.B(:,:,n) = fM.YTXTilde + fM.dX;
          fM.D(:,:,:,n) = fM.dXTilde;
          fM.ddT(:,:,n) = squeeze(fM.d(:,:,n)) * squeeze(fM.d(:,:,n))';
          fM.YTdT(n,:) = Ycut(:,n)' * squeeze(fM.d(:,:,n))';
        end
      end
    end
    function initialize(fM)
      K = fM.fS.K; N = fM.N; KSpat = fM.fS.KSpat;
      
      fM.iter = 0;
      fM.marglik = 0;
      fM.marglikVec = zeros(fM.fS.maxiter+1,1);
      fM.QkList = cell(fM.fS.K,1);
      fM.computeQTilde;
      fM.reoQ = amd(fM.QTilde);fM.ireoQ(fM.reoQ) = 1:(N*K);
      fM.QTildeReo = fM.QTilde(fM.reoQ,fM.reoQ);   
      fM.muTilde = zeros(N*K,1);fM.muTildeReo = zeros(N*K,1);
      
      if strcmp(fM.fS.solveMethod,'PCG')
        fM.startIcholEps1 = 1e-15;fM.startIcholEps2 = 1e-15;
        if fM.fS.useHutchRandomPool
            fM.randPool = HutchRandomPool(fM.N*fM.fS.K,fM.fS.SHutchRandomPool,fM.fS.Ns);
        end
      end            
      
    end
    function gradientStep(fM)
      
      if fM.fS.doPrint; fM.dispStatus; end;
      if fM.fS.doPlot; fM.plotConvergence; end;
      K = fM.fS.K; N = fM.N; T = fM.fS.T; P = fM.fS.P; KSpat = fM.fS.KSpat;
      
      fM.iter = fM.iter + 1;
      fM.computeQTilde;
      fM.bReo = fM.b(fM.reoQ);
      fM.QTildeReo = fM.QTilde(fM.reoQ,fM.reoQ);
      
      % Compute muTilde and traces needed for the gradients      
      if strcmp(fM.fS.solveMethod,'Chol') || fM.fS.doMargLikComp
        fM.RQTildeReo = chol(fM.QTildeReo);
        iQTildeReo = Qinv(fM.RQTildeReo);
        fM.iQTilde = iQTildeReo(fM.ireoQ,fM.ireoQ);
        fM.muTildeReo = fM.RQTildeReo \ (fM.bReo' / fM.RQTildeReo)';
      else
        
        if strcmp(fM.fS.solveMethod,'PCG')
          if fM.fS.updatePrecond && (fM.iter < 5 || ...
              mod(fM.iter,fM.fS.precondUpdateFreq) == 0)
            try [fM.icholQTildeReo,fM.startIcholEps1] = ...
                icholSafe(fM.QTildeReo,fM.startIcholEps1,fM.fS.maxIcholEps);
            catch
              fM.icholQTildeReo = spdiags(1./sqrt(diag(fM.QTildeReo)),0,N*K,N*K);
            end
          end
%           tic
          [fM.muTildeReo,flag,relres] = pcg(fM.QTildeReo,fM.bReo,fM.fS.PCGtol,...
              fM.fS.PCGmaxiter,fM.icholQTildeReo,fM.icholQTildeReo',fM.muTildeReo);
          if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for mean.']); end

        elseif strcmp(fM.fS.solveMethod,'MG')
          
        end
        if fM.fS.useHutchRandomPool
          [fM.Vs,lastiQTildeVs] = fM.randPool.getSamples;
        else
          fM.Vs = 2*(round(rand(N*K,fM.fS.Ns))-.5);
          lastiQTildeVs = zeros(N*K,fM.fS.Ns);
        end
        VsReo = fM.Vs(fM.reoQ,:);
        fM.iQTildeVs = zeros(N*K,fM.fS.Ns);
        if strcmp(fM.fS.solveMethod,'PCG')
          if fM.fS.doParallel
            samps = parallelPCG(fM.QTildeReo,VsReo,fM.fS.PCGtol,...
              fM.fS.PCGmaxiter,fM.icholQTildeReo,lastiQTildeVs(fM.reoQ,:));
            fM.iQTildeVs = samps(fM.ireoQ,:);
          else
            for i = 1:fM.fS.Ns
              [samp,flag,relres] = pcg(fM.QTildeReo,VsReo(:,i),fM.fS.PCGtol,...
                  fM.fS.PCGmaxiter,fM.icholQTildeReo,fM.icholQTildeReo',...
                  lastiQTildeVs(fM.reoQ,i));                
              if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for QTilde trace.']); end
              fM.iQTildeVs(:,i) = samp(fM.ireoQ);
            end
          end
        elseif strcmp(fM.fS.solveMethod,'MG')
          
        end
        if fM.fS.useHutchRandomPool;fM.randPool.storeSolutions(fM.iQTildeVs);end
        
      end
      
      fM.muTilde = fM.muTildeReo(fM.ireoQ);
      fM.M = reshape(fM.muTilde,[N,K])';

      if fM.fS.useRobbinsMonro
        if fM.iter >= fM.fS.startIterRobbinsMonro
          fM.hessAlpha = fM.fS.hessAlphaStart/(.1*(fM.iter-fM.fS.startIterRobbinsMonro)+1);
          fM.gam2 = 1;
        else
          fM.hessAlpha = fM.fS.hessAlphaStart;   
          fM.gam2 = fM.fS.gam2Start;        
        end
      else
        fM.hessAlpha = fM.fS.hessAlphaStart;   
        fM.gam2 = fM.fS.gam2Start;      
      end
      for k = 1:K
        fM.priorList{k}.gradientStep(fM);
      end
      fM.nM.gradientStep(fM);      
      if fM.fS.doMargLikComp; fM.computeMargLik(fM);end;
      
    end
    function doPolyakAveraging(fM)
      K = fM.fS.K;
      for k = 1:K
        fM.priorList{k}.doPolyakAveraging(fM);
      end
      fM.nM.doPolyakAveraging(fM);      
    end
    function computePostMeanAndStd(fM)
      K = fM.fS.K; N = fM.N; T = fM.fS.T; P = fM.fS.P;
      fM.computeQTilde;
      fM.bReo = fM.b(fM.reoQ);
      fM.QTildeReo = fM.QTilde(fM.reoQ,fM.reoQ);
      
      if strcmp(fM.fS.solveMethod,'Chol')
        fM.RQTildeReo = chol(fM.QTildeReo);
        iQTildeReo = Qinv(fM.RQTildeReo);
        fM.iQTilde = iQTildeReo(fM.ireoQ,fM.ireoQ);
        fM.muTildeReo = fM.RQTildeReo \ (fM.bReo' / fM.RQTildeReo)';
        iQTildeKN = fM.iQTilde(fM.PKN,fM.PKN);
        fM.wCovMat = sparse(fM.iQData,fM.jQData,iQTildeKN(fM.indQData));
        
      elseif strcmp(fM.fS.solveMethod,'PCG')        
        if fM.fS.updatePrecond && (fM.iter < 5 || ...
            mod(fM.iter,fM.fS.precondUpdateFreq) == 0)
          try [fM.icholQTildeReo,fM.startIcholEps1] = ...
              icholSafe(fM.QTildeReo,fM.startIcholEps1,fM.fS.maxIcholEps);
          catch
            fM.icholQTildeReo = spdiags(1./sqrt(diag(fM.QTildeReo)),0,N*K,N*K);
          end
        end
        [fM.muTildeReo,flag,relres] = pcg(fM.QTildeReo,fM.bReo,fM.fS.PCGtol,...
            fM.fS.PCGmaxiter,fM.icholQTildeReo,fM.icholQTildeReo',fM.muTildeReo);
        if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for mean.']); end
        
        % Estimation by simple RBMC
        try [M1,fM.startIcholEps1] = icholSafe(fM.QTildeReo,fM.startIcholEps1,fM.fS.maxIcholEps);
        catch; M1 = spdiags(1./sqrt(diag(fM.QTildeReo)),0,N*K,N*K); end;        
        inversePriorSamples = zeros(N*K,fM.fS.NsRBMC);
        for k = 1:K
          inversePriorSamples((k-1)*N+1:k*N,:) = fM.priorList{k}.getInversePriorSamples(fM.fS.NsRBMC);
        end          
       
        QDataLambda = sparse(fM.iQData,fM.jQData,fM.QnTildeLambda(:));
        cholQnTildeLambda = chol(QDataLambda);
        bSamp = inversePriorSamples + ...
          cholQnTildeLambda(fM.PNK,fM.PNK)' * randn(K*N,fM.fS.NsRBMC);
        bSamp = bSamp(fM.reoQ,:);
        demeanSampMat = zeros(K*N,fM.fS.NsRBMC);
        if fM.fS.doParallel
          demeanSamps = parallelPCG(fM.QTildeReo,bSamp,fM.fS.PCGtol,...
            fM.fS.PCGmaxiter,M1,zeros(N*K,fM.fS.NsRBMC));
          demeanSampMat = demeanSamps(fM.ireoQ,:);
        else
          for i = 1:fM.fS.NsRBMC
            [demeanSamp,flag,relres] = pcg(fM.QTildeReo,bSamp(:,i),...
                          fM.fS.PCGtol,fM.fS.PCGmaxiter,M1,M1',zeros(N*K,1));                        
            if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for RBMC.']); end
            demeanSampMat(:,i) = demeanSamp(fM.ireoQ);
          end
        end
        
        QTildeKN = fM.QTilde(fM.PKN,fM.PKN);
        QmDQ = spdiags(zeros(N*K,2*K-1),(-K+1):(K-1),QTildeKN);
        QmDQX = QmDQ * demeanSampMat(fM.PKN,:);
        DQ = QTildeKN - QmDQ;
        invCholDQ = chol(DQ) \ speye(N*K);
        invDQ = invCholDQ * invCholDQ';
        xi = invDQ * QmDQX;
        xi2 = permute(reshape(xi,K,N,fM.fS.NsRBMC),[2 3 1]);
        xi3 = repmat(reshape(xi2,N*fM.fS.NsRBMC,K),[1,1,K]);
        xi4 = reshape(xi3 .* permute(xi3,[1 3 2]),[N,fM.fS.NsRBMC,K,K]);
        xi5 = permute(squeeze(sum(xi4,2)),[2,3,1]);
        [i2,j2] = ind2sub([N*K,N*K],find(kron(speye(N),ones(K))==1));
        fM.wCovMat = invDQ +1/fM.fS.NsRBMC*sparse(i2,j2,xi5(:),N*K,N*K);
      end
      
      fM.muTilde = fM.muTildeReo(fM.ireoQ);
      fM.M = reshape(fM.muTilde,[N,K])';
      
    end
    function computeQTilde(fM)
      K = fM.fS.K; N = fM.N; KSpat = fM.fS.KSpat;
      
      for k = 1:K
        fM.QkList{k} = fM.priorList{k}.computeQk;
      end
      fM.Q = blkdiag(fM.QkList{:});
      [fM.QnTilde,fM.QnTildeLambda,fM.b] = fM.nM.computeQnTilde(fM);
      QDataLambda = sparse(fM.iQData,fM.jQData,fM.QnTildeLambda(:));
      fM.QTilde = QDataLambda(fM.PNK,fM.PNK) + fM.Q;
            
    end   
    function computeMargLik(fM) % Works only for iid?
      K = fM.fS.K;
      % Marginal likelihood
      RQkReo = cell(K,1); logDetQk = zersos(K,1);
      for k = 1:K
        reoQk = amd(fM.QkList{k});
        RQkReo{k} = chol(fM.QkList{k}(reoQk{k},reoQk{k}));
        logDetQk(k) = 2*sum(log(diag(RQkReo{k})));
      end
      fM.marglik = .5*(fM.fS.T-fM.fS.P)*sum(log(fM.nM.lambda)) - ...
                .5*fM.nM.lambda'*fM.nM.quadFormLambda1 + ...
                .5*sum(logDetQk) - .5*fM.muTilde'*fM.Q*fM.muTilde - ...
                sum(log(diag(fM.RQTildeReo)));
      fM.marglikVec(fM.iter) = fM.marglik;
    end
    function dispStatus(fM)
      dispVec = ['iter: ',num2str(fM.iter),', ',fM.priorList{1}.getParameterString,...
                 ', ',fM.nM.getParameterString,',',char(10),char(9),'logMargLike: ',...
                 num2str(fM.marglik),', Time: ',datestr(now,'dd-mmm-yyyy HH:MM:SS')];
      disp(dispVec);
    end
    function plotConvergence(fM)
      if fM.fS.doPlot && (mod(fM.iter,fM.fS.plotUpdateFreq) == 1 || fM.iter == fM.fS.maxiter)
        figure(10)
        plotInd = max(1,fM.iter-fM.fS.plotLength+1):fM.iter;
        nbrOfVoxelsToPlot = 100; stepMNOVTP = ceil(fM.N / nbrOfVoxelsToPlot);
        tau2Vec = fM.priorList{1}.tau2Vec'; kappa2Vec = fM.priorList{1}.kappa2Vec';
        for k = 2:fM.fS.K
          if k == 3 || k == 5 || k == 7
            tau2Vec = [tau2Vec;fM.priorList{k}.tau2Vec'];
            kappa2Vec = [kappa2Vec;fM.priorList{k}.kappa2Vec'];
          end
        end
        subplot(321); semilogy(plotInd,tau2Vec(:,plotInd)); title('tau2');
        subplot(322); semilogy(plotInd,kappa2Vec(:,plotInd)); title('kappa2');
        subplot(323); plot(plotInd,fM.nM.lambdaVec(1:stepMNOVTP:fM.N,plotInd)'); title('lambda');
        subplot(324); plot(plotInd,fM.marglikVec(plotInd)); title('log marglik');
        if fM.fS.isGLMNoiseAR
          subplot(325); plot(plotInd,squeeze(fM.nM.aVec(1,1:stepMNOVTP:fM.N,plotInd))'); title('a');
          subplot(326); hist(fM.nM.a(:)); title('a');
        end
        if strcmp(fM.priorList{1}.priorName,'Matern2Aniso')
          hxVec = fM.priorList{1}.hxVec'; hyVec = fM.priorList{1}.hyVec';
          for k = 2:fM.fS.K
            if k == 3 || k == 5 || k == 7              
              hxVec = [hxVec;fM.priorList{k}.hxVec'];
              hyVec = [hyVec;fM.priorList{k}.hyVec'];
            end
          end
          subplot(325); semilogy(plotInd,hxVec(:,plotInd)); title('hx');
          subplot(326); semilogy(plotInd,hyVec(:,plotInd)); title('hy');
        end
        drawnow
      end
    end
    function fMThin = saveModelThin(fM)
      fMThin.fS = fM.fS; 
      fMThin.priorList = cell(fM.fS.K,1);
      for k = 1:fMThin.fS.K
        fMThin.priorList{k} = {fM.priorList{k}.priorName,fM.priorList{k}.getParameters};
      end
      fMThin.lambda = fM.nM.lambda;
      if fMThin.fS.isGLMNoiseAR; fMThin.a = fM.nM.a; end
      
    end 
    function [MOOS,wCovMatOOS] = predictOutOfSample(fM,lostInd,computeVariance)
      K = fM.fS.K; N = fM.N; KSpat = fM.fS.KSpat;
      if size(lostInd,1) < size(lostInd,2); lostInd = lostInd'; end;
      lostIndLong = lostInd + N*(0:K-1);
      lostIndLong = lostIndLong(:); 
      
      for k = 1:K
        fM.QkList{k} = fM.priorList{k}.computeQk;
      end
      fM.Q = blkdiag(fM.QkList{:});
      [fM.QnTilde,QnTildeLambda,b] = fM.nM.computeQnTilde(fM);
      QnTildeLambda(:,:,lostInd) = 0;
      b(lostIndLong) = 0;
      QDataLambda = sparse(fM.iQData,fM.jQData,QnTildeLambda(:));
      QTilde = QDataLambda(fM.PNK,fM.PNK) + fM.Q;
      
      bReo = b(fM.reoQ);
      QTildeReo = QTilde(fM.reoQ,fM.reoQ);
      
      if strcmp(fM.fS.solveMethod,'Chol')
        RQTildeReo = chol(QTildeReo);
        muTildeReo = RQTildeReo \ (bReo' / RQTildeReo)';    
        if computeVariance
          iQTildeReo = Qinv(RQTildeReo');
          iQTilde = iQTildeReo(fM.ireoQ,fM.ireoQ);
          iQTildeKN = iQTilde(fM.PKN,fM.PKN);
          wCovMatOOS = sparse(fM.iQData,fM.jQData,iQTildeKN(fM.indQData));
        end

      elseif strcmp(fM.fS.solveMethod,'PCG')        
        if fM.fS.updatePrecond && (fM.iter < 5 || ...
            mod(fM.iter,fM.fS.precondUpdateFreq) == 0)
          try [fM.icholQTildeReo,fM.startIcholEps1] = ...
              icholSafe(QTildeReo,fM.startIcholEps1,fM.fS.maxIcholEps);
          catch
            fM.icholQTildeReo = spdiags(1./sqrt(diag(QTildeReo)),0,N*K,N*K);
          end
        end
        [muTildeReo,flag,relres] = pcg(QTildeReo,bReo,fM.fS.PCGtol,...
            fM.fS.PCGmaxiter,fM.icholQTildeReo,fM.icholQTildeReo',fM.muTildeReo);
        if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for mean.']); end
        
        if computeVariance
                    
          % Estimation by simple RBMC 
          try [M1,fM.startIcholEps1] = icholSafe(QTildeReo,fM.startIcholEps1,fM.fS.maxIcholEps);
          catch; M1 = spdiags(1./sqrt(diag(QTildeReo)),0,N*K,N*K); end;
          inversePriorSamples = zeros(N*K,fM.fS.NsRBMC);
          for k = 1:K
            inversePriorSamples((k-1)*N+1:k*N,:) = fM.priorList{k}.getInversePriorSamples(fM.fS.NsRBMC);
          end          

          lostIndLongKN = K*(lostInd-1)+(1:K);
          lostIndLongKN = lostIndLongKN(:)';  
          diagEps = sparse(lostIndLongKN,lostIndLongKN,1e-40*ones(length(lostIndLongKN),1),N*K,N*K);
          QDataLambdaTemp = QDataLambda + diagEps;
          cholQnTildeLambda = chol(QDataLambdaTemp);
          bSamp = inversePriorSamples + ...
            cholQnTildeLambda(fM.PNK,fM.PNK)' * randn(K*N,fM.fS.NsRBMC);
          bSamp = bSamp(fM.reoQ,:);
          demeanSampMat = zeros(K*N,fM.fS.NsRBMC);
          if fM.fS.doParallel
            demeanSamps = parallelPCG(QTildeReo,bSamp,fM.fS.PCGtol,...
              fM.fS.PCGmaxiter,M1,zeros(N*K,fM.fS.NsRBMC));
            demeanSampMat = demeanSamps(fM.ireoQ,:);
          else
            for i = 1:fM.fS.NsRBMC
              [demeanSamp,flag,relres] = pcg(QTildeReo,bSamp(:,i),...
                            fM.fS.PCGtol,fM.fS.PCGmaxiter,M1,M1',zeros(N*K,1));                        
              if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for RBMC.']); end
              demeanSampMat(:,i) = demeanSamp(fM.ireoQ);
            end
          end

          QTildeKN = QTilde(fM.PKN,fM.PKN);
          QmDQ = spdiags(zeros(N*K,2*K-1),(-K+1):(K-1),QTildeKN);
          QmDQX = QmDQ * demeanSampMat(fM.PKN,:);
          DQ = QTildeKN - QmDQ;
          invCholDQ = chol(DQ) \ speye(N*K);
          invDQ = invCholDQ * invCholDQ';
          xi = invDQ * QmDQX;
          xi2 = permute(reshape(xi,K,N,fM.fS.NsRBMC),[2 3 1]);
          xi3 = repmat(reshape(xi2,N*fM.fS.NsRBMC,K),[1,1,K]);
          xi4 = reshape(xi3 .* permute(xi3,[1 3 2]),[N,fM.fS.NsRBMC,K,K]);
          xi5 = permute(squeeze(sum(xi4,2)),[2,3,1]);
          [i2,j2] = ind2sub([N*K,N*K],find(kron(speye(N),ones(K))==1));
          wCovMatOOS = invDQ +1/fM.fS.NsRBMC*sparse(i2,j2,xi5(:),N*K,N*K);
        end
      end      
      
      muTilde = muTildeReo(fM.ireoQ);
      MOOS = reshape(muTilde,[N,K])';
      
      if ~computeVariance
        wCovMatOOS = 0;
      end
      
    end
    function sampMat = sampleFromHrfRegPrior(fM,Ns,KMax)

      N = fM.N; KSpat = fM.fS.KSpat;
      
      for k = 1:KMax
        QkList{k} = fM.priorList{k}.computeQk;
      end
      Q = blkdiag(QkList{:});
      reoQ = amd(Q); ireoQ(reoQ) = 1:(N*KMax);
      QReo = Q(reoQ,reoQ);
      
      try [M1,fM.startIcholEps1] = icholSafe(QReo,fM.startIcholEps1,fM.fS.maxIcholEps);
      catch; M1 = spdiags(1./sqrt(diag(QReo)),0,N*KMax,N*KMax); end       
      inversePriorSamples = zeros(N*KMax,Ns);
      for k = 1:KMax
        inversePriorSamples((k-1)*N+1:k*N,:) = fM.priorList{k}.getInversePriorSamples(Ns);
      end      
            
      bSamp = inversePriorSamples;
      bSamp = bSamp(reoQ,:);
      sampMat = zeros(KMax*N,Ns);
      if fM.fS.doParallel
        demeanSamps = parallelPCG(QReo,bSamp,fM.fS.PCGtol,...
          fM.fS.PCGmaxiter,M1,zeros(N*KMax,Ns));
        sampMat = demeanSamps(ireoQ,:);
      else
        for i = 1:Ns
          [samp,flag,relres] = pcg(QReo,bSamp(:,i),...
                        fM.fS.PCGtol,10*fM.fS.PCGmaxiter,M1,M1',zeros(N*KMax,1));                        
          if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for prior sampling.']); end
          sampMat(:,i) = samp(ireoQ);
        end
      end
    end
  end  
end