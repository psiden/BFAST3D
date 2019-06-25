classdef NoiseModel < handle
  %NoiseModel fMRI Noise Model
  %
  %
  % AUTHOR:       Per Siden
  %               Division of Statistics and Machine Learning
  %               Department of Computer and Information Science
  %               Linkoping University
  %
  % FIRST VER.:   2017-09-07
  
  properties
    
    P                   % AR order
    isGLMNoiseAR
    lambda              % Name of the prior
    a                   % PxN matrix of AR coefficients
    u1, u2              % lambda gamma prior hyperparameters
    tauA2               % AR parameter normal prior hyperparameter
    lambdaHB, aHB       % Hard bounds
    lambdaVec, aVec     % Results storage
    SA, ASA, RA, DA
    QTildeMuMat, dQTildedAKN
    quadFormLambda1
    traceA1, quadFormA1
    BM, MRM, ddTA, DMA, MDA
    MSMA, MSM, DM, RM
    dlambda0Avg, d2lambda0Avg, deltaLambda0
    unstableInd
  end
  
  methods
    
    function ob = NoiseModel(fM)
      K = fM.fS.K; N = fM.N;
      ob.P = fM.fS.P; P = ob.P;
      ob.isGLMNoiseAR = fM.fS.isGLMNoiseAR;
      
      % Set priors
      ob.u1 = 10; ob.u2 = 0.1;
      ob.tauA2 = 1e-3;
      
      % Hard bounds
      ob.lambdaHB = [1e-6,1e6];
      ob.aHB = (1-1e-4) * [-1,1];
      
      % Set start values
      betaOLS = fM.X \ fM.Y;
      if ~fM.fS.isGLMNoiseAR
        sigma2OLS = mean((fM.Y-fM.X*betaOLS).^2)';
      else
        sigma2OLS = mean((fM.Y-fM.X*betaOLS).^2)';
        ob.a = zeros(fM.fS.P,fM.N);
      end
      if fM.fS.preOptNoiseParams
        ob.initializeParameters(fM);
      end
      
      % Results storage
      ob.lambda = 1./sigma2OLS; % 0.5*ones(N,1); %
      ob.lambdaVec = zeros(fM.N,fM.fS.maxiter+1);ob.lambdaVec(:,1) = ob.lambda;
      if fM.fS.isGLMNoiseAR
        ob.aVec = zeros(fM.fS.P,fM.N,fM.fS.maxiter+1);ob.aVec(:,:,1) = ob.a;
      end
      
      % Allocate
      if fM.fS.isGLMNoiseAR
        ob.SA = zeros(P,K,K,N); ob.ASA = zeros(K,K,N);
        ob.RA = zeros(K,K,N); ob.DA = zeros(P,K,N);
        ob.QTildeMuMat = zeros(N,K); ob.dQTildedAKN = cell(P,1);
        ob.traceA1 = zeros(P,N); ob.quadFormA1 = zeros(P,N);
        ob.BM  = zeros(P,N); ob.MRM = zeros(P,N); ob.ddTA = zeros(P,N);
        ob.DMA = zeros(P,N); ob.MDA = zeros(P,N); ob.MSMA = zeros(P,N);
        ob.MSM = zeros(P,P,N); ob.DM = zeros(P,P,N); ob.RM = zeros(K,P,N);
        if fM.fS.checkARStability
          ob.unstableInd = zeros(N,1);
        end
      end
    end
    function [QnTilde,QnTildeLambda,b] = computeQnTilde(ob,fM)%,XTX,R,S,fS,N)
      K = fM.fS.K; P = ob.P; N = fM.N;
      if ~ob.isGLMNoiseAR
        QnTilde = repmat(fM.XTX,[1,1,N]);
        b = reshape(ob.lambda.*fM.YTX,[N*K,1]);
      else
        for p = 1:P
          for j = 1:K
            ob.SA(p,j,:,:) = (squeeze(fM.S(p,j,:,:)) * ob.a);
          end
        end
        %%% Make computation of ASA simpler using SA at some point
        for j = 1:K
          for k = 1:K; ob.ASA(j,k,:) = dot(ob.a,squeeze(fM.S(:,j,k,:)) * ob.a,1)'; end
          ob.RA(:,j,:) = permute(squeeze(fM.R(:,j,:)) * ob.a,[1,3,2]);
        end
        QnTilde = repmat(fM.XTX,[1,1,N]) - ob.RA - ...
          permute(ob.RA,[2,1,3]) + ob.ASA;
        
      end
      QnTildeLambda = zeros(K,K,N);
      for j = 1:K; QnTildeLambda(:,j,:) = ...
          (ob.lambda .* squeeze(QnTilde(:,j,:))')'; end
      if ob.isGLMNoiseAR
        if P <= 1
          for j = 1:K; ob.DA(1,j,:) = squeeze(fM.D(1,j,:,:)) .* ob.a'; end
          for j = 1:K; ob.QTildeMuMat(:,j) = ob.lambda .* ...
              (fM.YTX(:,j) - ob.a' .* squeeze(fM.B(:,j,:) - ob.DA(:,j,:))); end
        else
          for p = 1:P
            for j = 1:K; ob.DA(p,j,:) = dot(squeeze(fM.D(p,j,:,:)),ob.a); end
          end
          for j = 1:K; ob.QTildeMuMat(:,j) = ob.lambda .* ...
              (fM.YTX(:,j) - dot(ob.a,squeeze(fM.B(:,j,:) - ob.DA(:,j,:)))'); end
        end
        b = reshape(ob.QTildeMuMat,[N*K,1]);
      end
    end
    function gradientStep(ob,fM)
      N = fM.N; P = ob.P; K = fM.fS.K;
      
      if fM.fS.isGLMNoiseAR
        
        % a:
        dQnTildedAVec = -repmat(permute(fM.R,[3,1,2]) + permute(fM.R,[3,2,1]),[1,1,1,N]) + ...
          ob.SA + permute(ob.SA,[1,3,2,4]);
        dQnTildedALambdaVec = permute(ob.lambda,[4,3,2,1]) .* dQnTildedAVec;
        for p = 1:P
          dQtemp = squeeze(dQnTildedALambdaVec(p,:,:,:));
          ob.dQTildedAKN{p} = sparse(fM.iQData,fM.jQData,dQtemp(:));
        end
      end
      
      if strcmp(fM.fS.solveMethod,'Chol')
        iQTildeKN = fM.iQTilde(fM.PKN,fM.PKN);
        traceLambda1 = sum(reshape(iQTildeKN(fM.indQData).*...
          fM.QnTilde(:),[K*K,N]))';
        
        if fM.fS.isGLMNoiseAR
          for p = 1:P
            ob.traceA1(p,:) = sum(reshape(iQTildeKN(fM.indQData).* ...
              ob.dQTildedAKN{p}(fM.indQData),[K*K,N]));
          end
        end
      elseif strcmp(fM.fS.solveMethod,'PCG')
        VsKN = fM.Vs(fM.PKN,:);
        QData = sparse(fM.iQData,fM.jQData,fM.QnTilde(:));
        traceLambda1 = 1/fM.fS.Ns*sum(reshape(dot(fM.iQTildeVs(fM.PKN,:),...
          QData*VsKN,2),[K,N]))';
        if fM.fS.isGLMNoiseAR
          for p = 1:P
            ob.traceA1(p,:) = 1/fM.fS.Ns*sum(reshape(dot(fM.iQTildeVs(fM.PKN,:),...
              ob.dQTildedAKN{p}*VsKN,2),[K,N]));
          end
        end
      end
      
      % lambda and a: compute quadratic form
      if ~fM.fS.isGLMNoiseAR
        
        ob.quadFormLambda1 = (fM.YTY - 2*dot(fM.YTX',fM.M) + dot(fM.M,fM.XTX*fM.M))';
        
      else
        
        for p = 1:P
          for q = 1:P
            ob.MSM(p,q,:) = dot(fM.M,squeeze(fM.S(p,:,:,q)) * fM.M)';
            ob.DM(p,q,:) = dot(squeeze(fM.D(p,:,q,:)),fM.M);
          end
        end
        
        for p = 1:P
          ob.RM(:,p,:) = permute(squeeze(fM.R(:,:,p)) * fM.M,[1,3,2]);
        end
        for p = 1:P
          ob.BM(p,:) = dot(squeeze(fM.B(p,:,:)),fM.M);
          ob.MRM(p,:) = dot(fM.M,squeeze(ob.RM(:,p,:)));
        end
        if P <= 1
          ob.ddTA = squeeze(fM.ddT)' .* ob.a;
          ob.DMA = squeeze(ob.DM)' .* ob.a;
          ob.MDA = squeeze(ob.DM)' .* ob.a;
          ob.MSMA = squeeze(ob.MSM)' .* ob.a;
        else
          for p = 1:P
            ob.ddTA(p,:) = dot(squeeze(fM.ddT(p,:,:)),ob.a);
            ob.DMA(p,:) = dot(squeeze(ob.DM(p,:,:)),ob.a);
            ob.MDA(p,:) = dot(squeeze(ob.DM(:,p,:)),ob.a);
            ob.MSMA(p,:) = dot(squeeze(ob.MSM(p,:,:)),ob.a);
          end
        end
        ob.quadFormLambda1 = (fM.YTY - 2*dot(fM.YTX',fM.M) + dot(fM.M,fM.XTX*fM.M) ...
          - 2*dot(fM.YTdT',ob.a,1) + 2*dot(ob.BM,ob.a,1) - 2*dot(ob.MRM,ob.a,1) ...
          + dot(ob.a,ob.ddTA,1) - 2*dot(ob.a,ob.DMA,1) + dot(ob.a,ob.MSMA,1))';
        
        % a:
        ob.quadFormA1 = ob.lambda' .* (2*(-fM.YTdT' + ob.BM - ob.MRM + ...
          ob.ddTA - ob.DMA - ob.MDA + ob.MSMA));
        
      end
      
      % prior contributions
      dPriorLambda = (ob.u2-1)./ob.lambda - 1/ob.u1;
      if fM.fS.useHessNoise; d2PriorLambda = 0; end
      if fM.fS.isGLMNoiseAR
        dPriorA = -ob.tauA2*ob.a;
      end
      
      % Gradient
      dlambda = .5*((fM.fS.T-P)./ob.lambda - traceLambda1 - ...
        ob.quadFormLambda1) + dPriorLambda;
      lambda0 = log(ob.lambda); dlambda0 = dlambda .* ob.lambda;
      
      % Hessian approximation
      if fM.fS.useHessNoise
        d2lambda0 = -.5*ob.lambda.*(traceLambda1 + ob.quadFormLambda1);
        d2lambda0 = d2lambda0.*sign(-d2lambda0);
        
        if fM.iter <= 1
          ob.dlambda0Avg = dlambda0;
          ob.d2lambda0Avg = d2lambda0;
        else
          gam1 = fM.fS.gam1; gam2 = fM.gam2;
          ob.dlambda0Avg = gam1*ob.dlambda0Avg + (1-gam1)*dlambda0;
          ob.d2lambda0Avg = gam2*ob.d2lambda0Avg + (1-gam2)*d2lambda0;
        end
        
        if fM.iter >= fM.fS.nWarmupIterations
          ob.deltaLambda0 = fM.fS.momNoise*ob.deltaLambda0 - ...
            (fM.hessAlpha/ob.d2lambda0Avg)*ob.dlambda0Avg;
        else
          ob.deltaLambda0 = fM.hessAlpha*fM.fS.gradGainNoise*dlambda0;
        end
      else
        ob.deltaLambda0 = fM.hessAlpha*fM.fS.gradGainNoise*dlambda0;
      end
      
      bigDeltaInd = find(abs(ob.deltaLambda0) > fM.fS.maxDelta);
      ob.deltaLambda0(bigDeltaInd) = fM.fS.maxDelta*sign(ob.deltaLambda0(bigDeltaInd));
      
      % Take step
      lambda0 = lambda0 + ob.deltaLambda0; ob.lambda = exp(lambda0);
      ob.lambda(ob.lambda > ob.lambdaHB(2)) = ob.lambdaHB(2);
      ob.lambda(ob.lambda < ob.lambdaHB(1)) = ob.lambdaHB(1);
      ob.lambdaVec(:,fM.iter) = ob.lambda;
      
      if fM.fS.isGLMNoiseAR

        da = -.5*(ob.quadFormA1 + ob.traceA1) + dPriorA;
        a0 = log((ob.a+1)./2) - log((1-ob.a)./2);
        da0 = da .* (1-ob.a.^2)./2;
        a0 = a0 + fM.hessAlpha*fM.fS.gradGainNoise.*da0; 
        if fM.fS.checkARStability
          aPrel = 2*exp(a0)./(exp(a0)+1) - 1;
          unstableIncrease = 0;
          for i = 1:N
            if any(abs(roots([1,-aPrel(:,i)'])) >= (1 - ...
                                            fM.fS.ARCheckUnitCircleOffset))
              if ob.unstableInd(i) == 0; unstableIncrease = 1; end
              ob.unstableInd(i) = 1;
            else
              ob.unstableInd(i) = 0;              
            end
            
          end
          if unstableIncrease
            disp(['Warning: Near AR model instability for ',...
              num2str(sum(ob.unstableInd)), ' voxels.']); 
          end
          ob.a(:,~ob.unstableInd) = aPrel(:,~ob.unstableInd);
        else
          ob.a = 2*exp(a0)./(exp(a0)+1) - 1;          
        end
        ob.a(ob.a > ob.aHB(2)) = ob.aHB(2);
        ob.a(ob.a < ob.aHB(1)) = ob.aHB(1);
        ob.aVec(:,:,fM.iter) = ob.a;
      end
      
    end
    function initializeParameters(ob,fM)
      % Initialize parameters using non-spatial model
      fMTemp = FMRIModel(fM.fS); % could this be replaced with copy(fM)?
      fMTemp.X = fM.X; fMTemp.Y = fM.Y;
      fMTemp.maskInd = fM.maskInd; fMTemp.N = fM.N;
      fMTemp.fS.solveMethod = 'Chol';
      fMTemp.fS.maxiter = fM.fS.preOptNoiseMaxiter;
      fMTemp.fS.preOptNoiseParams = 0;
      fMTemp.fS.useRobbinsMonro = 0;
      fMTemp.fS.doPlot = 0;
      fMTemp.fS.doPrint = 0;
      fMTemp.priorList = cell(fMTemp.fS.K,1);
      for k = 1:fMTemp.fS.K
        fMTemp.priorList{k} = Eye(fMTemp,k,[1e-15],[1]);
      end
      fMTemp.doPrecalculations;
      fMTemp.nM = NoiseModel(fMTemp);
      fMTemp.initialize;
      if fM.fS.doPrint; disp(['Initialize noise parameters, Time: ',datestr(now,'dd-mmm-yyyy HH:MM:SS')]); end
      while(fMTemp.iter <= fMTemp.fS.maxiter)
        fMTemp.gradientStep;
      end
      ob.lambda = fMTemp.nM.lambda;
      if ob.isGLMNoiseAR; ob.a = fMTemp.nM.a; end;
    end
    function doPolyakAveraging(ob,fM)
      ob.lambda = mean(ob.lambdaVec(:,fM.iter-fM.fS.nPolyakVal+1:fM.iter),2);
      ob.a = mean(ob.aVec(:,:,fM.iter-fM.fS.nPolyakVal+1:fM.iter),3);
      
    end
    
    function ARVar = computeARVariances(ob,fM)
      N = fM.N; P = ob.P;
      ARVar = zeros(N,1);
      for n = 1:N
        an = ob.a(:,n)';
        An = tril(toeplitz([1 -an])) + [zeros(P,1), ...
              rot90(tril(toeplitz(-an(P:-1:1))),3);zeros(1,P+1)];
        r = An \ [1;zeros(P,1)];
        ARVar(n) = r(1);
      end
      ARVar = ARVar ./ ob.lambda;
    end
    function str = getParameterString(ob)
      str = ['lambda_1: ' ,num2str(ob.lambda(1))];
      if ob.isGLMNoiseAR; str = [str,', a_1: ' ,num2str(ob.a(1,1))]; end
    end
    
  end
  
end

