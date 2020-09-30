classdef Matern1Iso < handle
  %Matern2Iso Matern regression coefficient prior with
  %   alpha = 1, isotropic
  %
  % AUTHOR:       Per Siden
  %               Division of Statistics and Machine Learning
  %               Department of Computer and Information Science
  %               Linkoping University
  %
  % FIRST VER.:   2017-10-07
  % REVISED.:     2020-09-30
  
  properties
    
    priorName           % Name of the prior
    startParams         % 1x2 vector giving start values of parameters [tau2,kappa2]
    fixParams           % 1x2 vector stating fixed parameters as 1s, e.g. [0,1]
                        % for non-fixed tau2, fixed kappa2
    Qk                  % Prior precision matrix
    k                   % Prior number
    C                   % SPDE C matrix
    Ci                  % Diagonal (mass lumped) version of C
    cholCi
    G                   % SPDE G
    pCholG              % pseudo Chol factor of G s.t. pCholG'*pCholG=G
    paramPriorName      % Prior for tau2 and kappa2
    isIntrinsic         % If kappa == 0
    tau2, kappa2
    tau2Vec, kappa2Vec  % Results storage
    N, ndim
    PCPalpha            % PC prior alpha
    PCPsigma0, PCPlambda1, PCPlambda2, PCPlambda3, Cnu  
    logNPmu, logNPsigma % logN prior parameters
    q1, q2              % Gamma prior parameters
    kappaPriorMean      % Hard bounds
    kappa2HB, tau2HB
    Kk, reoK, ireoK
    lastiKCVsPool
    dtau0Avg, dkappa0Avg
    d2tau0Avg, d2kappa0Avg
    deltaTau0, deltaKappa0
    
  end
  
  methods
    
    function ob = Matern1Iso(fM,k,startParams,fixParams)
      ob.priorName = 'Matern1Iso';
      ob.k = k;
      ob.N = fM.N;
      ob.ndim = fM.fS.ndim;
      [QList,GList] = setupPrecMats({'LI'},ob.N,fM.fS.sz,fM.maskInd,ob.ndim);
      ob.G = QList{1};
      ob.C = speye(ob.N);
      ob.pCholG = GList{1}';
      ob.Ci = spdiags(1./diag(ob.C),0,ob.N,ob.N);
      ob.cholCi = chol(ob.Ci);
      ob.lastiKCVsPool = zeros(ob.N,fM.fS.SHutchRandomPool);
                  
      % Start values
      if nargin < 3
        ob.startParams = [1e-2,.1]; % ob.kappaPriorMean^2;
      else
        ob.startParams = startParams;
      end
      if nargin < 4
        ob.fixParams = [0,0];
      else
        ob.fixParams = fixParams;
      end      
      
      ob.isIntrinsic = (ob.startParams(2) == 0 && ob.fixParams(2) == 1);
      ob.PCPalpha = 0.05;
      
      if ob.isIntrinsic
        ob.paramPriorName = 'PC'; % 'PC_cond'; % 'Gamma'; % 
        
        if strcmp(ob.paramPriorName,'PC')
          % tau2 PC prior
          ob.PCPsigma0 = 2 / fM.fS.sf;
          if ob.ndim == 2
            error('No support for this ICAR hyperparameter prior in 2D.')
          elseif ob.ndim == 3
            ob.PCPlambda2 = -log(ob.PCPalpha) * sqrt(0.29) / ob.PCPsigma0;          
          end
        elseif strcmp(ob.paramPriorName,'PC_cond')
          % tau2 PC prior
          ob.PCPsigma0 = 0.5 / fM.fS.sf;
          if ob.ndim == 2
            ob.PCPlambda2 = -log(ob.PCPalpha) / (sqrt(4)*ob.PCPsigma0);
          elseif ob.ndim == 3
            ob.PCPlambda2 = -log(ob.PCPalpha) / (sqrt(6)*ob.PCPsigma0);          
          end
        elseif strcmp(ob.paramPriorName,'Gamma')
          ob.q1 = 10; ob.q2 = 0.1;
        end
        ob.tau2HB = [1e-5,10];        
      else
        ob.paramPriorName = 'logN';
        ob.logNPmu = log([1e-2,1e-1]);
        ob.logNPsigma = [4,1];
        
        % Hard bounds
        ob.tau2HB = [1e-5 10];
        ob.kappa2HB = [1e-5 10];
      end
            
      ob.tau2 = ob.startParams(1);
      ob.kappa2 = ob.startParams(2); 
      
      % Results storage matrices
      ob.tau2Vec = zeros(fM.fS.maxiter+1,1);ob.tau2Vec(1) = ob.tau2;
      ob.kappa2Vec = zeros(fM.fS.maxiter+1,1);ob.kappa2Vec(1) = ob.kappa2;
      
    end
    function Qk = computeQk(ob)
      ob.Kk = ob.G+ob.kappa2*ob.C;
      ob.reoK = amd(ob.Kk);ob.ireoK(ob.reoK) = 1:ob.N;
      ob.Qk = ob.tau2*ob.Kk;
      Qk = ob.Qk;
    end
    function samples = getInversePriorSamples(ob,Ns)
      samples = sqrt(ob.tau2) * (ob.pCholG * randn(size(ob.pCholG,2),Ns) + ...
        sqrt(ob.kappa2) * diag(sqrt(diag(ob.C))) * randn(ob.N,Ns));
    end
    function gradientStep(ob,fM)
      if ~all(ob.fixParams)
        N = ob.N; k = ob.k;
        if strcmp(fM.fS.solveMethod,'Chol')
          iQTildek = fM.iQTilde(N*(k-1)+1:N*k,N*(k-1)+1:N*k); % CAN THIS BE DONE FASTER??
          iQtimesK = iQTildek.*ob.Kk;    
          
          if ~ob.fixParams(1)    
            traceTau1 = sum(sum(iQtimesK));
          end
          
          if ~ob.fixParams(2)
            iQtimesC = iQTildek.*ob.C;
            RKReo = chol(ob.Kk(ob.reoK,ob.reoK));
            iKReo = Qinv(RKReo);
            traceKappa1 = sum(sum(iKReo.*ob.C(ob.reoK,ob.reoK)));
            traceKappa2 = sum(sum(iQtimesC));
            if fM.fS.useHess
              KCiKk = ob.Kk*ob.Ci*ob.Kk;
              reoKCiKk = amd(KCiKk);ireoKCiKk(reoKCiKk) = 1:N;
              cholKCiKk = chol(KCiKk(reoKCiKk,reoKCiKk));
              iKCiKkReo = Qinv(cholKCiKk); 
              iKCiKk = iKCiKkReo(ireoKCiKk,ireoKCiKk);
              traceKappa3 = sum(sum(iKCiKk.*ob.C));
            end
          end          

        elseif strcmp(fM.fS.solveMethod,'PCG')
          
          if ~ob.fixParams(1)
            traceTau1 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                 ob.Kk*fM.Vs(N*(k-1)+1:N*k,:)));
          end
                    
          if ~ob.fixParams(2)
            traceKappa2 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                   ob.C*fM.Vs(N*(k-1)+1:N*k,:)));  
            cholKkCi = ob.cholCi' * ob.Kk * ob.cholCi;
            KkCiReo = cholKkCi(ob.reoK,ob.reoK);
            if fM.fS.useHutchRandomPool
              lastiKCVs = ob.lastiKCVsPool(:,fM.randPool.lastInd);
            else
              lastiKCVs = zeros(N,fM.fS.Ns);
            end
            VsK = fM.Vs((k-1)*N+1:k*N,:);
            VsReo = VsK(ob.reoK,:);
            iKCVs = zeros(N,fM.fS.Ns);
            try [M1,~] = icholSafe(KkCiReo,fM.startIcholEps2,fM.fS.maxIcholEps);
            catch; M1 = spdiags(1./sqrt(diag(KkCiReo)),0,N,N); end;
            for i = 1:fM.fS.Ns
              [sampK,flag,relres] = pcg(KkCiReo,VsReo(:,i),fM.fS.PCGtol,...
                          fM.fS.PCGmaxiter,M1,M1',lastiKCVs(ob.reoK,i));
              iKCVs(:,i) = sampK(ob.ireoK);
            end
            if fM.fS.useHutchRandomPool; ob.lastiKCVsPool(:,fM.randPool.lastInd) = iKCVs;end
            traceKappa1 = mean(dot(iKCVs,VsK));

            if fM.fS.useHess
              traceKappa3 = mean(dot(iKCVs,iKCVs));
            end
          end
        end

        % Prior contributions
        if ob.isIntrinsic
          if strcmp(ob.paramPriorName,'PC') || strcmp(ob.paramPriorName,'PC_cond')
            dPriorTau2 = -3./(2*ob.tau2) + ob.PCPlambda2/2.*ob.tau2.^(-3/2);
            if fM.fS.useHess
              d2PriorTau0 = -ob.PCPlambda2/4*ob.tau2^(-1/2);
            end
          elseif strcmp(ob.paramPriorName,'Gamma')
            dPriorTau2 = (ob.q2-1)./ob.tau2 - 1/ob.q1;            
            if fM.fS.useHess
              d2PriorTau0 = - 1/ob.q1*ob.tau2;
            end
          end
        else
          if strcmp(ob.paramPriorName,'logN')
            dPriorTau2 = - (log(ob.tau2) - ob.logNPmu(1)) / ...
                           (ob.logNPsigma(1)^2*ob.tau2);
            dPriorKappa2 = - (log(ob.kappa2) - ob.logNPmu(2)) / ...
                             (ob.logNPsigma(2)^2*ob.kappa2);
            if fM.fS.useHess
              d2PriorTau0 = - 1/ob.logNPsigma(1)^2;
              d2PriorKappa0 = - 1/ob.logNPsigma(2)^2;
            end
          else
            dPriorTau2 = 0;
            dPriorKappa2 = 0;
            if fM.fS.useHess
              d2PriorTau0 = 0;
              d2PriorKappa0 = 0;
            end
          end
        end

        % Reparameterize
        tau0 = log(ob.tau2); kappa0 = log(ob.kappa2); 
        
        % Compute gradient and Hessian approximation
        if ~ob.fixParams(1)
          dtau2 = N/(2*ob.tau2) - 0.5*traceTau1 - ...
                    0.5*fM.M(k,:)*ob.Kk*fM.M(k,:)' + dPriorTau2;
          dtau0 = dtau2 .* ob.tau2;
          if fM.fS.useHess
            d2tau0 = -0.5*fM.M(k,:)*ob.Qk*fM.M(k,:)' - 0.5*ob.tau2*traceTau1 ...
              + d2PriorTau0;
            d2tau0 = d2tau0.*sign(-d2tau0);
          end
        else
          dtau0 = 0; d2tau0 = -1;
        end
        
        if ~ob.fixParams(2)
          dkappa2 = 0.5*traceKappa1 - 0.5*ob.tau2*traceKappa2 - ...
                        0.5*ob.tau2*fM.M(k,:)*ob.C*fM.M(k,:)' + dPriorKappa2;
          dkappa0 = dkappa2 .* ob.kappa2;
          if fM.fS.useHess
            d2kappa0 = 0.5*ob.kappa2*(traceKappa1 - ob.kappa2*traceKappa3 ...
                -ob.tau2*fM.M(k,:)*ob.C*fM.M(k,:)' - ob.tau2*traceKappa2) + d2PriorKappa0;  
            d2kappa0 = d2kappa0.*sign(-d2kappa0);      
          end        
        else
          dkappa0 = 0; d2kappa0 = -1;
        end

        if fM.fS.useHess
          if fM.iter <= 1
            ob.dtau0Avg = dtau0; ob.dkappa0Avg = dkappa0;
            ob.d2tau0Avg = d2tau0; ob.d2kappa0Avg = d2kappa0;
          else
            gam1 = fM.fS.gam1; gam2 = fM.gam2;
            ob.dtau0Avg = gam1*ob.dtau0Avg + (1-gam1)*dtau0;
            ob.dkappa0Avg = gam1*ob.dkappa0Avg + (1-gam1)*dkappa0;
            ob.d2tau0Avg = gam2*ob.d2tau0Avg + (1-gam2)*d2tau0;
            ob.d2kappa0Avg = gam2*ob.d2kappa0Avg + (1-gam2)*d2kappa0;     
          end

          if fM.iter >= fM.fS.nWarmupIterations
            ob.deltaTau0 = fM.fS.momPrior*ob.deltaTau0 - ...
                           (fM.hessAlpha/ob.d2tau0Avg)*ob.dtau0Avg;
            ob.deltaKappa0 = fM.fS.momPrior*ob.deltaKappa0 - ...
                             (fM.hessAlpha/ob.d2kappa0Avg)*ob.dkappa0Avg;
          else      
            ob.deltaTau0 = fM.fS.gradGain*dtau0;
            ob.deltaKappa0 = fM.fS.gradGain*dkappa0;
          end
        else
          ob.deltaTau0 = fM.fS.gradGain*dtau0;
          ob.deltaKappa0 = fM.fS.gradGain*dkappa0;        
        end
        if abs(ob.deltaTau0) > fM.fS.maxDelta;ob.deltaTau0 = ...
                                      fM.fS.maxDelta*sign(ob.deltaTau0);end
        if abs(ob.deltaKappa0) > fM.fS.maxDelta;ob.deltaKappa0 = ...
                                      fM.fS.maxDelta*sign(ob.deltaKappa0);end

        % Take step and check hard bounds
        if ~ob.fixParams(1)
          tau0 = tau0 + ob.deltaTau0; ob.tau2 = exp(tau0);
          ob.tau2(ob.tau2 > ob.tau2HB(:,2)) = ob.tau2HB(ob.tau2 > ob.tau2HB(:,2),2);
          ob.tau2(ob.tau2 < ob.tau2HB(:,1)) = ob.tau2HB(ob.tau2 < ob.tau2HB(:,1),1);
        end
        if ~ob.fixParams(2)
          kappa0 = kappa0 + ob.deltaKappa0; ob.kappa2 = exp(kappa0);
          ob.kappa2(ob.kappa2 > ob.kappa2HB(:,2)) = ob.kappa2HB(ob.kappa2 > ob.kappa2HB(:,2),2);
          ob.kappa2(ob.kappa2 < ob.kappa2HB(:,1)) = ob.kappa2HB(ob.kappa2 < ob.kappa2HB(:,1),1);
        end
      end
      
      % Save in vector
      ob.tau2Vec(fM.iter) = ob.tau2; ob.kappa2Vec(fM.iter) = ob.kappa2; 
      
    end
    function doPolyakAveraging(ob,fM)
      ob.tau2 = mean(ob.tau2Vec(fM.iter-fM.fS.nPolyakVal+1:fM.iter));  
      ob.kappa2 = mean(ob.kappa2Vec(fM.iter-fM.fS.nPolyakVal+1:fM.iter));    
    end
    function params = getParameters(ob)
      params = [ob.tau2,ob.kappa2];
    end
    function str = getParameterString(ob)
      str = ['tau2_',num2str(ob.k),': ',num2str(ob.tau2), ', kappa2_',...
             num2str(ob.k),': ' ,num2str(ob.kappa2)];
    end
  end
end

