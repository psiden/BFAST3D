classdef Eye < handle
  %Eye Global shrinkage regression coefficient prior
  %
  % AUTHOR:       Per Siden
  %               Division of Statistics and Machine Learning
  %               Department of Computer and Information Science
  %               Linkoping University
  %
  % FIRST VER.:   2017-09-06
  
  properties
    
    priorName           % Name of the prior
    startParams         % 1x1 vector giving start values of parameters [tau2]
    fixParams           % 1x1 vector stating fixed parameters as 1s, e.g. [1]
                        % for fixed tau2
    Qk                  % Prior precision matrix
    k                   % Prior number
    pCholQk             % pseudo Cholesky factor of Q so that
                        % pCholQk*pCholQk' = Qk    
    paramPriorName   
    tau2
    tau2Vec             % Results storage
    N, ndim
    PCPalpha, PCPsigma0, PCPlambda2 
    tau2HB
    
    kappa2, kappa2Vec   % Solely for plotting purposes
        
  end
  
  methods
    
    function ob = Eye(fM,k,startParams,fixParams)
      ob.priorName = 'Eye';
      ob.k = k;
      ob.N = fM.N;
      ob.ndim = fM.fS.ndim;
      
      % Start values
      if nargin < 3
        ob.startParams = 1e-3;
      else
        ob.startParams = startParams;
      end
      if nargin < 4
        ob.fixParams = 0;
      else
        ob.fixParams = fixParams;
      end
      ob.tau2 = ob.startParams(1);
      
      ob.paramPriorName = 'PC';
      ob.PCPalpha = 0.05;
      
      % tau2 PC prior
      ob.PCPsigma0 = 100;
      ob.PCPlambda2 = -log(ob.PCPalpha) / ob.PCPsigma0;
      ob.tau2HB = [1e-5,10];
      
      % Hard bounds
      ob.tau2HB = [1e-12,Inf];      
      
      % Results storage matrices
      ob.tau2Vec = zeros(fM.fS.maxiter+1,1);ob.tau2Vec(1) = ob.tau2;
      
      % Solely for plotting purposes
      ob.kappa2 = 0;ob.kappa2Vec = zeros(fM.fS.maxiter+1,1);
      
    end
    
    function Qk = computeQk(ob)
      ob.Qk = ob.tau2*speye(ob.N);
      Qk = ob.Qk;
    end
    
    function samples = getInversePriorSamples(ob,Ns)
      ob.pCholQk = sqrt(ob.tau2) * speye(ob.N);
      samples = ob.pCholQk * randn(ob.N,Ns);
    end
    
    function gradientStep(ob,fM)
      if ~all(ob.fixParams)
        N = ob.N; k = ob.k;
        if strcmp(fM.fS.solveMethod,'Chol')
          iQTildek = fM.iQTilde(N*(k-1)+1:N*k,N*(k-1)+1:N*k); % CAN THIS BE DONE FASTER??
          traceTau1 = sum(diag(iQTildek));

        elseif strcmp(fM.fS.solveMethod,'PCG')

          traceTau1 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                               fM.Vs(N*(k-1)+1:N*k,:)));
        end

        % Prior contributions
        dPriorTau2 = -3./(2*ob.tau2) + ob.PCPlambda2/2.*ob.tau2.^(-3/2);
        if fM.fS.useHess
          d2PriorTau0 = -ob.PCPlambda2/4*ob.tau2^(-1/2);
        end

        % Reparameterize
        tau0 = log(ob.tau2); 
        
        % Compute gradient and Hessian approximation
        dtau2 = N/(2*ob.tau2) - 0.5*traceTau1 - 0.5*fM.M(k,:)*fM.M(k,:)' + dPriorTau2;
        dtau0 = dtau2 .* ob.tau2;
        if fM.fS.useHess
          d2tau0 = -0.5*fM.M(k,:)*ob.Qk*fM.M(k,:)' - 0.5*ob.tau2*traceTau1 ...
            + d2PriorTau0;
          d2tau0 = d2tau0.*sign(-d2tau0);
        end

        % Take step
        tau0 = tau0 + fM.fS.gradGain.*dtau0; ob.tau2 = exp(tau0);

        % Hard bounds
        ob.tau2(ob.tau2 > ob.tau2HB(:,2)) = ob.tau2HB(ob.tau2 > ob.tau2HB(:,2),2);
        ob.tau2(ob.tau2 < ob.tau2HB(:,1)) = ob.tau2HB(ob.tau2 < ob.tau2HB(:,1),1);
      end
      % Save in vector
      ob.tau2Vec(fM.iter) = ob.tau2;
    end    
    function doPolyakAveraging(ob,fM)
      ob.tau2 = mean(ob.tau2Vec(fM.iter-fM.fS.nPolyakVal+1:fM.iter));     
    end
    function params = getParameters(ob)
      params = ob.tau2;
    end
    function str = getParameterString(ob)
      str = ['tau2_',num2str(ob.k),': ',num2str(ob.tau2)];
    end
  end
  
end

