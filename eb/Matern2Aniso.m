classdef Matern2Aniso < handle
  %Matern2Iso Matern regression coefficient prior with
  %   alpha = 2, anisotropic
  %   Code not safe for 2D nor for intrinsic prior
  %
  % AUTHOR:       Per Siden
  %               Division of Statistics and Machine Learning
  %               Department of Computer and Information Science
  %               Linkoping University
  %
  % FIRST VER.:   2017-11-14
  % REVISED.:     2020-09-30
  
  properties
    
    priorName           % Name of the prior
    startParams         % 1x4 vector giving start values of parameters [tau2,kappa2,hx,hy]
    fixParams           % 1x2 vector stating fixed parameters as 1s, e.g. [0,1,1,1]
                        % for non-fixed tau2, fixed kappa2
    Qk                  % Prior precision matrix
    k                   % Prior number
    pCholQk             % pseudo Cholesky factor of Q so that
                        % pCholQk'*pCholQk = Qk    
    C                   % SPDE C matrix
    Ci                  % Diagonal (mass lumped) version of C
    cholCi
    Gx, Gy, Gz          % SPDE G
    Lx, Ly, Lz          % Pseudo-Cholesky factors such that e.g. Gx = Lx'*Lx
    nu                  % alpha-d/2
    paramPriorName      % Prior for tau2 and kappa2
    isIntrinsic         % If kappa == 0
    tau2, kappa2
    hx, hy
    tau2Vec, kappa2Vec  
    hxVec, hyVec        % Results storage
    gradGainAniso       % Addtional gradient gain factor for h parameters
    N, ndim
    PCPalpha, PCPrho0   % PC prior alpha
    PCPsigma0, PCPlambda1, PCPlambda2, PCPlambda3, Cnu    
    rhoHardBound        % Hard bounds
    sigmaHardBound
    kappaPriorMean
    tau2HB, kappa2HB
    sigmah, hxHB, hyHB % h prior parameters
    Kk, KCiKk, reoK, ireoK
    lastiKCVsPool, lastiKdKdh0xVsPool, lastiKdKdh0yVsPool
    dtau0Avg, dkappa0Avg
    d2tau0Avg, d2kappa0Avg
    deltaTau0, deltaKappa0
    dh0xAvg, dh0yAvg, d2h0xAvg, d2h0yAvg, deltah0x, deltah0y
    dtau0Vec, dkappa0Vec, dh0xVec, dh0yVec
    dtau0AvgVec, dkappa0AvgVec, dh0xAvgVec, dh0yAvgVec % Gradient storage
    d2tau0Vec, d2kappa0Vec, d2h0xVec, d2h0yVec
    d2tau0AvgVec, d2kappa0AvgVec, d2h0xAvgVec, d2h0yAvgVec % Hessian storage
    
  end
  
  methods
    
    function ob = Matern2Aniso(fM,k,startParams,fixParams)
      ob.priorName = 'Matern2Aniso';
      ob.gradGainAniso = .1; % .04;
      ob.k = k;
      ob.N = fM.N;
      ob.ndim = fM.fS.ndim;
      [QxList,QyList,QzList,GxList,GyList,GzList] = setupAnisoPrecMats(fM.fS.K,ob.N,fM.fS.sz,fM.maskInd,ob.ndim);
      ob.Gx = QxList{1}; ob.Gy = QyList{1}; ob.Gz = QzList{1};
      ob.Lx = GxList{1}; ob.Ly = GyList{1}; ob.Lz = GzList{1};
      ob.C = speye(ob.N);
      ob.Ci = spdiags(1./diag(ob.C),0,ob.N,ob.N);
      ob.cholCi = chol(ob.Ci);
      ob.nu = 2 - ob.ndim/2;
      ob.lastiKCVsPool = zeros(ob.N,fM.fS.SHutchRandomPool);
      ob.lastiKdKdh0xVsPool = zeros(ob.N,fM.fS.SHutchRandomPool);
      ob.lastiKdKdh0yVsPool = zeros(ob.N,fM.fS.SHutchRandomPool);
      
      % Start values
      if nargin < 3
        ob.startParams = [1e-2,.1,1,1];
      else
        ob.startParams = startParams;
      end
      if nargin < 4
        ob.fixParams = [0,0,0,0];
      else
        ob.fixParams = fixParams;
      end
      
      ob.isIntrinsic = (ob.startParams(2) == 0 && ob.fixParams(2) == 1);
      ob.paramPriorName = 'PC';
      ob.PCPalpha = 0.05;
      
      if ob.isIntrinsic
        % This case should be avoided as inv(K) is not defined for kappa2=0
%         % tau2 PC prior
%         ob.PCPsigma0 = 4;
%         if ob.ndim == 2
%           ob.PCPlambda2 = -log(ob.PCPalpha) / (sqrt(20)*ob.PCPsigma0);
%         elseif ob.ndim == 3
%           ob.PCPlambda2 = -log(ob.PCPalpha) / (sqrt(42)*ob.PCPsigma0);          
%         end
%         ob.tau2HB = [1e-5,10];
      else
        % kappa2 and tau2 PC priors Fuglstad2016 notation
        ob.PCPrho0 = 2;
        ob.PCPsigma0 = 2 / fM.fS.sf;
        ob.PCPlambda1 = -log(ob.PCPalpha)*...
          (ob.PCPrho0/sqrt(8*ob.nu))^(ob.ndim/2);
        ob.Cnu = gamma(ob.nu) / (gamma(ob.nu+ob.ndim/2)*...
          ((4*pi)^(ob.ndim/2)));
        ob.PCPlambda3 = -sqrt(ob.Cnu) * log(ob.PCPalpha) / ob.PCPsigma0;
        
        % Hard bounds
        ob.rhoHardBound = [1 40];
        ob.sigmaHardBound = [.1 50];
        ob.kappaPriorMean = ob.PCPlambda1^(-2/ob.ndim);
        ob.kappa2HB = 8*ob.nu ./ (ob.rhoHardBound([2,1]).^2);
        ob.tau2HB = ob.Cnu ./ (ob.kappaPriorMean.^(2*ob.nu) .* ...
          ob.sigmaHardBound(:,[2,1]).^2);
      end
      % hx, hy logN prior 
      if fM.fS.ndim == 3
        ob.sigmah = 0.1;
      end
      ob.hxHB = [.05,20];
      if fM.fS.ndim == 3; ob.hyHB = [.05,20]; end
        
      ob.tau2 = ob.startParams(1);
      ob.kappa2 = ob.startParams(2); 
      ob.hx = ob.startParams(3); 
      if fM.fS.ndim == 3; ob.hy = ob.startParams(4); end
      
      % Results storage matrices
      ob.tau2Vec = zeros(fM.fS.maxiter+1,1);ob.tau2Vec(1) = ob.tau2;
      ob.kappa2Vec = zeros(fM.fS.maxiter+1,1);ob.kappa2Vec(1) = ob.kappa2;
      ob.hxVec = zeros(fM.fS.maxiter+1,1);ob.hxVec(1) = ob.hx;
      if fM.fS.ndim == 3
        ob.hyVec = zeros(fM.fS.maxiter+1,1);ob.hyVec(1) = ob.hy;
      end
      
      % Gradient storage matrices
      ob.dtau0AvgVec = zeros(fM.fS.maxiter+1,1);
      ob.dkappa0AvgVec = zeros(fM.fS.maxiter+1,1);
      ob.dh0xAvgVec = zeros(fM.fS.maxiter+1,1);
      ob.dh0yAvgVec = zeros(fM.fS.maxiter+1,1);
      
      % Hessian storage matrices
      ob.d2tau0AvgVec = zeros(fM.fS.maxiter+1,1);
      ob.d2kappa0AvgVec = zeros(fM.fS.maxiter+1,1);
      ob.d2h0xAvgVec = zeros(fM.fS.maxiter+1,1);
      ob.d2h0yAvgVec = zeros(fM.fS.maxiter+1,1);      
      
    end
    function Qk = computeQk(ob)
      ob.Kk = ob.hx*ob.Gx+ob.hy*ob.Gy+(1/(ob.hx*ob.hy))*ob.Gz+ob.kappa2*ob.C;
      ob.reoK = amd(ob.Kk);ob.ireoK(ob.reoK) = 1:ob.N;
      ob.KCiKk = ob.Kk*ob.Ci*ob.Kk;
      ob.Qk = ob.tau2*ob.KCiKk;
      Qk = ob.Qk;
    end
    function samples = getInversePriorSamples(ob,Ns)
      ob.Kk = ob.hx*ob.Gx+ob.hy*ob.Gy+(1/(ob.hx*ob.hy))*ob.Gz+ob.kappa2*ob.C;
      ob.pCholQk = sqrt(ob.tau2) * ob.Kk * diag(sqrt(diag(ob.Ci)));
      samples = ob.pCholQk * randn(size(ob.pCholQk,1),Ns);
    end
    function gradientStep(ob,fM)
      if ~all(ob.fixParams)
        N = ob.N; k = ob.k;
        
        if ~ob.fixParams(3)
          dKdh0x = ob.hx*ob.Gx - (1/(ob.hx*ob.hy))*ob.Gz;
          H1x = dKdh0x*ob.Ci*ob.Kk;
          if fM.fS.useHess
            dKdh0xSq = dKdh0x*dKdh0x;
            d2Kdh0x2 = ob.hx*ob.Gx + (1/(ob.hx*ob.hy))*ob.Gz;
            H2x = dKdh0xSq + ob.Kk*d2Kdh0x2; 
          end
        end
        
        if ~ob.fixParams(4)
          dKdh0y = ob.hy*ob.Gy - (1/(ob.hx*ob.hy))*ob.Gz;
          H1y = dKdh0y*ob.Ci*ob.Kk;
          if fM.fS.useHess
            dKdh0ySq = dKdh0y*dKdh0y;
            d2Kdh0y2 = ob.hy*ob.Gy + (1/(ob.hx*ob.hy))*ob.Gz;
            H2y = dKdh0ySq + ob.Kk*d2Kdh0y2; 
          end
        end
        
        if strcmp(fM.fS.solveMethod,'Chol')
          iQTildek = fM.iQTilde(N*(k-1)+1:N*k,N*(k-1)+1:N*k); % CAN THIS BE DONE FASTER??
          
          if ~ob.fixParams(1)
            iQtimesKCiK = iQTildek.*ob.KCiKk;          
            traceTau1 = sum(sum(iQtimesKCiK));
          end
          
          if any(~ob.fixParams(2:4))
            iQtimesK = iQTildek.*ob.Kk;
            RKReo = chol(ob.Kk(ob.reoK,ob.reoK));
            iKReo = Qinv(RKReo);
          end
            
          if ~ob.fixParams(2)
            traceKappa1 = sum(sum(iKReo.*ob.C(ob.reoK,ob.reoK)));
            traceKappa2 = sum(sum(iQtimesK));
            if fM.fS.useHess
              iQtimesC = iQTildek.*ob.C;
              reoQprior = amd(ob.Qk);ireoQprior(reoQprior) = 1:N;
              cholQk = chol(ob.Qk(reoQprior,reoQprior));
              iQkReo = Qinv(cholQk);
              iQk = iQkReo(ireoQprior,ireoQprior);
              traceKappa3 = sum(sum(iQk.*ob.C));
              traceKappa4 = sum(sum(iQtimesC));
            end
          end     
          
          if ~ob.fixParams(3)
            iQtimesH1x = iQTildek.*H1x;      
            traceh0x1 = sum(sum(iKReo.*dKdh0x(ob.reoK,ob.reoK)));
            traceh0x2 = sum(sum(iQtimesH1x));
            if fM.fS.useHess
              disp('Warning: solveMethod Chol is not compatible with useHess.');
              traceh0x3 = 0;         
              traceh0x4 = sum(sum(iKReo.*d2Kdh0x2(ob.reoK,ob.reoK)));    
              traceh0x5 = sum(sum(iQTildek.*H2x));
            end
          end

          if ~ob.fixParams(4)
            iQtimesH1y = iQTildek.*H1y;      
            traceh0y1 = sum(sum(iKReo.*dKdh0y(ob.reoK,ob.reoK)));
            traceh0y2 = sum(sum(iQtimesH1y));
            if fM.fS.useHess
              disp('Warning: solveMethod Chol is not compatible with useHess.');
              traceh0y3 = 0;         
              traceh0y4 = sum(sum(iKReo.*d2Kdh0y2(ob.reoK,ob.reoK)));    
              traceh0y5 = sum(sum(iQTildek.*H2y));
            end
          end

        elseif strcmp(fM.fS.solveMethod,'PCG')
          
          if ~ob.fixParams(1)
            traceTau1 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                 ob.KCiKk*fM.Vs(N*(k-1)+1:N*k,:)));
          end
                    
          if any(~ob.fixParams(2:4))
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
              if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for KCiK trace.']); end
              iKCVs(:,i) = sampK(ob.ireoK);
            end
            if fM.fS.useHutchRandomPool; ob.lastiKCVsPool(:,fM.randPool.lastInd) = iKCVs;end
          end
          
          if ~ob.fixParams(2)
            traceKappa1 = mean(dot(iKCVs,VsK));
            traceKappa2 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                   ob.Kk*fM.Vs(N*(k-1)+1:N*k,:)));  
            if fM.fS.useHess
              traceKappa3 = mean(dot(iKCVs,iKCVs));
              traceKappa4 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                   ob.C*fM.Vs(N*(k-1)+1:N*k,:)));
            end
          end
          
          if ~ob.fixParams(3)
            assert(all(diag(ob.C)==1),'Assuming C is an identity matrix');
            traceh0x1 = mean(dot(iKCVs,dKdh0x*VsK));
            traceh0x2 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                   H1x*fM.Vs(N*(k-1)+1:N*k,:)));  
            if fM.fS.useHess
              if fM.fS.useHutchRandomPool
                lastiKdKdh0xVs = ob.lastiKdKdh0xVsPool(:,fM.randPool.lastInd);
              else
                lastiKdKdh0xVs = zeros(N,fM.fS.Ns);
              end
              dKdh0xVs = dKdh0x * VsK;
              dKdh0xVsReo = dKdh0xVs(ob.reoK,:);
              iKdKdh0xVs = zeros(N,fM.fS.Ns);
              for i = 1:fM.fS.Ns
                [sampK,flag,relres] = pcg(KkCiReo,dKdh0xVsReo(:,i),fM.fS.PCGtol,...
                            fM.fS.PCGmaxiter,M1,M1',lastiKdKdh0xVs(ob.reoK,i));
                if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for iKdKdh0x trace.']); end
                iKdKdh0xVs(:,i) = sampK(ob.ireoK);
              end
              if fM.fS.useHutchRandomPool; ob.lastiKdKdh0xVs(:,fM.randPool.lastInd) = iKdKdh0xVs;end            
              traceh0x3 = mean(dot(iKCVs,dKdh0x*iKdKdh0xVs));         
              traceh0x4 = mean(dot(iKCVs,dKdh0xSq*VsK));       
              traceh0x5 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                     H2x*fM.Vs(N*(k-1)+1:N*k,:))); 
            end
          end
          
          if ~ob.fixParams(4)
            assert(all(diag(ob.C)==1),'Assuming C is an identity matrix');
            traceh0y1 = mean(dot(iKCVs,dKdh0y*VsK));
            traceh0y2 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                   H1y*fM.Vs(N*(k-1)+1:N*k,:))); 
            if fM.fS.useHess
              if fM.fS.useHutchRandomPool
                lastiKdKdh0yVs = ob.lastiKdKdh0yVsPool(:,fM.randPool.lastInd);
              else
                lastiKdKdh0yVs = zeros(N,fM.fS.Ns);
              end
              dKdh0yVs = dKdh0y * VsK;
              dKdh0yVsReo = dKdh0yVs(ob.reoK,:);
              iKdKdh0yVs = zeros(N,fM.fS.Ns);
              for i = 1:fM.fS.Ns
                [sampK,flag,relres] = pcg(KkCiReo,dKdh0yVsReo(:,i),fM.fS.PCGtol,...
                            fM.fS.PCGmaxiter,M1,M1',lastiKdKdh0yVs(ob.reoK,i));
                if flag ~= 0; disp(['Warning: PCG flag=',num2str(flag),' for iKdKdh0y trace.']); end
                iKdKdh0yVs(:,i) = sampK(ob.ireoK);
              end
              if fM.fS.useHutchRandomPool; ob.lastiKdKdh0yVs(:,fM.randPool.lastInd) = iKdKdh0yVs;end            
              traceh0y3 = mean(dot(iKCVs,dKdh0y*iKdKdh0yVs));         
              traceh0y4 = mean(dot(iKCVs,dKdh0ySq*VsK));       
              traceh0y5 = mean(dot(fM.iQTildeVs(N*(k-1)+1:N*k,:),...
                                     H2y*fM.Vs(N*(k-1)+1:N*k,:))); 
            end             
          end
        end

        % Prior contributions
        if ob.isIntrinsic
%           dPriorTau2 = -3./(2*ob.tau2) + ob.PCPlambda2/2.*ob.tau2.^(-3/2);
%           if fM.fS.useHess
%             d2PriorTau0 = -ob.PCPlambda2/4*ob.tau2^(-1/2);
%           end
        else
          ndim = ob.ndim; nu = ob.nu;
          dPriorTau2 = -3./(2*ob.tau2) + ...
                      ob.PCPlambda3.*sqrt(ob.kappa2).^(-nu)./2.*ob.tau2.^(-3/2);
          dPriorKappa2 = 1./(2*sqrt(ob.kappa2)) .* ((ndim/2-1-nu)./sqrt(ob.kappa2) - ...
                ob.PCPlambda1.*ndim/2.*sqrt(ob.kappa2).^(ndim/2-1) + ...
                ob.PCPlambda3.*nu.*sqrt(ob.kappa2).^(-nu-1).*ob.tau2.^(-1/2));
          if fM.fS.useHess
            d2PriorTau0 = -ob.PCPlambda3/4*ob.kappa2^(-nu/2)*ob.tau2^(-1/2);
            d2PriorKappa0 = 1/4 * (-(ndim/2-1-nu)*ob.kappa2^(-1/2) - ...
                  (ndim/2-1)*ob.PCPlambda1*ndim/2*ob.kappa2^(1/2*(ndim/2-1)) + ...
                  (nu-1)*ob.PCPlambda3*nu*ob.kappa2^(1/2*(nu-1))*ob.tau2^(-1/2));            
          end
        end
        dPriorh0x = -2/(3*ob.sigmah^2)*(log(ob.hx)+log(ob.hy)/2);
        dPriorh0y = -2/(3*ob.sigmah^2)*(log(ob.hy)+log(ob.hx)/2);
        if fM.fS.useHess
          d2Priorh0x = -2/(3*ob.sigmah^2);
          d2Priorh0y = -2/(3*ob.sigmah^2);
        end
        
        % Reparameterize
        tau0 = log(ob.tau2); kappa0 = log(ob.kappa2); 
        h0x = log(ob.hx); h0y = log(ob.hy); 
        
        % Compute gradient and Hessian approximation
        if ~ob.fixParams(1)
          dtau2 = N/(2*ob.tau2) - 0.5*traceTau1 - ...
                    0.5*fM.M(k,:)*ob.KCiKk*fM.M(k,:)' + dPriorTau2;
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
          dkappa2 = traceKappa1 - ob.tau2*traceKappa2 - ...
                        ob.tau2*fM.M(k,:)*ob.Kk*fM.M(k,:)' + dPriorKappa2;
          dkappa0 = dkappa2 .* ob.kappa2;
          if fM.fS.useHess
            d2kappa0 = ob.kappa2*(traceKappa1 - ob.kappa2*ob.tau2*traceKappa3 ...
                -ob.tau2*fM.M(k,:)*ob.Kk*fM.M(k,:)' - ob.tau2*traceKappa2 ...
                -ob.kappa2*ob.tau2*fM.M(k,:)*ob.C*fM.M(k,:)' - ...
                ob.kappa2*ob.tau2*traceKappa4) + d2PriorKappa0;  
    %         dtau0dkappa0 = -ob.kappa2*ob.tau2*(traceKappa2 + fM.M(k,:)*ob.Kk*fM.M(k,:)');
    %         Hess = [d2tau0,dtau0dkappa0;dtau0dkappa0,d2kappa0];
            d2kappa0 = d2kappa0.*sign(-d2kappa0);      
          end        
        else
          dkappa0 = 0; d2kappa0 = -1;
        end
        
        if ~ob.fixParams(3)
          dh0x = traceh0x1 - ob.tau2*traceh0x2 - ...
                    ob.tau2*fM.M(k,:)*H1x*fM.M(k,:)' + dPriorh0x;
          if fM.fS.useHess
            d2h0x = -traceh0x3 + traceh0x4 - ob.tau2*fM.M(k,:)*H2x*fM.M(k,:)' ...
                        -ob.tau2*traceh0x5 + d2Priorh0x;
            d2h0x = d2h0x.*sign(-d2h0x);
          end
        else
          dh0x = 0; d2h0x = -1;
        end
        
        if ~ob.fixParams(4)
          dh0y = traceh0y1 - ob.tau2*traceh0y2 - ...
                    ob.tau2*fM.M(k,:)*H1y*fM.M(k,:)' + dPriorh0y;
          if fM.fS.useHess
            d2h0y = -traceh0y3 + traceh0y4 - ob.tau2*fM.M(k,:)*H2y*fM.M(k,:)' ...
                        -ob.tau2*traceh0y5 + d2Priorh0y;
            d2h0y = d2h0y.*sign(-d2h0y);
          end
        else
          dh0y = 0; d2h0y = -1;
        end

        if fM.fS.useHess
          if fM.iter <= 1
            ob.dtau0Avg = dtau0; ob.dkappa0Avg = dkappa0;
            ob.d2tau0Avg = d2tau0; ob.d2kappa0Avg = d2kappa0;
            ob.dh0xAvg = dh0x; ob.dh0yAvg = dh0y;
            ob.d2h0xAvg = d2h0x; ob.d2h0yAvg = d2h0y;
          else
            gam1 = fM.fS.gam1; gam2 = fM.gam2;
            ob.dtau0Avg = gam1*ob.dtau0Avg + (1-gam1)*dtau0;
            ob.dkappa0Avg = gam1*ob.dkappa0Avg + (1-gam1)*dkappa0;
            ob.dh0xAvg = gam1*ob.dh0xAvg + (1-gam1)*dh0x;
            ob.dh0yAvg = gam1*ob.dh0yAvg + (1-gam1)*dh0y;
            ob.d2tau0Avg = gam2*ob.d2tau0Avg + (1-gam2)*d2tau0;
            ob.d2kappa0Avg = gam2*ob.d2kappa0Avg + (1-gam2)*d2kappa0;  
            ob.d2h0xAvg = gam2*ob.d2h0xAvg + (1-gam2)*d2h0x;
            ob.d2h0yAvg = gam2*ob.d2h0yAvg + (1-gam2)*d2h0y;   
          end

          if fM.iter >= fM.fS.nWarmupIterations
            ob.deltaTau0 = fM.fS.momPrior*ob.deltaTau0 - ...
                           (fM.hessAlpha/ob.d2tau0Avg)*ob.dtau0Avg;
            ob.deltaKappa0 = fM.fS.momPrior*ob.deltaKappa0 - ...
                             (fM.hessAlpha/ob.d2kappa0Avg)*ob.dkappa0Avg;
%             ob.deltah0x = fM.fS.gradGain*ob.gradGainAniso*ob.dh0xAvg;
%             ob.deltah0y = fM.fS.gradGain*ob.gradGainAniso*ob.dh0yAvg;
            ob.deltah0x = fM.fS.momPrior*ob.deltah0x - ...
                           (fM.hessAlpha/ob.d2h0xAvg)*ob.dh0xAvg;
            ob.deltah0y = fM.fS.momPrior*ob.deltah0y - ...
                           (fM.hessAlpha/ob.d2h0yAvg)*ob.dh0yAvg;
          else      
            ob.deltaTau0 = fM.fS.gradGain*dtau0;
            ob.deltaKappa0 = fM.fS.gradGain*dkappa0;
%             ob.deltah0x = fM.fS.gradGain*ob.gradGainAniso*dh0x;
%             ob.deltah0y = fM.fS.gradGain*ob.gradGainAniso*dh0y;
            ob.deltah0x = fM.fS.gradGain*dh0x;
            ob.deltah0y = fM.fS.gradGain*dh0y;
          end
        else
          ob.deltaTau0 = fM.fS.gradGain*dtau0;
          ob.deltaKappa0 = fM.fS.gradGain*dkappa0;
%           ob.deltah0x = fM.fS.gradGain*ob.gradGainAniso*dh0x;
%           ob.deltah0y = fM.fS.gradGain*ob.gradGainAniso*dh0y;  
          ob.deltah0x = fM.fS.gradGain*dh0x;
          ob.deltah0y = fM.fS.gradGain*dh0y;        
        end
        if abs(ob.deltaTau0) > fM.fS.maxDelta;ob.deltaTau0 = ...
                                      fM.fS.maxDelta*sign(ob.deltaTau0);end
        if abs(ob.deltaKappa0) > fM.fS.maxDelta;ob.deltaKappa0 = ...
                                      fM.fS.maxDelta*sign(ob.deltaKappa0);end
        if abs(ob.deltah0x) > fM.fS.maxDelta;ob.deltah0x = ...
                                      fM.fS.maxDelta*sign(ob.deltah0x);end
        if abs(ob.deltah0y) > fM.fS.maxDelta;ob.deltah0y = ...
                                      fM.fS.maxDelta*sign(ob.deltah0y);end

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
        if ~ob.fixParams(3)
          h0x = h0x + ob.deltah0x; ob.hx = exp(h0x);
          ob.hx(ob.hx > ob.hxHB(:,2)) = ob.hxHB(ob.hx > ob.hxHB(:,2),2);
          ob.hx(ob.hx < ob.hxHB(:,1)) = ob.hxHB(ob.hx < ob.hxHB(:,1),1);
        end
        if ~ob.fixParams(4)
          h0y = h0y + ob.deltah0y; ob.hy = exp(h0y);
          ob.hy(ob.hy > ob.hyHB(:,2)) = ob.hyHB(ob.hy > ob.hyHB(:,2),2);
          ob.hy(ob.hy < ob.hyHB(:,1)) = ob.hyHB(ob.hy < ob.hyHB(:,1),1);
        end
      end
      
      % Save in vector
      ob.tau2Vec(fM.iter) = ob.tau2; ob.kappa2Vec(fM.iter) = ob.kappa2; 
      ob.hxVec(fM.iter) = ob.hx; ob.hyVec(fM.iter) = ob.hy; 
      
      % Store gradients
      ob.dtau0Vec(fM.iter) = dtau0';
      ob.dkappa0Vec(fM.iter) = dkappa0';
      ob.dh0xVec(fM.iter) = dh0x'; 
      ob.dh0yVec(fM.iter) = dh0y';
      ob.dtau0AvgVec(fM.iter) = ob.dtau0Avg;
      ob.dkappa0AvgVec(fM.iter) = ob.dkappa0Avg;
      ob.dh0xAvgVec(fM.iter) = ob.dh0xAvg; 
      ob.dh0yAvgVec(fM.iter) = ob.dh0yAvg;
      
      % Store Hessians
      if fM.fS.useHess
        ob.d2tau0Vec(fM.iter) = d2tau0';
        ob.d2kappa0Vec(fM.iter) = d2kappa0';
        ob.d2h0xVec(fM.iter) = d2h0x';
        ob.d2h0yVec(fM.iter) = d2h0y';
        ob.d2tau0AvgVec(fM.iter) = ob.d2tau0Avg;
        ob.d2kappa0AvgVec(fM.iter) = ob.d2kappa0Avg;
        ob.d2h0xAvgVec(fM.iter) = ob.d2h0xAvg;
        ob.d2h0yAvgVec(fM.iter) = ob.d2h0yAvg;
      end
      
      
    end    
    function doPolyakAveraging(ob,fM)
      ob.tau2 = mean(ob.tau2Vec(fM.iter-fM.fS.nPolyakVal+1:fM.iter));  
      ob.kappa2 = mean(ob.kappa2Vec(fM.iter-fM.fS.nPolyakVal+1:fM.iter));
      ob.hx = mean(ob.hxVec(fM.iter-fM.fS.nPolyakVal+1:fM.iter));    
      ob.hy = mean(ob.hyVec(fM.iter-fM.fS.nPolyakVal+1:fM.iter));    
    end
    function params = getParameters(ob)
      params = [ob.tau2,ob.kappa2,ob.hx,ob.hy];
    end
    function str = getParameterString(ob)
      str = ['tau2_',num2str(ob.k),': ',num2str(ob.tau2),', kappa2_',...
             num2str(ob.k),': ' ,num2str(ob.kappa2),...
             ' hx_',num2str(ob.k),': ',num2str(ob.hx),...
             ' hy_',num2str(ob.k),': ',num2str(ob.hy)];
    end
  end
end

