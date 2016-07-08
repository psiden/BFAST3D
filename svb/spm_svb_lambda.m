function [block] = spm_svb_lambda(Y,block)
% Modified version for SVB, Per Siden 2016-06-09
%
% Variational Bayes for GLM-AR models - Update lambda
% FORMAT [block] = spm_vb_lambda (Y,block)
%
% Y             [T x N] time series 
% block         data structure containing the following fields:
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_lambda.m 2451 2008-11-10 16:20:32Z lee $

if block.verbose
    disp('Updating lambda');
end

P=block.p;
K=block.k;
N=block.N;

if block.SVB.inWarmup
    Ns = block.SVB.NsWarmup;
else
    Ns = block.SVB.Ns;
end

if P <= 0
    
    GnEst = zeros(Ns,N);
    for i = 1:Ns
        wTemp = reshape(block.SVB.wSampSim(:,i)',K,N);
        GnEst(i,:) = block.SVB.YTY - 2*dot(block.SVB.YTX',wTemp) + dot(wTemp,block.SVB.XTX*wTemp);
    end
    Gn = mean(GnEst);
    block.b_lambda    = 1./(Gn'./2 + 1./block.b_lambda_prior);
    block.mean_lambda = block.c_lambda.*block.b_lambda;

else
        
    GnEst = zeros(Ns,N);
    for i = 1:Ns
        wTemp = reshape(block.SVB.wSampSim(:,i)',K,N);
        aTemp = reshape(block.SVB.aSampSim(:,i)',P,N);
        for p = 1:P
            for q = 1:P
                block.SVB.WSW(p,q,:) = dot(wTemp,squeeze(block.SVB.S(p,:,:,q)) * wTemp)';
                block.SVB.DW(p,q,:) = dot(squeeze(block.SVB.D(p,:,q,:)),wTemp);
            end
        end
        for p = 1:P
            block.SVB.RW(:,p,:) = permute(squeeze(block.SVB.R(:,:,p)) * wTemp,[1,3,2]);
        end
        for p = 1:P
            block.SVB.BW(p,:) = dot(squeeze(block.SVB.B(p,:,:)),wTemp);
            block.SVB.WRW(p,:) = dot(wTemp,squeeze(block.SVB.RW(:,p,:)));
        end
        if P <= 1
            block.SVB.ddTA = squeeze(block.SVB.ddT)' .* aTemp;
            block.SVB.DWA = squeeze(block.SVB.DW)' .* aTemp;
            block.SVB.WSWA = squeeze(block.SVB.WSW)' .* aTemp;
        else
            for p = 1:P
                block.SVB.ddTA(p,:) = dot(squeeze(block.SVB.ddT(p,:,:)),aTemp);
                block.SVB.DWA(p,:) = dot(squeeze(block.SVB.DW(p,:,:)),aTemp);
                block.SVB.WSWA(p,:) = dot(squeeze(block.SVB.WSW(p,:,:)),aTemp);
            end
        end
        GnEst(i,:) = block.SVB.YTY - 2*dot(block.SVB.YTX',wTemp) + dot(wTemp,block.SVB.XTX*wTemp) ...
                - 2*dot(block.SVB.YTdT',aTemp,1) + 2*dot(block.SVB.BW,aTemp,1) ...
                - 2*dot(block.SVB.WRW,aTemp,1) + dot(aTemp,block.SVB.ddTA,1) ...
                - 2*dot(aTemp,block.SVB.DWA,1) + dot(aTemp,block.SVB.WSWA,1);
    end
    Gn = mean(GnEst);
    block.b_lambda    = 1./(Gn'./2 + 1./block.b_lambda_prior);
    block.mean_lambda = block.c_lambda.*block.b_lambda;
        
end
