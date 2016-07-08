%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Sample from GMRF in parallel
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [samples,iterSave] = spm_svb_parallel_GMRF_sampling(Q,b,tol,icholQ,xStart)
%%

n = size(Q,1);
Ns = size(xStart,2);
samples = zeros(n,Ns);
iterSave = zeros(1,Ns);
parfor i = 1:Ns
    [x,flag,relres,iter] = pcg(Q,b(:,i),tol,500,icholQ,icholQ',xStart(:,i));
    iterSave(i) = iter;
    samples(:,i) = x;
end