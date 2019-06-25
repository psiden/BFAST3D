function x = parallelPCG(Q,b,tol,maxiter,icholQ,xStart)
%%

N = size(Q,1);
Ns = size(xStart,2);
x = zeros(N,Ns);
iterSave = zeros(1,Ns);
parfor i = 1:Ns
    [xTemp,flag,relres,iter] = pcg(Q,b(:,i),tol,maxiter,icholQ,icholQ',xStart(:,i));
    iterSave(i) = iter;
    x(:,i) = xTemp;
end