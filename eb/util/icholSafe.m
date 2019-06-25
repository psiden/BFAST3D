function [R,epsA] = icholSafe(A,starteps,maxeps)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Help function that calls Matlab's ichol for sparse matrix A
%               in try/catch to
%               prevent numerical errors when the matrix is just close to
%               positive definite. If chol fails, retry after adding maxeps
%               (or less) to the diagonal.
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2017-06-30
% REVISED:      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
try
    R = ichol(A);
    epsA = starteps;
catch
    disp('Warning: ichol failed,');
    tryCatchDone = 0;
    epsA = starteps;
    while epsA <= maxeps && ~tryCatchDone
        try
            R = ichol(A + epsA*spdiags(diag(A),0,size(A,1),size(A,1)));
            tryCatchDone = 1;            
        catch
            epsA = epsA * 10;
        end
        if epsA > maxeps
            disp(['and adding a diagonal eps <= ',num2str(maxeps),' did not help.']);
            return
        end
    end
    disp(['but worked after adding a diagonal eps = ',num2str(epsA),'.']);
    epsA = .1*epsA;    
end
