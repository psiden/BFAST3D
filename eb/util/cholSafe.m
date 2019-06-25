function R = cholSafe(A,maxeps)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Help function that calls Matlab's chol for sparse matrix A
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
% FIRST VER.:   2017-01-10
% REVISED:      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
try
    R = chol(A);
catch
    disp('Warning: chol failed,');
    epsA = 1e-15;
    tryCatchDone = 0;
    while epsA <= maxeps && ~tryCatchDone
        try
            R = chol(A + epsA*speye(size(A,1)));
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
end