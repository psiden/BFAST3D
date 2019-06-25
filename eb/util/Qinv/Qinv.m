function iQ = Qinv(R,nthreads)
%Qinv - Compute sparse elements of inv(Q).
%
% iQ = Qinv(R,nthreads=8)
%
% Here R=chol(Q) and nthreads is the number of threads to be used by openMP
% in the C+ implementation. For small matrices, selecting 1 thread might be
% faster due to the openMP overhead.
%
% Uses precompiled mex-file if possible.

% $Id: Qinv.m 2759 2011-09-07 12:35:07Z johanl $

warning('Using matlab code, not pre-compiled C-code');

%compute the D and U matrices
D = diag(R); 
if nnz(R(:,end))==1
  %assume lower triangular matrix
  U = R*spdiags(1./D,0,size(R,1),size(R,1));
else
  U = spdiags(1./D,0,size(R,1),size(R,1))*R;
  U=U';
end
D = 1./(D.*D);

%find number of non-zero elements in the final matrix
N = size(R,1);

%create initial iQ matrix using the lower right element.
iQ = sparse(N, N, D(N));

%loop over the columns of the matrix
for i=(N-1):-1:1
  %find non-zero elements in Uik
  Uik = U(:,i)';
  j = find(Uik);
  j = j(2:end);
  %extract the relevant parts of U and Z
  %compute sum(Uik*Zkj) for the non-zero elements
  UZ = Uik*iQ(:,j);
  Uik = Uik(j);
  %add these elements to iQ
  i2 = i*ones(1,length(j));
  ii = [i i2 j];
  jj = [i j i2];
  s = [D(i) + sum(UZ.*Uik) -UZ -UZ];
  %Add the new elements to the existing iQ-matrix
  iQ = iQ + sparse(ii,jj,s,N,N);
end
