function [G,Ce,C]=ndfem(FV,P,dirichlet)
% [G,Ce,C]=ndfem(FV,P,dirichlet)
%	Function that computes mass matrix and stiffness matrices for nD
%	triangulations
%
%	INPUTS:
%	FV		: Facet-vertex matrix 
%	P     : mesh nodes
%	dirichlet	: nodes with dirichlet boundary conditions
%	OUTPUTS:
%	G	: Stiffness matrix
%	Ce		: Mass matrix
%	C		: Mass lumped version of mass matrix.

d = size(FV,2)-1;
if(size(P,2)~=d)
  P = P';
end
nV = size(P,1); 
nF = size(FV,1); 
Gi = zeros(nF*(d+1),d+1); 
Gj = Gi; Gz = Gi; Ci = Gi; Cj = Gi; Cz = Gi;
    
for f=1:nF
  m = [ones(1,d+1);P(FV(f,:),:)'] \ [zeros(1,d);eye(d)];
  ddet = abs(det([ones(1,d+1);P(FV(f,:),:)'])); 
  
  Gi((d+1)*(f-1)+(1:d+1),:) = FV(f,:)'*ones(1,d+1);
  Gj((d+1)*(f-1)+(1:d+1),:) = ones(d+1,1)*FV(f,:);
  Gz((d+1)*(f-1)+(1:d+1),:) = ddet * (m*m') / prod(1:d);
  
  Ci((d+1)*(f-1)+(1:d+1),:) = FV(f,:)'*ones(1,d+1);
  Cj((d+1)*(f-1)+(1:d+1),:) = ones(d+1,1)*FV(f,:);
  Cz((d+1)*(f-1)+(1:d+1),:) = ddet*(ones(d+1)+eye(d+1)) ./ prod(1:d+2);
  
end
G = sparse(Gi(:),Gj(:),Gz(:),nV,nV);
Ce = sparse(Ci(:),Cj(:),Cz(:),nV,nV);
if(nargout>2)
  C = spdiags(sum(Ce)',0,nV,nV);
end
if (nargin>2)
  G = G(~dirichlet,~dirichlet);
  Ce = Ce(~dirichlet,~dirichlet);
  if(nargout>2)
    C = C(~dirichlet,~dirichlet);  
  end
end