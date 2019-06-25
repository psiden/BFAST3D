%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Setup precision matrices
%               Q = G'*G
%               G not properly updated for all dimensions/types
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University
%
% FIRST VER.:   2016-01-15 (replaces setupQs.m)
% REVISED:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [QxList,QyList,QzList,GxList,GyList,GzList] = setupAnisoPrecMats(K,N,sz,bmask,ndim)

QxList = cell(K,1);
GxList = cell(K,1);
QyList = cell(K,1);
GyList = cell(K,1);
QzList = cell(K,1);
GzList = cell(K,1);

for k = 1:K
  if ndim == 2
    Dr = spdiags([-1*ones(sz(1)-1,1),1*ones(sz(1)-1,1)],[0,1],sz(1)-1,sz(1));
    Dc =  spdiags([-1*ones(sz(2)-1,1),1*ones(sz(2)-1,1)],[0,1],sz(2)-1,sz(2));
    Drs = kron(speye(sz(2)),Dr);
    Dcs = kron(Dc,speye(sz(1)));
    bmask_mat = zeros(sz); bmask_mat(bmask) = 1;
    hasColNeigh = bmask_mat(1:(end-1),:) & ~diff(bmask_mat);
    hasRowNeigh = bmask_mat(:,1:(end-1)) & ~diff(bmask_mat')';
    GxList{k} = Drs(hasColNeigh(:),bmask);
    QxList{k} = GxList{k}'*GxList{k};
    GyList{k} = Dcs(hasRowNeigh(:),bmask);
    QyList{k} = GyList{k}'*GyList{k};
    GzList{k} = 0;
    QzList{k} = 0;
  elseif ndim == 3
    Dx = spdiags([-1*ones(sz(1)-1,1),1*ones(sz(1)-1,1)],[0,1],sz(1)-1,sz(1));
    Dy = spdiags([-1*ones(sz(2)-1,1),1*ones(sz(2)-1,1)],[0,1],sz(2)-1,sz(2));
    Dz = spdiags([-1*ones(sz(3)-1,1),1*ones(sz(3)-1,1)],[0,1],sz(3)-1,sz(3));
    Dxs = kron(speye(sz(3)),kron(speye(sz(2)),Dx));
    Dys = kron(speye(sz(3)),kron(Dy,speye(sz(1))));
    Dzs = kron(Dz,kron(speye(sz(2)),speye(sz(1))));
    bmask_mat = zeros(sz); bmask_mat(bmask) = 1;
    hasXNeigh = bmask_mat(1:(end-1),:,:) & ~diff(bmask_mat);
    hasYNeigh = bmask_mat(:,1:(end-1),:) & ~diff(bmask_mat,1,2);
    hasZNeigh = bmask_mat(:,:,1:(end-1)) & ~diff(bmask_mat,1,3);
    GxList{k} = Dxs(hasXNeigh(:),bmask);
    QxList{k} = GxList{k}'*GxList{k};
    GyList{k} = Dys(hasYNeigh(:),bmask);
    QyList{k} = GyList{k}'*GyList{k};
    GzList{k} = Dzs(hasZNeigh(:),bmask);
    QzList{k} = GzList{k}'*GzList{k};
  end
end


