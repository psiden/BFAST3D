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
% FIRST VER.:   2016-06-09
% REVISED:      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [QList,GList] = spm_svb_setupPrecMats(QTypes,N,sz,bmask,ndim)

K = length(QTypes);
QList = cell(K,1);
GList = cell(K,1);

for k = 1:K
    if strcmp(QTypes{k},'eye') % Identity
        GList{k} = speye(N);
        QList{k} = speye(N);
        
    elseif strcmp(QTypes{k},'L') % Laplacian
        if ndim == 2
            Rr = spdiags([-1*ones(sz(1),1),2*ones(sz(1),1),-1*ones(sz(1),1)],[-1,0,1],sz(1),sz(1));
            Rc = spdiags([-1*ones(sz(2),1),2*ones(sz(2),1),-1*ones(sz(2),1)],[-1,0,1],sz(2),sz(2));
            Limage = kron(Rc,spdiags(ones(sz(1),1),0,sz(1),sz(1))) + ...
                     kron(spdiags(ones(sz(2),1),0,sz(2),sz(2)),Rr);
            L = Limage(bmask,bmask);
            QList{k} = L;
            GList{k} = '';
        elseif ndim == 3
            Rx = spdiags([-1*ones(sz(1),1),2*ones(sz(1),1),-1*ones(sz(1),1)],[-1,0,1],sz(1),sz(1));
            Ry = spdiags([-1*ones(sz(2),1),2*ones(sz(2),1),-1*ones(sz(2),1)],[-1,0,1],sz(2),sz(2));
            Rz = spdiags([-1*ones(sz(3),1),2*ones(sz(3),1),-1*ones(sz(3),1)],[-1,0,1],sz(3),sz(3));
            Lxy = kron(Ry,spdiags(ones(sz(1),1),0,sz(1),sz(1))) + ...
                  kron(spdiags(ones(sz(2),1),0,sz(2),sz(2)),Rx);
            Limage = kron(Rz,spdiags(ones(sz(1)*sz(2),1),0,sz(1)*sz(2),sz(1)*sz(2))) + ...
                     kron(spdiags(ones(sz(3),1),0,sz(3),sz(3)),Lxy);
            L = Limage(bmask,bmask);
            QList{k} = L;
            GList{k} = '';
        end
        
    elseif strcmp(QTypes{k},'L2') % squared Laplacian
        if ndim == 2
            Rr = spdiags([-1*ones(sz(1),1),2*ones(sz(1),1),-1*ones(sz(1),1)],[-1,0,1],sz(1),sz(1));
            Rc = spdiags([-1*ones(sz(2),1),2*ones(sz(2),1),-1*ones(sz(2),1)],[-1,0,1],sz(2),sz(2));
            Limage = kron(Rc,spdiags(ones(sz(1),1),0,sz(1),sz(1))) + ...
            kron(spdiags(ones(sz(2),1),0,sz(2),sz(2)),Rr);
            L = Limage(bmask,bmask);
            QList{k} = L'*L;
            GList{k} = L;
        elseif ndim == 3
            Rx = spdiags([-1*ones(sz(1),1),2*ones(sz(1),1),-1*ones(sz(1),1)],[-1,0,1],sz(1),sz(1));
            Ry = spdiags([-1*ones(sz(2),1),2*ones(sz(2),1),-1*ones(sz(2),1)],[-1,0,1],sz(2),sz(2));
            Rz = spdiags([-1*ones(sz(3),1),2*ones(sz(3),1),-1*ones(sz(3),1)],[-1,0,1],sz(3),sz(3));
            Lxy = kron(Ry,spdiags(ones(sz(1),1),0,sz(1),sz(1))) + ...
                  kron(spdiags(ones(sz(2),1),0,sz(2),sz(2)),Rx);
            Limage = kron(Rz,spdiags(ones(sz(1)*sz(2),1),0,sz(1)*sz(2),sz(1)*sz(2))) + ...
                     kron(spdiags(ones(sz(3),1),0,sz(3),sz(3)),Lxy);
            L = Limage(bmask,bmask);
            QList{k} = L'*L;
            GList{k} = L;
        end
        
    elseif strcmp(QTypes{k},'LI') % intrinsic Laplacian
        if ndim == 2
            Dr = spdiags([-1*ones(sz(1)-1,1),1*ones(sz(1)-1,1)],[0,1],sz(1)-1,sz(1));
            Dc =  spdiags([-1*ones(sz(2)-1,1),1*ones(sz(2)-1,1)],[0,1],sz(2)-1,sz(2));
            Drs = kron(speye(sz(2)),Dr);
            Dcs = kron(Dc,speye(sz(1)));
            bmask_mat = zeros(sz); bmask_mat(bmask) = 1;
            hasColNeigh = bmask_mat(1:(end-1),:) & ~diff(bmask_mat);
            hasRowNeigh = bmask_mat(:,1:(end-1)) & ~diff(bmask_mat')';
            GList{k} = [Drs(hasColNeigh(:),bmask);Dcs(hasRowNeigh(:),bmask)];
            QList{k} = GList{k}'*GList{k};
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
            GList{k} = [Dxs(hasXNeigh(:),bmask);Dys(hasYNeigh(:),bmask);Dzs(hasZNeigh(:),bmask)];
            QList{k} = GList{k}'*GList{k};
        end
        
    elseif strcmp(QTypes{k},'L2I') % intrinsic squared Laplacian
        if ndim == 2
            Rr = spdiags([-1*ones(sz(1),1),2*ones(sz(1),1),-1*ones(sz(1),1)],[-1,0,1],sz(1),sz(1));
            Rc = spdiags([-1*ones(sz(2),1),2*ones(sz(2),1),-1*ones(sz(2),1)],[-1,0,1],sz(2),sz(2));
            Limage = kron(Rc,spdiags(ones(sz(1),1),0,sz(1),sz(1))) + ...
                     kron(spdiags(ones(sz(2),1),0,sz(2),sz(2)),Rr);
            L = Limage(bmask,bmask);
            L = L - diag(sum(L));
            QList{k} = L'*L;
            GList{k} = L;
        elseif ndim == 3
            Rx = spdiags([-1*ones(sz(1),1),2*ones(sz(1),1),-1*ones(sz(1),1)],[-1,0,1],sz(1),sz(1));
            Ry = spdiags([-1*ones(sz(2),1),2*ones(sz(2),1),-1*ones(sz(2),1)],[-1,0,1],sz(2),sz(2));
            Rz = spdiags([-1*ones(sz(3),1),2*ones(sz(3),1),-1*ones(sz(3),1)],[-1,0,1],sz(3),sz(3));
            Lxy = kron(Ry,spdiags(ones(sz(1),1),0,sz(1),sz(1))) + ...
                  kron(spdiags(ones(sz(2),1),0,sz(2),sz(2)),Rx);
            Limage = kron(Rz,spdiags(ones(sz(1)*sz(2),1),0,sz(1)*sz(2),sz(1)*sz(2))) + ...
                     kron(spdiags(ones(sz(3),1),0,sz(3),sz(3)),Lxy);
            L = Limage(bmask,bmask);
            L = L - diag(sum(L));
            QList{k} = L'*L;
            GList{k} = L;
        end
    else
        display('ERROR: Invalid Precision Matrix Type!');
    end
    
end


