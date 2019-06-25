function x = plot3d(clplot,fnbr,snbr)

if nargin < 2 
	figure
else
	figure(fnbr)
end
if nargin ==3
	subplot(snbr)
end

[m n l] = size(clplot);

[XX YY ZZ] = meshgrid(1:m,1:n,1:l);
zmid = round(l/2);
ymid = round(n/2);
xmid = round(m/2);


surf(squeeze(XX(:,:,zmid)),squeeze(YY(:,:,zmid)),squeeze(ZZ(:,:,zmid)),...
		squeeze(clplot(:,:,zmid)),'EdgeColor','none')
hold on
surf(squeeze(XX(:,ymid,:)),squeeze(YY(:,ymid,:)),squeeze(ZZ(:,ymid,:)),...
		squeeze(clplot(:,ymid,:)),'EdgeColor','none')
surf(squeeze(XX(xmid,:,:)),squeeze(YY(xmid,:,:)),squeeze(ZZ(xmid,:,:)),...
		squeeze(clplot(xmid,:,:)),'EdgeColor','none')
axis image
x = 0;

ind = ~isnan(clplot(:,ymid,zmid));
ymin = min(YY(ind,ymid,zmid));
ymax = max(YY(ind,ymid,zmid));
plot3([ymid ymid],[ymin ymax], [zmid zmid],'k')

ind = ~isnan(clplot(xmid,:,zmid));
xmin = min(XX(xmid,ind,zmid));
xmax = max(XX(xmid,ind,zmid));
plot3([xmin xmax],[ymid ymid], [zmid zmid],'k')

ind = ~isnan(clplot(xmid,ymid,:));
zmin = min(ZZ(xmid,ymid,ind));
zmax = max(ZZ(xmid,ymid,ind));
plot3([xmid xmid],[ymid ymid], [zmin zmax],'k')


