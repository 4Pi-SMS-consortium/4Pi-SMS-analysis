function [mephi_ini,cpeak]=find_mephi_ini(dmap,srange,msz)

% normalization helps to emphasize on center map
normdmap=double(dmap.*gaussianblob(newim(size(dmap)),[size(dmap,1)/2 size(dmap,2)/2],size(dmap,1)/4,1));
[I,row,col]=findmax(normdmap);
cpeak=dmap(row,col);
[mephi_ini]=coords2mephi([col row],msz);
mephi_ini(1,3)=I;

xx=linspace(-pi,pi,msz);
yy=linspace(-pi,pi,msz);
[xgrid ygrid]=meshgrid(xx,yy);
rmap=sqrt((xgrid-mephi_ini(end,1)).^2+wrapToPi(ygrid-mephi_ini(end,2)).^2);
anglemap=angle((xgrid-mephi_ini(end,1))+(ygrid-mephi_ini(end,2))*1i);

anglemask=anglemap<-pi./4&anglemap>-pi.*7/16;
directmask=(xgrid-mephi_ini(end,1))>=0.01;
maskr=rmap<(srange(2)./3)&rmap>(srange(1)./3);
% dipshow(maskr.*dmap.*directmask.*anglemask);
[I,row,col]=findmax(maskr.*dmap.*directmask.*anglemask);
mephi_ini(end+1,2)=yy(row);
mephi_ini(end,1)=xx(col);
mephi_ini(end,3)=I;