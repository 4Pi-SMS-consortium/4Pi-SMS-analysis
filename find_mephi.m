function [mephi,uwmephi,mephipp]=find_mephi(dmap,mephi_ini,maxangle,srange,stopval)

[out1]=mephigrow(mephi_ini,dmap,maxangle,srange,stopval,1);
[out2]=mephigrow(mephi_ini([2 1],:),dmap,maxangle,srange,stopval,1);
if sum(out1(:,1))>sum(out2(:,1))
    [mephi]=cat(1,flipud(out2),out1(3:end,:));
else
    [mephi]=cat(1,flipud(out1),out2(3:end,:));
end
%     figure
%     scatter(mephi(:,1),mephi(:,2))
    mask=mephi(:,1)~=-pi&mephi(:,1)~=pi;
%
% uwmephi=mephi;
uwmephi(:,1)=mephi(mask,1);
uwmephi(:,2)=unwrap(mephi(mask,2));

mephi=mephi(mask,:);
mephipp = interp1(uwmephi(:,1),uwmephi(:,2),'linear','pp');
% mephipp=spline(uwmephi(:,1),uwmephi(:,2));
