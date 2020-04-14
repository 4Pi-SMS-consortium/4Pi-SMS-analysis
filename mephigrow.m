function out=mephigrow(mephi,dmap,maxangle,srange,stopval,invflag)

I=1e37;
msz=size(dmap,1);
while I>stopval
    v1=(mephi(end,1)-mephi(end-1,1))+wrapToPi(mephi(end,2)-mephi(end-1,2))*1i;
    xx=linspace(-pi,pi,msz);
    yy=linspace(-pi,pi,msz);
    [xgrid ygrid]=meshgrid(xx,yy);
    rmap=sqrt((xgrid-mephi(end,1)).^2+wrapToPi(ygrid-mephi(end,2)).^2);
    anglemap=wrapToPi(acos(((mephi(end,1)-mephi(end-1,1)).*(xgrid-mephi(end,1))+...
        wrapToPi(mephi(end,2)-mephi(end-1,2))*wrapToPi(ygrid-mephi(end,2)))./abs(v1)./abs(rmap)));
    
    maska=anglemap<maxangle;
    maskr=rmap<srange(2)&rmap>srange(1);
    directmask=(xgrid-mephi(end,1))./(mephi(end,1)-mephi(end-1,1))>=invflag.*0.5;    
    mask=maska.*maskr.*directmask;
    
    vismap=mask.*dmap;
    [I row col]=findmax(vismap);
    mephi(end+1,2)=yy(row);
    mephi(end,1)=xx(col);
    mephi(end,3)=I;
end

out=mephi;