function maskout=reject_beads(V,fnum,sz)

V(:,4)=1;
L=ceil(max(V(:,3))/fnum);
for i=1:L
    st=(i-1)*fnum+1;
    et=i*fnum;
    id=V(:,3)>=st&V(:,3)<=et;
    V1=V(id,:);
    im=cHistRecon(sz,sz,V1(:,1),V1(:,2),0);
    [cx,cy]=find(im>300);
    l=length(cx);
    for j=1:l
        xx=cx(j);
        yy=cy(j);
        dist=sqrt((V1(:,1)-xx).^2+(V1(:,2)-yy).^2);
        ix=dist<3;
        V1(ix,4)=0;
    end
    V(id,4)=V1(:,4);
end

maskout=V(:,4);