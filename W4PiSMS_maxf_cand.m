function coords=W4PiSMS_maxf_cand(imunif,sz,thresh)

im_max=(imunif>=maxf(imunif,[sz sz 0],'rectangular'))&(imunif>thresh);
imsz=size(imunif,1);
a=find(im_max);
z=floor(a/imsz/imsz);
pnum=mod(a,imsz*imsz);
y=mod(pnum,imsz);
x=floor(pnum/imsz);
coords=[x y z];