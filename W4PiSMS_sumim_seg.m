function [subims tlz locmaxc subvar_g]=W4PiSMS_sumim_seg(im,thresh,subsz,varim,sumvg,pick_flag,fit_flag)

if fit_flag==2
    sz2=5;
else
    sz2=9;
end

%% modified Apr27,2016
fcin=im;
imsz=size(fcin,1);
varim=repmat(varim,[1 1 size(fcin,3)]);
sz=5;
[filteredim1]=varunif(squeeze(fcin),squeeze(varim),sz);
sz=9;
[filteredim2]=varunif(squeeze(filteredim1),squeeze(varim),sz);
im_unif=filteredim1-filteredim2;
L=size(im_unif,3);
for i=1:L
    A=im_unif(:,:,i);
    deta=median(median(abs(A-median(median(A)))));
    th=thresh*deta/0.67;
    A=(A>=th).*A;
    im_unif(:,:,i)=A;
end
sz=sz2;
im_max=(im_unif==maxf(im_unif,[sz sz 0],'rectangular'))&(im_unif>0);
a=find(im_max);
z=floor(a/imsz/imsz);
pnum=mod(a,imsz*imsz);
y=mod(pnum,imsz);
x=floor(pnum/imsz);
locmaxc=[x y z];

%%
r=(subsz-1)/2;
rangemin=[r-1 r-1];
rangemax=[imsz-r-1 imsz-r-1];
if pick_flag==1
    mask=locmaxc(:,1)<rangemax(1)&locmaxc(:,1)>rangemin(1)...
        &locmaxc(:,2)<rangemax(2)&locmaxc(:,2)>rangemin(2);
    locmaxc=locmaxc(mask,:);
end

if isempty(locmaxc)
    subims=[];
    tlz=[];
else
    [subims t l]=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(im,[1 2 3])));
    [subvar_g]=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(repmat(permute(sumvg,[1 2 3]),[1 1 size(im,3)])));
    tlz=[t l locmaxc(:,3)];
end