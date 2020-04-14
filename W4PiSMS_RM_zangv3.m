function [z_ang,ang_contrast]=W4PiSMS_RM_zangv3(anglefile,q1,q2,q3,q4,locmaxc_f,subsz,xf,yf,offsetim,varim,gainim,fmfile,tlz_f,bg_f)
ldtmp=load(fmfile);
ldtmp.R(:,:,1)=eye(3);
offsetim=single(repmat(permute(offsetim,[1 2 4 3]),[1 1 size(q1,3) 1]));
varim=single(repmat(permute(varim,[1 2 4 3]),[1 1 size(q1,3) 1]));
gainim=single(repmat(permute(gainim,[1 2 4 3]),[1 1 size(q1,3) 1]));

xc=double(xf+tlz_f(:,2));
yc=double(yf+tlz_f(:,1));
tc=double(tlz_f(:,3));

L=length(xf);
subims=zeros(subsz,subsz,L,4);
xcenter=zeros(L,4);
ycenter=zeros(L,4);
suboffset=zeros(subsz,subsz,L,4);
subvar=zeros(subsz,subsz,L,4);
subgain=zeros(subsz,subsz,L,4);

for ii=1:1:4
    eval(['im=q' num2str(ii) ';']);
    cocurr=[xc' ; yc'; ones(1,numel(xc))];
    transco=ldtmp.R(:,:,ii)*cocurr;
    trans_x=transco(1,:);
    trans_y=transco(2,:);
    [subims(:,:,:,ii),t,l]=cMakeSubregions(round(trans_y(:)),round(trans_x(:)),tc,subsz,single(permute(im,[1 2 3])));
    xcenter(:,ii)=trans_x(:)-l;
    ycenter(:,ii)=trans_y(:)-t;
    
    % crop other regions
    [suboffset(:,:,:,ii)]=cMakeSubregions(round(trans_y(:)),round(trans_x(:)),tc,subsz,single(permute(offsetim(:,:,:,ii),[1 2 3])));  % cmake subregions should crop noise maps
    [subvar(:,:,:,ii)]=cMakeSubregions(round(trans_y(:)),round(trans_x(:)),tc,subsz,single(permute(varim(:,:,:,ii),[1 2 3])));  % cmake subregions should crop noise maps
    [subgain(:,:,:,ii)]=cMakeSubregions(round(trans_y(:)),round(trans_x(:)),tc,subsz,single(permute(gainim(:,:,:,ii),[1 2 3])));  % cmake subregions should crop noise maps
end

sigma=0.9;
sz=size(subims,1);
Model=single(zeros(sz,sz,size(subims,3),4));
for ii=1:4
    yf=ycenter(:,ii);
    xf=xcenter(:,ii);
    ROI=finiteGaussPSFerf(sz,sigma,1,0,[xf,yf]);
    Model(:,:,:,ii)=single(dip_array(ROI));
end

int_est=zeros(size(subims,3),size(subims,4));
for ii=1:1:size(subims,3)
    for jj=1:1:size(subims,4)
        model=Model(:,:,ii,jj);
        varim=subvar(:,:,ii,jj);
        gainim=subgain(:,:,ii,jj);
        subim=subims(:,:,ii,jj);
        bg=bg_f(ii);
        var_all=subim+varim./gainim./gainim;
        nom=(subim-bg./4).*model./var_all;
        denom=model.^2./var_all;
        int_est(ii,jj)=sum(nom(:))./sum(denom(:));                        
    end
    if mod(ii,10000)==0
        disp([num2str(ii) ' out of ' num2str(size(subims,3)) ' is done...']);
    end
end

rmp=(int_est(:,1)-int_est(:,3))./(int_est(:,1)+int_est(:,3));
rms=(int_est(:,4)-int_est(:,2))./(int_est(:,4)+int_est(:,2));

tmpload=load(anglefile);
ang2=MyCalPahse_new(rms,rmp,tmpload.phi_s,tmpload.phi_p);
z_ang=wrapToPi(ang2(:,2));
ang_contrast=ang2(:,1);

id=ang_contrast<0 | ang_contrast>sqrt(2) | abs(rms)>1 | abs(rmp)>1;
z_ang(id)=0;
ang_contrast(id)=NaN;

figure;
scatter(rms(~id),rmp(~id),5,ang_contrast(~id));
pause(eps)