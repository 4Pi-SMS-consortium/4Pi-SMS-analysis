function [zest_f,zerr_f mephimask]=mephi_zest(mephipp,uwmephi,mtc_seg,zang_seg,lambda,Nmedia,centermtc,freq)
% center phi_0 can not be obtained by option 'cofocal' which will perform a
% estimate of the phi_0 based on the zang_seg data which is more
% accurate than mephipp.

% center phi_0 is the distribution center of the mephi plot
% phi_0=ppval(mephipp,(max(uwmephi(:,1))+min(uwmephi(:,1)))/2); %phi_mid can be set
if isempty(centermtc)
    centermtc=0;
    phi_0=ppval(mephipp,centermtc);
%     phi_0=wrapToPi(phi_0);
elseif strcmp(class(centermtc),'double')==1
    tmpmask=mtc_seg>(centermtc-0.02)&mtc_seg<(centermtc+0.02);
    tmp=zang_seg(tmpmask);
    [tmphis xout]=hist(tmp,50);
    inis=[max(tmphis(:)) mean(tmp(:)) 0.5];
    [c]=fit1dgaussian(tmphis,xout,inis);
    phi_0=c(2);
elseif strcmp(centermtc,'mid')
    [nout xout]=hist(mtc_seg,20);
    [I ind]=max(nout);
    centermtc2=xout(ind);
    phi_0=ppval(mephipp,centermtc2);   
elseif strcmp(centermtc,'cofocal')
    tmpmask=mtc_seg>-0.02&mtc_seg<0.02;
    tmp=zang_seg(tmpmask);
    [tmphis xout]=hist(tmp,50);
    inis=[max(tmphis(:)) mean(tmp(:)) 0.5];
    [c]=fit1dgaussian(tmphis,xout,inis);
    phi_0=c(2);
end

% if phi_0>=2*pi||phi_0<=-2*pi
%     display(['The phi_0 estimation did not converge, setting this value to zero;']);
%     phi_0=0;
% end

% phi_0=ppval(mephipp,centermtc); %phi_mid set to 0 for segment align

display(['phi_0 is ' num2str(phi_0)]);

% tmpld=load('frequency_simulate4Pi_570_671.mat');

tmp=-10*2*pi:2*pi:10*2*pi;
pispace=repmat(tmp,[size(zang_seg,1) 1]);
angle1space=repmat(zang_seg,[1 size(pispace,2)]);
angle1span=(angle1space+pispace);
phispan=ppval(mephipp,mtc_seg);
phispan=repmat(phispan(:),[1 size(pispace,2)]);
phidiff1=abs(phispan-angle1span);

phidiff=phidiff1;
[C, I]=min(phidiff,[],2);
mindiffspan=repmat(C,[1 size(pispace,2)]);
mask=(mindiffspan==phidiff);
uwangle=sum((double(mask).*angle1span),2);
angle_error=C;

% figure;hold on
% scatter(mtc_seg,uwangle,9)
% plot(-pi:0.1:pi,ppval(mephipp,-pi:0.1:pi),'r--','linewidth',6)

norm_uwangle=(uwangle-phi_0);

% f=2.*pi./(lambda./2./Nmedia);
% global freqstr;
% eval(['f=tmpld.' freqstr ';']);
f = freq;
zest=norm_uwangle./f;
zerr=angle_error./f;

%filtering
mephimask=mtc_seg<max(uwmephi(:,1))&mtc_seg>min(uwmephi(:,1));

zest_f=zest;
zerr_f=zerr;
