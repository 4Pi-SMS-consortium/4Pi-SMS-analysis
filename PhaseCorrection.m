function [zangresult,mtcresult]=PhaseCorrection(zangresult,mtcresult,tresult,stacktot,fnum,reverseflag)
% fnum=1500;
frmnum=fnum;
msz=256;
cutmeth='nocut';
pixelsz=1; % nm
thresh=2; % nm

order=(1:length(stacktot))';
Tout=tresult;
unistack=unique(stacktot);
L=numel(unistack) ;
EA={};

for k=1:L
    k
    close all
    currst=unistack(k);
    maskst=(stacktot==currst);
    zangresult1=zangresult(maskst);
    mtcresult1=mtcresult(maskst);
    tt=Tout(maskst);
    cycle=floor(tt/(frmnum*numel(unistack)));
    remn=rem(tt,frmnum*numel(unistack));
    tout=cycle*frmnum+remn-(k-1)*frmnum;
    OD=order(maskst);
    E=[];
    i=0;
%     if k==1 
%         thre1=0.1;
%         thre2=0.1;
%     else
%         thre1=0.05;
%         thre2=0.05;    
%     end
    while 1
        i=i+1;
        xout=(zangresult1+pi).*msz./2./pi;
        yout=(mtcresult1+pi)*msz./2./pi;
        [shiftx,shifty]=correct_drift_LS(single(xout),single(yout),tout,frmnum,pixelsz,thresh,cutmeth);
%         x=xout*10;
%         y=yout*10;
%         z=yout*10;
%         frame=tout;
%         p.framestart=0;
%         p.framestop=ceil((max(frame)/frmnum))*frmnum-1;
%         points=ceil((max(frame)/frmnum));
%         p.correctxy=true;
%         p.correctz=false;
%         p.drift_timepoints=points;
%         p.drift_maxdrift=500;
%         close all
%         [drift,driftinfo,fieldc]=driftcorrection3D_so(x,y,z,frame,p);
%         shiftx=driftinfo.xy.dx;
%         shifty=driftinfo.xy.dy;
%         shiftx=diff(shiftx)/10;
%         shifty=diff(shifty)/10;    

        E(i,1)=sum(shiftx.^2);
        E(i,2)=sum(shifty.^2);
        disp(sum(shiftx.^2))
        disp(sum(shifty.^2))
        close all; figure; plot(shiftx,'b-'); hold on; plot(shifty,'g-'); pause(0.5);
        
        if i>1
%             if (E(i,1)>=E(i-1,1) && E(i,2)>=E(i-1,2) && E(i-1,1)<thre1 && E(i-1,2)<thre2) || i>100 %|| (E(i,1)<1 && E(i,2)<0.1)
            if (max(abs(shiftx))<1 && max(abs(shifty))<0.5)|| i>50 %&& (E(i,1)>=E(i-1,1) && E(i,2)>=E(i-1,2)) %|| i>100 
                i
                break
            end
        end
        [xout2]=shiftcoords_LS(xout,shiftx,tout,frmnum,reverseflag);
        [yout2]=shiftcoords_LS(yout,shifty,tout,frmnum,reverseflag);
        zangresult1=xout2*2*pi/msz-pi;
        mtcresult1=yout2*2*pi/msz-pi;
        zangresult1=wrapToPi(zangresult1);
    end
    zangresult(OD)=zangresult1;
    mtcresult(OD)=mtcresult1;
    E
    EA{1}=E;
end
% str='E:\4Pi_two_color\2018-3-2\Pahse_E.mat';
% save(str,'EA');

% segnum=ceil((max(tresult))/fnum);
% zest=[];
% zerr=[];
% zmask=[];
% sigma=2;
% msz=256;
% for ii=1:1:segnum
%     st=(ii-1)*fnum;
%     if ii==segnum
%         ed=max(tresult);
%     else
%         ed=(ii)*fnum-1;
%     end
%     maskt=tresult>=st&tresult<=ed;
%     mtc_seg=mtcresult(maskt);
%     zang_seg=zangresult(maskt);    
% %     tmpmask=mtc_seg>-0.05&mtc_seg<0.05;
% %     tmp=zang_seg(tmpmask);
% %     [tmphis xout]=hist(tmp,64);
% %     inis=[max(tmphis(:)) mean(tmp(:)) 0.5];
% %     [c]=fit1dgaussian(tmphis,xout,inis);
% %     phi_0(ii,1)=c(2);        
%     [dmap]=build_dmap(mtc_seg,zang_seg,msz,sigma); 
%     tmp=mean(dmap(:,128:129),2);
%     xout=(1:256)';
%     inis=[max(tmp(:)) 128 20];
%     [c]=fit1dgaussian(tmp,xout,inis);
%     phi_0(ii,1)=c(2); 
% end
% 
% phi_0=phi_0-phi_0(end);
% phi_0=phi_0*2*pi/256;
% 
% for ii=1:1:segnum
%     st=(ii-1)*fnum;
%     if ii==segnum
%         ed=max(tresult);
%     else
%         ed=(ii)*fnum-1;
%     end
%     maskt=tresult>=st&tresult<=ed;
% %     mtc_seg=mtcresult(maskt);
% %     zang_seg=zangresult(maskt);   
%     zangresult(maskt)=zangresult(maskt)-phi_0(ii);
% end
% zangresult=wrapToPi(zangresult);