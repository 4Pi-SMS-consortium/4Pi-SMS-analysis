function [zest zerr zmask]=Mephi_z_4PiSMS(zangresult,mtcresult,tresult,fnum,mephiBin,stopval,centermtc,freq)
% setup
stopval0=stopval;
msz=256;
srange=[0.4 0.6];
maxangle=0.2;
% stopval=0.03;
lambda=680;
Nmedia=1.33;
sigma=5; % from 10, for Spermatocytes to 5 for cy3b
%
% fnum=max(tresult)+1;
% fnum=fnum*2;
% fnum=3000;
fnum = fnum*mephiBin;
segnum=ceil((max(tresult))/fnum);
zest=[];
zerr=[];
zmask=[];
h=figure(1);
ii =1;
while ii<=segnum
% for ii=1:1:segnum
%     if ii<=38
%         stopval=0.15;
%     else
%         stopval=0.25;
%     end
    st=(ii-1)*fnum;
    if ii==segnum
        ed=max(tresult);
    else
        ed=(ii)*fnum-1;
    end
    
    maskt=tresult>=st&tresult<=ed;
    mtc_seg=mtcresult(maskt);
    %     currzast_err=zangresult(maskt);
    zang_seg=zangresult(maskt);
    %     currtf=tf(maskt);
    
    [dmap]=build_dmap(mtc_seg,zang_seg,msz,sigma);
    dmap=dmap/max(dmap(:));
    %     [coords]=pick0spot(dmap);1
%     dmap(148-5:148+5,112-5:112+5)=0;  
%     if ii==15
%     dmap(203:250,87:109)=0;  
%     end
    dipshow(h,dmap,'lin');
%     figure(h)
%     imagesc(dmap)
%     colormap(gray)
%     axis equal;
%     axis off
    hold on
%     keyboard
    [mephi_ini,cpeak]=find_mephi_ini(dmap,srange,msz);
%     stopval=max(stopval,cpeak./4);
    
%     keyboard
try
    if stopval0==0
        [mephi,uwmephi,mephipp]=find_mephi(dmap,mephi_ini,maxangle,srange,0);
        I=mephi(:,3);
        N=length(I);
        xx=(1:N)';
        inis=[max(I), length(I)/2, 5];
        [c]=fit1dgaussian(I,xx,inis);
        figure(2);plot(xx,I,'b-'); hold on
        Fy = c(1).* exp(-(xx-c(2)).^2./2./c(3)./c(3))+c(4);
        plot(xx,Fy,'r-');hold off
        if c(2)<N/2
            stopval=max(I(min(round(c(2)+2*c(3)),N)),Fy(min(round(c(2)+2*c(3)),N)));
        else
            stopval=max(I(max(round(c(2)-2*c(3)),1)),Fy(max(round(c(2)-2*c(3)),1)));
        end
        stopval
        stopval=min(stopval,0.50);
        stopval=max(stopval,0.20);       
    end
    [mephi,uwmephi,mephipp]=find_mephi(dmap,mephi_ini,maxangle,srange,stopval);
catch
    keyboard
end
    [zest_f,zerr_f,mephimask]=mephi_zest(mephipp,uwmephi,mtc_seg,zang_seg,lambda,Nmedia,centermtc,freq);

    % scatter(uwmephi(:,1),uwmephi(:,2))
    zest=cat(1,zest,zest_f(:));
    zerr=cat(1,zerr,zerr_f(:));
    zmask=cat(1,zmask,mephimask(:));
    
    mep1=(mephi(:,1)+pi).*msz./2./pi;
    mep2=(mephi(:,2)+pi).*msz./2./pi;
    figure(1);scatter(mep1,mep2,'r')
%     keyboard
% %     size(zest)
    % plot
%     figure
%     scatter(mtc_seg,zang_seg,9)
%     hold on
%     scatter(mephi(:,1),mephi(:,2),'r')
%     plot(-pi:0.1:pi,ppval(mephipp,-pi:0.1:pi),'r--','linewidth',6)
%     xlim([min(mephi(:,1)) max(mephi(:,1))]);
%     xlim([-pi pi]);
%     ylim([-2*pi 2*pi]);
    hold off
    pause(eps)
    %% check stopVal 
%     if ii ==1& continueAll ==0;
%         choice = questdlg('Change Stop Value?','','Continue','Change','Use Same Value for All','Exit','Exit')
%         switch choice
%             case 'Continue'
%             case 'Change'
%                 options.WindowStyle = 'normal';
%                 x = inputdlg('Enter Stop Value for Dmap','Stopval',1,{'0.2'},options);
%                 if isempty(x)
%                     return
%                 end
%                 stopval = str2num(x{1})
%                 ii = ii-1;
%             case 'Use Same Value for All'
%                 continueAll = 1;
%             case 'Exit'
%                 return
%         end
%     end

    ii = ii+1;
%     keyboard
%     
end
% close(h)
  