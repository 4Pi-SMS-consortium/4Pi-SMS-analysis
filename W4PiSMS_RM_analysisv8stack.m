function [xresult yresult tresult zfresult zangresult llresult CRLBresult Iresult zast_err_result stacktot num_images imagesz zangctrresult bgresult subimstot sxtot sytot]=W4PiSMS_RM_analysisv8stack(foldername,nameroot,fmfile,astfile,anglefile,llrthresh,scmos_cali_file,det_thresh,CRLB_thresh,sigma_thresh,centers,handles)
%% parameters
fit_flag=str2double(get(handles.fit_flag,'string')); % 1: 4Pi, 2: 2D, 3: 3D
if fit_flag==2
    subsz=7;   % sub region size
else
    subsz=13;
end

%% camera calibration file
load(handles.scmos_cali_file);
caliims=cat(3,offsetim,varim,gainim);

%% 
xresult=[];
yresult=[];
zfresult=[];
zangresult=[];
tresult=[];
zast_err_result=[];
llresult=[];
CRLBresult=[];
Iresult=[];
stacktot=[];
zangctrresult=[];
bgresult=[];
sxtot=[];
sytot=[];

%% 
files=dir([foldername nameroot]);
lastf=0;
subimstot=[];
for ff=1:1:numel(files)
    disp(['Reading dcimg file: ',num2str(ff)]);
    filename=files(ff).name;
    tic
    [~,qds,calicrops]=W4PiSMS_readdcimg([foldername filename],centers,caliims);
    toc
    tmpim=calicrops(:,:,3,:);
    tmpim(tmpim<1.3|tmpim>3.5)=mean(mean(mean(calicrops(:,:,3,:))));
    calicrops(:,:,3,:)=tmpim;
    offsetim=squeeze(calicrops(:,:,1,:));
    varim=squeeze(calicrops(:,:,2,:));
    gainim=squeeze(calicrops(:,:,3,:));
    stacknum=getstepnum(filename);
    qd1=(qds(:,:,:,1)-repmat(offsetim(:,:,1),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,1),[1 1 size(qds,3) 1]);
    qd2=(qds(:,:,:,2)-repmat(offsetim(:,:,2),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,2),[1 1 size(qds,3) 1]);
    qd3=(qds(:,:,:,3)-repmat(offsetim(:,:,3),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,3),[1 1 size(qds,3) 1]);
    qd4=(qds(:,:,:,4)-repmat(offsetim(:,:,4),[1 1 size(qds,3) 1]))./repmat(gainim(:,:,4),[1 1 size(qds,3) 1]);
    num_images=size(qds,3);
    imagesz=size(qds,1);
    qds=[];
        
    %% rotate and align 
    tic
    [q1 q2 q3 q4]=W4PiSMS_RotAlign_FMT(qd1,qd2,qd3,qd4,fmfile);
    [v1 v2 v3 v4]=W4PiSMS_RotAlign_FMT(varim(:,:,1),varim(:,:,2),varim(:,:,3),varim(:,:,4),fmfile);
    [g1 g2 g3 g4]=W4PiSMS_RotAlign_FMT(gainim(:,:,1),gainim(:,:,2),gainim(:,:,3),gainim(:,:,4),fmfile);
    toc
    
    maskg=g1<=1|g2<=1|g3<=1|g4<=1;
    v1(maskg)=1e7;
    v2(maskg)=1e7;
    v3(maskg)=1e7;
    v4(maskg)=1e7;
    g1(maskg)=1e-7;
    g2(maskg)=1e-7;
    g3(maskg)=1e-7;
    g4(maskg)=1e-7;

    sumim1=q1+q2+q3+q4;
    sumim1(sumim1<=1e-6)=1e-6;
    q1=[];q2=[];q3=[];q4=[];
    
    %% save to tif images
    if ff==1
        filestr=[foldername filename(1:end-5) 'tif'];
        if ~exist(filestr,'file')
            tiffwrite(sumim1,filestr);
        end
    end
    
    %%
    sumv=v1+v2+v3+v4;
    sumvg=v1./g1./g1+v2./g2./g2+v3./g3./g3+v4./g4./g4;   
    pick_flag=1; % use range
    tic
    [sub_regions tlz locmaxc subvar_g]=W4PiSMS_sumim_seg(sumim1,det_thresh,subsz,sumv,sumvg,pick_flag,fit_flag);
    toc
    disp(['A total of ' num2str(size(sub_regions,3)) ' subregions were detected. Start sCMOS_sigmaxy fitting']);
    
    %% show detection result
    f=100;
    id=locmaxc(:,3)==f-1;
    V=locmaxc(id,:)+1;
    close all
    figure;imshow(sumim1(:,:,f),[0 500],'InitialMagnification',250);hold on;plot(V(:,1),V(:,2),'bo');pause(eps);
        
    %%
    sumim1=[];
    tic
    Nt=size(sub_regions,3);
    Nf=ceil(Nt/100000);
    P=[];
    CRLB=[];
    LL=[];
    for k=1:Nf
        st=(k-1)*100000+1;
        et=min(k*100000,Nt);
        P1=[];
        CRLB1=[];
        LL1=[];
        if fit_flag==2
            [P1,CRLB1,LL1]=mleFit_LM(single(sub_regions(:,:,st:et)),2,50,1.4,0,0,0);
            P1(:,6)=P1(:,5);
        else
            [P1,CRLB1,LL1]=mleFit_LM(single(sub_regions(:,:,st:et)),4,50,1.4,0,0,0);
        end         
        P=cat(1,P,P1);
        CRLB=cat(1,CRLB,CRLB1);
        LL=cat(1,LL,LL1);
    end
    toc
    disp('Fitting finished...Start Z_ast estimation');
    
    sub_regions=[];subvar_g=[];  
    xco=P(:,2);
    yco=P(:,1);
    I=P(:,3);
    llr=-2*LL;
    sigmax=P(:,5);
    sigmay=P(:,6);
    bg=P(:,4);

    %% filter
    rr=(subsz-1)/2;
    maskxy=abs(xco-rr)<3&abs(yco-rr)<3;
    maskothers=CRLB(:,1)>0&CRLB(:,2)>0&CRLB(:,1)<1&CRLB(:,2)<1;
    mask=maskxy&maskothers;%&masks&maskll&maskothers;
    
    xf=xco(mask);
    yf=yco(mask);
    llr_f=llr(mask);
    CRLB(:,1:2)=sqrt(CRLB(:,1:2));
    CRLB_f=CRLB(mask,:);
    I_f=I(mask);
    tlz_f=tlz(mask,:);
    sigmaxf=sigmax(mask);
    sigmayf=sigmay(mask);
    locmaxc_f=locmaxc(mask,:);
    stacknum_f=repmat(stacknum(1),size(xf));
    bg_f=bg(mask);
    
    %% z ast initial guess
    ast_config=astfile;
    zf=[];
    zerr=[];
    
    %% z estimation by astigmatism by look up table
    if fit_flag==3  % astigmatism
        if ff==1
            tmp=load(ast_config);
            zdata=(0:1200)';
            Yf=[];
            params=tmp.estx;
            w=params(1);
            c=params(2);
            d=params(3);
            A=params(4);
            B=params(5);
            Yf(:,1) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
            params=tmp.esty;
            w=params(1);
            c=params(2);
            d=params(3);
            A=params(4);
            B=params(5);
            Yf(:,2) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
            dif=abs(Yf(:,1)-Yf(:,2));
            id=find(dif==min(dif(:)));

            % center
            zdata=id-800:id+800;
            Yf=[];
            params=tmp.estx;
            w=params(1);
            c=params(2);
            d=params(3);
            A=params(4);
            B=params(5);
            Yf(:,1) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
            params=tmp.esty;
            w=params(1);
            c=params(2);
            d=params(3);
            A=params(4);
            B=params(5);
            Yf(:,2) = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
            M=length(zdata);
        end
        
        N=length(sigmaxf);
        zf=zeros(N,1);
        zerr=zeros(N,1);
        tic
        parfor ii=1:N
            w0=[sigmaxf(ii),sigmayf(ii)];
            dist=Yf.^0.5-(ones(M,1)*w0).^0.5;
            dist=sum(dist.^2,2);
            id=dist==min(dist);
            zerr(ii,1)=min(dist);
            zf(ii,1)=mean(zdata(id));
        end
        toc
    else
        zerr=zeros(length(xf),1);
        zf=zeros(length(xf),1);
    end
    disp('Z_ast fitting finished...Start phase estimation');
    
    %% determine the phase    
    if isempty(locmaxc_f)
        continue
    end
    
    if fit_flag==1
        tic
        [z_ang,ang_ctr]=W4PiSMS_RM_zangv3(anglefile,qd1,qd2,qd3,qd4,locmaxc_f,subsz,xf,yf,offsetim,varim,gainim,fmfile,tlz_f,bg_f);
        toc
    else
        z_ang=zeros(length(locmaxc_f),1);
        ang_ctr=ones(length(locmaxc_f),1);
    end

    qd1=[];qd2=[];qd3=[];qd4=[];
    disp('Phase estimation finished...Start next file...');
    
    %% collect data
    xest=xf+tlz_f(:,2);
    yest=yf+tlz_f(:,1);    
    sxtot=cat(1,sxtot,sigmaxf(:));
    sytot=cat(1,sytot,sigmayf(:));
    xresult=cat(1,xresult,xest);
    yresult=cat(1,yresult,yest);
    bgresult=cat(1,bgresult,bg_f);
    zfresult=cat(1,zfresult,zf(:));
    zangresult=cat(1,zangresult,z_ang(:));
    zangctrresult=cat(1,zangctrresult,ang_ctr(:));
    tresult=cat(1,tresult,tlz_f(:,3)+(ff-1)*num_images);
    zast_err_result=cat(1,zast_err_result,zerr(:));
    llresult=cat(1,llresult,llr_f);
    CRLBresult=cat(1,CRLBresult,CRLB_f);
    Iresult=cat(1,Iresult,I_f);
    lastf=lastf+tlz(end,3);
    stacktot=cat(1,stacktot,stacknum_f); 
end

%% filtering 
mask=zangctrresult>0&zangctrresult<sqrt(2)&Iresult>0&llresult>0;
sxtot=sxtot(mask);
sytot=sytot(mask);
xresult=xresult(mask);
yresult=yresult(mask);
bgresult=bgresult(mask);
llresult=llresult(mask);
CRLBresult=CRLBresult(mask,:);
Iresult=Iresult(mask);
tresult=tresult(mask);
zfresult=single(zfresult(mask));
zangresult=single(zangresult(mask));
zangctrresult=single(zangctrresult(mask));
zast_err_result=single(zast_err_result(mask));
stacktot=single(stacktot(mask));