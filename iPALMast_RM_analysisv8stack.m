% for v14 iPALMast_analysisv14scmos.m

function [xresult yresult tresult zfresult zangresult llresult CRLBresult Iresult zast_err_result stacktot num_images zangctrresult bgresult subimstot sxtot sytot]=iPALMast_RM_analysisv8stack(foldername,nameroot,fmfile,astfile,anglefile,llrthresh,scmos_cali_file,det_thresh,CRLB_thresh,sigma_thresh,centers,handles)
% v7
% Thresh is an input from outside of the function

% v6
% Now loads sCMOS calibration data

% v5
% Now loads and analyze sCMOS data format (dcimg) file

% v4
% now code use q1+q2+q3+q4 as the sum channel

% v3
% code now output zf and z_ang for estimation of airphase outside this
% analysis code

% v2
% autophi

% close all
% addpath('iPALMast_calibration\');

%% parameters
%
warning('off','all')

tmpld=load(scmos_cali_file);
offsetim=tmpld.offsetim; % for scmos, yes for now
varim=tmpld.varim;
gainim=tmpld.g3;    % for current calibration file, the whole quadrant is not filed with valid number such that I have use the average gain to replace the real gain.
fit_flag=str2double(get(handles.fit_flag,'string')); % 1: 4pi, 2: 2D, 3: 3D
if fit_flag==2
    subsz=7;   % sub region size
else
    subsz=13;
end
% offsetim=ones(164,1108)*100;
% varim=ones(164,1108)*13;
% gainim=ones(164,1108)*2.2;

caliims=cat(3,offsetim,varim,gainim);
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
%% load one quadrant of the stack
% foldername=['E:\DATA_iPALM\Biological_Sample\Feb20th_2014\MT1\'];
% nameroot=['MT1_642_000_000.tif'];
files=dir([foldername nameroot]);
lastf=0;
subimstot=[];

for ff=1:1:numel(files)
    close all
    ff
    filename=files(ff).name;
%     centers=[279 454 1581 1756];
    display(['Reading dcimg files...']);
    [~,qds,calicrops]=iPALM_readdcimg([foldername filename],centers,caliims);
%     calicrops(:,:,3,:)=calicrops(:,:,3,:).*0+mean(mean(mean(calicrops(:,:,3,:))));      % gain image is currently constant
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
%     stacknum=0;
    num_images=size(qds,3);
    clear mex
    qds=[];
    
    %% median filtering
    if handles.MedFilter.Value
        tic
        [qd1,qd2,qd3,qd4]=medFilter(qd1,qd2,qd3,qd4,num_images);
        toc
    end

    %% find affine transformation
% if ff==1
%     tic
%     fmfile=GetTransform(qd1,qd2,qd3,qd4,foldername,filename);
%     toc
% end

    %% find affine transformation
%     if ff==1
%         qdall=cat(4,qd1,qd2,qd3,qd4);
%         qdall(qdall<=0)=1e-37;
%         se1=strel('disk',5);
%         for ii=1:4
%             for j=1:num_images
%                 ROI=qdall(:,:,j,ii);
% %                 [W2 W3] = waveletTransform(ROI,1,3);
%                 TB=imreconstruct(imerode(ROI,se1),ROI); % rolling ball filter
%                 ROI1=imsubtract(ROI,TB);
% %                 ROI1=W2;
%                 qdall(:,:,j,ii)=ROI1;
%             end
%         end
%         zm_all=[];
%         trans_all=[];
%         ang_all=[];
%         R=[];
%         invR=[];
%         for ii=1:1:4
%             im1=mean(qdall(:,:,:,1),3);
%             im2=mean(qdall(:,:,:,ii),3);
%             [zm,trans,ang] = fmmatch(im2,im1);
%             [out,R(:,:,ii)] = find_affine_trans(im2, im1, [[zm zm],trans,ang]);
%             zm_fin=out(1:2);
%             trans_fin=out(3:4);
%             ang_fin=out(5);
%             
%             [imout]=affine_trans(im2,zm_fin,trans_fin,ang_fin);
%             zm_all(ii,:)=zm_fin;
%             trans_all(ii,:)=trans_fin;
%             ang_all(ii,:)=ang_fin;
% %             joinchannels('RBG',im1,imout);
%             
%             [zm,trans,ang] = fmmatch(im1,im2);
%             [~,invR(:,:,ii)] = find_affine_trans(im1, im2, [[zm zm],trans,ang]);
%         end
%         filestr=[foldername filename(1:end-6) '_FMTtransform_' datestr(now,'yyyymmdd'),'.mat'];
%         save(filestr,'zm_all','trans_all','ang_all','R','invR');
%         fmfile=filestr;
%     end
    
%     %% split channel and flip left to right
%     flip_flag=[0 0 1 1];
%     padflag=0;
%     padsz=512;
%     [qd1 qd2 qd3 qd4]=iPALMast_chSplit(ims,flip_flag,padflag,padsz);
%     
    %% rotate and align
    %     addpath 'iPALMast_calibration'    
    tic
    [q1 q2 q3 q4]=iPALMast_RotAlign_FMT(qd1,qd2,qd3,qd4,fmfile);
    [v1 v2 v3 v4]=iPALMast_RotAlign_FMT(varim(:,:,1),varim(:,:,2),varim(:,:,3),varim(:,:,4),fmfile);
    [g1 g2 g3 g4]=iPALMast_RotAlign_FMT(gainim(:,:,1),gainim(:,:,2),gainim(:,:,3),gainim(:,:,4),fmfile);
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
    
    %% test
%     [shift1,shift2]=drift_correction_core(double(mean(q1,3)),double(mean(q2,3)),[0 0]);
%     [shift1,shift2]=drift_correction_core(double(mean(q1,3)),double(mean(q3,3)),[0 0]);
%     [shift1,shift2]=drift_correction_core(double(mean(q1,3)),double(mean(q4,3)),[0 0]);
    
    %% find local centers
    sumim1=q1+q2+q3+q4;
    q1=[];q2=[];q3=[];q4=[];
    
    %% save to tif images
    if ff==1
        filestr=[foldername filename(1:end-5) 'tif'];
        if ~exist(filestr,'file')
%             tiffwrite(sumim1,filestr);
        end
%         filestr=[foldername filename(1:end-6) '_q1.tif'];
%         if ~exist(filestr,'file')
%             tiffwrite(q1,filestr);
%         end
%         filestr=[foldername filename(1:end-6) '_q2.tif'];
%         if ~exist(filestr,'file')
%             tiffwrite(q2,filestr);
%         end
%         filestr=[foldername filename(1:end-6) '_q3.tif'];
%         if ~exist(filestr,'file')
%             tiffwrite(q3,filestr);
%         end
%         filestr=[foldername filename(1:end-6) '_q4.tif'];
%         if ~exist(filestr,'file')
%             tiffwrite(q4,filestr);
%         end
    end
    
    %%
    sumv=v1+v2+v3+v4;
    sumvg=v1./g1./g1+v2./g2./g2+v3./g3./g3+v4./g4./g4;
    sumim1(sumim1<=1e-6)=1e-6;
%     sumg=g1+g2+g3+g4;
%     thresh=thresh;
%     subsz=13;   % sub region size
    pick_flag=1; % use range
    tic
    [sub_regions tlz locmaxc subvar_g]=iPALMast_sumim_seg(sumim1,det_thresh,subsz,sumv,sumvg,pick_flag,fit_flag);
    toc
%     clear qd1 qd2 qd3 qd4
%     sum13=q1+q3;
%     sum24=q2+q4;
%     [sub_regions13]=iPALMast_sumim_seg(sum13,thresh,subsz,pick_flag,pick_sz);
%     [sub_regions24]=iPALMast_sumim_seg(sum24,thresh,subsz,pick_flag,pick_sz);
    %% fit sigmax sigmay
%     display(['dsdas' num2str(2) 'jhaa']);
    display(['A total of ' num2str(size(sub_regions,3)) ' subregions were detected. Start sCMOS_sigmaxy fitting']);
    %% show detection result
    f=1000;
    id=locmaxc(:,3)==f-1;
    V=locmaxc(id,:)+1;
    figure;imshow(sumim1(:,:,f),[]);hold on;plot(V(:,1),V(:,2),'bo');pause(1)
        
    %%
    sumim1=[];
    tic
    Nt=size(sub_regions,3);
    Nf=ceil(Nt/10000);
    P=[];
    CRLB=[];
    LL=[];
    for k=1:Nf
        st=(k-1)*10000+1;
        et=min(k*10000,Nt);
        P1=[];
        CRLB1=[];
        LL1=[];
%         [P1 CRLB1 LL1]=sCMOS_MLE_CUDA42(single(sub_regions(:,:,st:et)),1.4,200,4,single(subvar_g(:,:,st:et)),single(subvar_g(:,:,st:et).*0+1));
        if fit_flag==2
            [P1 CRLB1 LL1]=sCMOS_MLE_DB_YL_CUDA42(single(sub_regions(:,:,st:et)),1.4,100,6,single(subvar_g(:,:,st:et)),single(subvar_g(:,:,st:et).*0+1));
            P1(:,6)=P1(:,5);
        else
            [P1 CRLB1 LL1]=sCMOS_MLE_DB_YL_CUDA42(single(sub_regions(:,:,st:et)),1.4,200,7,single(subvar_g(:,:,st:et)),single(subvar_g(:,:,st:et).*0+1));
        end         
        P=cat(1,P,P1);
        CRLB=cat(1,CRLB,CRLB1);
        LL=cat(1,LL,LL1);
    end
    toc
    display(['Fitting finished...Start Z_ast estimation']);
    
    sub_regions=[];subvar_g=[];  
%     xstep=(1:1:length(P(:,1)))';
    xco=P(:,2);
    yco=P(:,1);
    I=P(:,3);
    llr=-2*LL;
    sigmax=P(:,5);
    sigmay=P(:,6);
    bg=P(:,4);
    % hist(llr(llr<1000),200);
    % test code
%     mask=llr<600;
%     sigx=sigmax(mask);
%     sigy=sigmay(mask);
%     scatter(sigx,sigy);
    %
    %% filter
    rr=(subsz-1)/2;
    maskxy=abs(xco-rr)<3&abs(yco-rr)<3;
%     maskxy=xco<subsz-2&xco>2&yco<subsz-2&yco>2;
%     masks=sigmax<subsz/2&sigmay<subsz/2;
%     masks=sigmax<sigma_thresh&sigmay<sigma_thresh;
%     maskll=llr<llrthresh; %<-- change
    maskothers=CRLB(:,1)>0&CRLB(:,2)>0&CRLB(:,1)<1&CRLB(:,2)<1;
    
    mask=maskxy&maskothers;%&masks&maskll&maskothers;
%     mask=maskxy&masks&maskll;
    
    xf=xco(mask);
    yf=yco(mask);
    llr_f=llr(mask);
    CRLB(:,1:2)=sqrt(CRLB(:,1:2));
    CRLB_f=CRLB(mask,:);
    I_f=I(mask);
    tlz_f=tlz(mask,:);
    %     xstepf=xstep(mask);
    sigmaxf=sigmax(mask);
    sigmayf=sigmay(mask);
    %     tlz_f=tlz(mask,:);
    locmaxc_f=locmaxc(mask,:);
    stacknum_f=repmat(stacknum(1),size(xf));
    bg_f=bg(mask);
%     subregion_sel=sub_regions(:,:,mask);
    
    %% z ast initial guess
    ast_config=astfile;
    zf=[];
    zerr=[];
%     addpath 'calibration_files\'
%     addpath 'Z:\Fang\MATLAB\Projects\iPALM\_dev_branch3_iPALMscmos_slope\calibration_files';
    
    %% z estimation by astigmatism by look up table
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
    
    if fit_flag==2
        zerr=zeros(length(xf),1);
        zf=zeros(length(xf),1);
    else
        tic
        N=length(sigmaxf);
        zf=zeros(N,1);
        zerr=zeros(N,1);
        parfor ii=1:N
            w0=[sigmaxf(ii),sigmayf(ii)];
            dist=Yf.^0.5-(ones(M,1)*w0).^0.5;
            dist=sum(dist.^2,2);
            id=dist==min(dist);
            zerr(ii,1)=min(dist);
            zf(ii,1)=mean(zdata(id));
        end
        toc
    end
    
%     tic
%     parfor ii=1:1:length(sigmaxf)
%         [zf(ii) zerr(ii)]=iPALMast_iniz_astfit(sigmaxf(ii),sigmayf(ii),ast_config);
%     end
%     toc
    
    display(['Z_ast fitting finished...Start phase estimation']);
    %% determine the phase
    
    if isempty(locmaxc_f)
        continue
    end
    
    if fit_flag==1
        tic
        [z_ang,~,~,ang_ctr subims z_angd]=iPALMast_RM_zangv3(anglefile,qd1,qd2,qd3,qd4,locmaxc_f,subsz,xf,yf,offsetim,varim,gainim,fmfile,tlz_f,bg_f);
        toc
    else
        z_ang=zeros(length(locmaxc_f),1);
        ang_ctr=ones(length(locmaxc_f),1);
    end

    qd1=[];qd2=[];qd3=[];qd4=[];
%     [z_ang]=iPALMast_RM_zang(q1,q2,q3,q4,locmaxc_f,subsz,xf,yf);
    display(['Phase estimation finished...Start the next file...']);
    %% collect data
    xest=xf+tlz_f(:,2);
    yest=yf+tlz_f(:,1);
    
    sxtot=cat(1,sxtot,sigmaxf(:));
    sytot=cat(1,sytot,sigmayf(:));
%     subimstot=cat(3,subimstot,subims);
    xresult=cat(1,xresult,xest);
    yresult=cat(1,yresult,yest);
    bgresult=cat(1,bgresult,bg_f);
    zfresult=cat(1,zfresult,zf(:));
    zangresult=cat(1,zangresult,z_ang(:));
%     zangoresult=cat(1,zangoresult,z_ang_o);
    zangctrresult=cat(1,zangctrresult,ang_ctr(:));
    tresult=cat(1,tresult,tlz_f(:,3)+(ff-1)*num_images);
    zast_err_result=cat(1,zast_err_result,zerr(:));
    llresult=cat(1,llresult,llr_f);
    CRLBresult=cat(1,CRLBresult,CRLB_f);
    Iresult=cat(1,Iresult,I_f);
    lastf=lastf+tlz(end,3);
    stacktot=cat(1,stacktot,stacknum_f); 
    
%     mtc=(sigmaxf.^3./sigmayf-sigmayf.^3./sigmaxf)./40*2*pi;
%     mask=abs(mtc)<pi & ang_ctr(:)>0;
%     dmap=build_dmap(mtc(mask),z_ang(mask),256,3);
%     dipshow(dmap,'lin');  
%     pause(1);
end

%% filtering 
mask=zangctrresult>0&zangctrresult<sqrt(2)&Iresult>0&llresult>0;
sxtot=sxtot(mask);
sytot=sytot(mask);
xresult=xresult(mask);
yresult=yresult(mask);
bgresult=bgresult(mask);
zfresult=zfresult(mask);
zangresult=zangresult(mask);
zangctrresult=zangctrresult(mask);
tresult=tresult(mask);
zast_err_result=zast_err_result(mask);
llresult=llresult(mask);
CRLBresult=CRLBresult(mask,:);
Iresult=Iresult(mask);
stacktot=stacktot(mask);

%% old code backup

%     mask=zerr<prctile(zerr,30);
%     zerrf=zerr(mask);
%     zff=zf(mask);
%     zangf=z_ang(mask);
%
%     %% determine the airphase
%     [phiest meanest]=iPALM_autophi(zff(:),zerrf(:),zangf,anglefile,1);
%     %     if ff==1
%     %         angfile=anglefile;
%     %         tmp=load(angfile);
%     %         scatter(zff,zangf,'*r');
%     %         dz=input('pick z (nm):','s');
%     %         dz=str2double(dz);
%     %         dphi=input('pick the corresponding phi (rad):','s');
%     %         dphi=str2double(dphi);
%     %
%     %         phi_air=wrapToPi(dphi-tmp.ang_coeff(1)*dz);
%     %     end
%     phi_air=phiest;
%     %% get z estimate
%     tmp=load(anglefile);
%     [zest z_error]=iPALMast_angle_z(z_ang,zf,phi_air,tmp.ang_coeff(1));
%
%     xest=xf+tlz_f(:,2);
%     yest=yf+tlz_f(:,1);
%
%     xresult=cat(1,xresult,xest);
%     yresult=cat(1,yresult,yest);
%     zresult=cat(1,zresult,zest);
%     tresult=cat(1,tresult,tlz_f(:,3)+lastf);
%     zerrresult=cat(1,zerrresult,z_error(:));
%     llresult=cat(1,llresult,llr_f);
%     Iresult=cat(1,Iresult,I_f);
%     zast_err_result=cat(1,zast_err_result,zerr(:));
%     phi_air_result=cat(1,phi_air_result,phi_air(:));
%     lastf=lastf+tlz(end,3);