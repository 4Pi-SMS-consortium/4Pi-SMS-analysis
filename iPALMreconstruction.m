function iPALMreconstruction(handles)

mainFolder = get(handles.pathMainfolder,'string');
index_selected = get(handles.pathSubfolder,'value');
filelist = get(handles.pathSubfolder,'string');
file_selected = filelist(index_selected);

photons=str2num(get(handles.photons,'string'));
llrthresh=str2num(get(handles.llrthresh2,'string'));
mephiBin = str2num(get(handles.mephiBin,'string'));
num_images = str2num(get(handles.numImages,'string'));
angctr=[0.4 1.2];
%     lambda=680;
%     Nmedia=1.33;
DC_cycle=1;
st_frame=3000*0;      %startframe
frmthresh=3000*10000;  %stopframe
stopval=str2double(get(handles.stopValue,'string'));     %stop value for Dmap
zthresh=[0 1200];
zasterr=[0 0.05];
phasehelper=0.6;%higher value makes the phase slope higher, recommend, 560: 1.0, 642: 0.6;
reverseflag=0; %642 ch does need to be reversed. from beginning to the end.
%%
% currentresultpath = get(handles.pathMainfolder,'string');
% files=dir([currentresultpath '*tmpresult*.mat']);

% global freqstr;
freqld=load('frequency_simulate4Pi_570_671.mat');

for ff=1:numel(file_selected)

    currentfile=file_selected{ff};
    I=find(currentfile=='\',1,'last');
    currentfolder=currentfile(1:I);
    fileName=currentfile(I+1:end);
    if ~isempty(strfind(fileName,'_642_'))
        channel = '642';
        reverseflag=0; %642 ch does need to be reversed. from beginning to the end.
        freqstr='f_690';       
        eval(['freq=freqld.' freqstr ';']);
    elseif ~isempty(strfind(fileName,'_561_'))
        channel = '561';
        reverseflag=0; %561 ch does not need to be reversed. from beginning to the end.
        freqstr='f_570';
         eval(['freq=freqld.' freqstr ';']);
    else
        channel = 'not 561 or 642'
    end
    savename=fileName(1:end-27);%name format foldername_642_tempresult_yymmdd.mat
    datestring = datestr(now,'yyyymmdd');
    set(handles.programStatus,'string',{'Current Analized File:',currentfile, ['Current Channel:' channel]})
    drawnow update
    %% focus center calibration
    %     ldast= load('Z:\Fang\MATLAB\Projects\iPALM\4PiSMS_sCMOS_astv19\calibration_files\bead9_642_Mar16th_2015_Astfit_031615.mat');
    %     ini=1000;
    %     [z_center]=findobj_center(ldast.estx,ldast.esty,ini);
    %     z_center
    %%
    %     tmpld=load([currentresultpath savename 'tmpresult']);
    tmpld=load(file_selected{ff});
    
    %% Drift correction by beads
    if handles.beadsDriftCorr.Value ==1
        Vin=[];
        Vin(:,1)=tmpld.xresult;
        Vin(:,2)=tmpld.yresult;
        Vin(:,3)=tmpld.tresult+1;
        Vin(:,5)=tmpld.Iresult;
        Vin(:,6)=tmpld.zangresult;
        Vin(:,7)=tmpld.zfresult;
        [Vout, drift, flag]=Beads_correction(Vin,reverseflag);
        if flag
            tmpld.xresult=Vout(:,1);
            tmpld.yresult=Vout(:,2);
            tmpld.zangresult=wrapToPi(Vout(:,3));
        end
    end
    
    %% reject beads
    Vin=[];
    Vin(:,1)=tmpld.xresult;
    Vin(:,2)=tmpld.yresult;
    Vin(:,3)=tmpld.tresult+1;
    maskbead=reject_beads(single(Vin),tmpld.num_images);
    
    %% rejection   
    mtc=(tmpld.sxtot.^3./tmpld.sytot-tmpld.sytot.^3./tmpld.sxtot)./40*2*pi; 
%     mtc=mtc/std(mtc)*pi;
%     mtc=(tmpld.sytot.^3./tmpld.sxtot-tmpld.sxtot.^3./tmpld.sytot)./40*2*pi;    
    % mtc mask
    maskmtc=mtc<pi&mtc>-pi &(tmpld.sxtot+tmpld.sytot)<6;% & tmpld.sxtot<2.5;
%     maskmtc=mtc<median(mtc)+1*std(mtc) & mtc>median(mtc)-1*std(mtc) & (tmpld.sxtot+tmpld.sytot)<6;
    % ll threshold
    maskll=tmpld.llresult<llrthresh;
    % crlb
    maskcrlb=tmpld.CRLBresult(:,1)<0.25 & tmpld.CRLBresult(:,2)<0.25;% & tmpld.CRLBresult(:,1)>0.045 & tmpld.CRLBresult(:,2)>0.045;
    % frame number threshold
    maskt=tmpld.tresult>=st_frame & tmpld.tresult<frmthresh;
%     maskt=ones(length(mtc),1);
%     id=tmpld.tresult>=204000&tmpld.tresult<206000;
%     maskt(id)=0;
    % contrast threshold
    maskzctr=tmpld.zangctrresult>angctr(1)&tmpld.zangctrresult<angctr(2);
    % background threshold
    maskbg=tmpld.bgresult<(median(tmpld.bgresult)+2*std(tmpld.bgresult))&tmpld.bgresult>(median(tmpld.bgresult)-2*std(tmpld.bgresult));
    % astigmatism fitting threshold
    maskzast_err=tmpld.zast_err_result<zasterr(2)&tmpld.zast_err_result>zasterr(1);
    % z range threshold
    maskzf=tmpld.zfresult<zthresh(2)&tmpld.zfresult>zthresh(1);
    % # photons range threshold
    maskI=tmpld.Iresult>photons&tmpld.Iresult<100000;
    % ROI mask
    maskxy=tmpld.xresult>5000/128&tmpld.xresult<16000/128&tmpld.yresult>4000/128&tmpld.yresult<15000/128;
    % mask ROI
%     str='E:\4Pi_two_color\2018-6-10\Cell12-0679-1072.roi';
%     roi=ReadImageJROI(str);
%     Rc=roi.mnCoordinates;
%     Rc(end+1,:)=Rc(1,:);
%     I=zeros(2150,2150);
%     BW=roipoly(I,Rc(:,1),Rc(:,2));
%     V=[tmpld.xresult*12.8,tmpld.yresult*12.8];
%     V=round(V);
%     L=length(V);
%     maskroi=zeros(L,1);
%     for jj=1:L
%         maskroi(jj)=BW(V(jj,2),V(jj,1));
%     end
    % apply
    mask=maskI&maskll&maskzctr&maskmtc&maskcrlb&maskt&maskbead;%&maskxy;%&maskroi;%&maskxy;%&maskroi;
    
    PR=(1:length(mask))';
    PR=PR(mask);
    rej=[];
    nl=length(mask);
    rej(1,1)=sum(maskmtc)/nl;
    rej(2,1)=sum(maskll)/nl;
    rej(3,1)=sum(maskcrlb)/nl;
    rej(4,1)=sum(maskzctr)/nl;
    rej(5,1)=sum(maskbg)/nl;
    rej(6,1)=sum(maskI)/nl;
    rej(7,1)=sum(maskxy)/nl;
    rej(8,1)=sum(mask)/nl;    
    
    %         subimstot2=tmpld.subimstot(:,:,mask,:);
    
    if num_images==0
        num_images=tmpld.num_images;
    end
    astfile=tmpld.astfile;
    anglefile=tmpld.anglefile;
    Iresult=tmpld.Iresult(mask);
    llresult=tmpld.llresult(mask);
    stacktot=tmpld.stacktot(mask);
    tresult=tmpld.tresult(mask)-st_frame;
    xresult=tmpld.xresult(mask);
    yresult=tmpld.yresult(mask);
    zangresult=tmpld.zangresult(mask);
    zast_err_result=tmpld.zast_err_result(mask);
    zfresult=tmpld.zfresult(mask);
    mtcresult=mtc(mask);
    sxresult=tmpld.sxtot(mask);
    syresult=tmpld.sytot(mask);
    zangctrresult=tmpld.zangctrresult(mask);
    crlbresult=tmpld.CRLBresult(mask,:);
    bgresult=tmpld.bgresult(mask,:);
       
    %%
    coords=[];
    coords(:,1)=xresult;
    coords(:,2)=yresult;
    szx=168;
    szy=168;
    im=cHistRecon(szx,szy,single(coords(:,2)),single(coords(:,1)),0);
    gaussim=gaussf(im,[1 1]);
    dipshow(gaussim);
    
    %% manual phase shift
%     [zangresult,mtcresult]=PhaseShift(zangresult,mtcresult,tresult,num_images);
     
    %% Phase drift correction
if handles.phaseDriftCorr.Value ==1
    str=[currentfolder savename '_phase.mat'];
    if exist(str,'file')
        load(str);
    else
        [zangresult,mtcresult]=PhaseCorrection(zangresult,mtcresult,tresult,stacktot,num_images,reverseflag);
        save(str,'zangresult','mtcresult');
    end
end

    %% dmap test
    fnum=num_images;
    segnum=ceil((max(tresult))/fnum);
    zest=[];
    zerr=[];
    zmask=[];
    Dmap0=[];
    Dmap=[];
    for ii=1:1:segnum
        close all
        st=(ii-1)*fnum;
        if ii==segnum
            ed=max(tresult);
        else
            ed=(ii)*fnum-1;
        end
        maskt=tresult>=st&tresult<=ed;
        mtc_seg=mtcresult(maskt);
        zang_seg=zangresult(maskt);
        [dmap]=build_dmap(mtc_seg,zang_seg,256,5);
        Dmap(:,:,ii)=dmap/max(dmap(:));
    end
    
    unistack=unique(stacktot);
    N=numel(unistack);
    if N>1
        OD=[];
        for j=1:N
            od=j:N:segnum;
            OD=cat(1,OD,od');
        end
        for ii=1:1:segnum
            Dmap0(:,:,ii)=Dmap(:,:,OD(ii));
        end
        Dmap=Dmap0;
    end
    
    if handles.phaseDriftCorr.Value ==1
        str=([currentfolder savename '_dmap_DC.tif']);
    else
        str=([currentfolder savename '_dmap.tif']);
    end
    tiffwrite(Dmap*10000,str);
    programStatus_string = get(handles.programStatus,'string');
    programStatus_string = [programStatus_string;['Save Dmap Image to:' str]];
    set(handles.programStatus,'string',programStatus_string)
%     continue
    
    %% phase estimate for every optical section separately
    unistack=unique(stacktot);
    phase_cyc=DC_cycle;
    xout=[];
    yout=[];
    zout=[];
    shifts=[];
    z_err_out=[];
    zf_out=[];
    t_out=[];
    ll_out=[];
    I_out=[];
    bg_out=[];
    pr_out=[];
    Z_ast_err_out=[];
    crlb_out=[];
    zangctrl_out=[];
    phaseslope=[];
    close all
    for ss=1:numel(unistack)
        close all
        currst=unistack(ss);
        maskst=(stacktot==currst);
        currt=tresult(maskst);
        currI=Iresult(maskst);
        currll=llresult(maskst);
        currx=xresult(maskst);
        curry=yresult(maskst);
        currzang=zangresult(maskst);
        currzangctrl=zangctrresult(maskst);
        currcrlb=sqrt((crlbresult(maskst,1).^2+crlbresult(maskst,2).^2)/2);
        currzasterr=zast_err_result(maskst);
        currzfresult=zfresult(maskst);
        currmtc=mtcresult(maskst);
        currbg=bgresult(maskst);
        currpr=PR(maskst);
        % get phase per cycle
        
        fnum=phase_cyc.*num_images;
%         fnum=phase_cyc.*numel(unistack).*num_images;
        
        if numel(unistack)>1
            currt=currt-(ss-1)*tmpld.num_images;
            fnum1=tmpld.num_images;
            cycle=floor(currt/(fnum1*numel(unistack)));
            remn=rem(currt,fnum1*numel(unistack));
            currt=cycle*fnum1+remn;
        end
        % %     fnum=phase_cyc.*tmpld.num_images;
        %     plotflag=1;
        %     airsmflag=0;
        %     [phi_ang_est]=iPALMast_phase_estimatev2(anglefile,currzfresult,...
        %         double(currzang),currt,currll,currzasterr,fnum,...
        %         plotflag,[],phasehelper);
        %     airphases=phi_ang_est(1:end-1);
        %     phaseslope(ss)=phi_ang_est(end);
        
        %     airphases=unwrap(airphases);
        %     alpha=0.05;
        %     [airp1 delind]=del_outlier(airphases,alpha);
        %     [airp2]=fillval(airp1);
        %     airphases=airp2;
        %     % plot
        %     % plot(airphases,'*');
        %     % hold on
        %     % plot(double(smooth(airphases)),'*r');
        %     if airsmflag==1
        %         airphases=double(smooth(airphases));
        %     end
        
        %     airphases=wrapToPi(airphases);
        
        %     [currzresult z_err]=iPALMast_zsegs2(currzang,currzfresult,currt,airphases,phaseslope(ss),fnum);
        %
        display(['processing segment:' num2str(ss)]);
        
        centermtc=[];
%         stopval=2;
%         SV=[0.2,0.2,0.2,0.3,0.4,0.4,0.5,0.6,0.45,0.6];
%         SV=[0.1,0.1,0.2,0.2,0.2];
%         SV=[0.1,0.1,0.15,0.1,0.1,0.1,0.15,0.15,0.2,0.25,0.25,0.2,0.15,0.2,0.1,0.1,0.1,0.1,0.1,0.1,0.01];
%         stopval=SV(ss);
%         if ss>=23
%             stopval1=0.05;
%         else
            stopval1=stopval;
%         end
        [currzresult z_err mephimask]=Mephi_z_4PiSMS(currzang,currmtc,currt,fnum,mephiBin,stopval1,centermtc,freq);
        %         close all
%         pause(0.5)
        %     figure
        %     hist(z_err,20)
                       
        % apply mask to every data

%         maskall=ones(length(z_err),1)>0;
%         currzresult=currzfresult(maskall>0);

        maskzerr=z_err<60;        
        maskall=maskzerr&(mephimask>0)&abs(currzresult)<700;
        currzresult=currzresult(maskall>0);        
        currx=currx(maskall>0);
        curry=curry(maskall>0);
        currt=currt(maskall>0);
        
%         if numel(unistack)>1
%             cycle=floor(currt/(fnum*numel(unistack)));
%             remn=rem(currt,fnum*numel(unistack));
%             currt=cycle*fnum+remn;
%         end
    
        z_err=z_err(maskall>0);
        currzfresult=currzfresult(maskall>0);
        currll=currll(maskall>0);
        currI=currI(maskall>0);
        currzangctrl=currzangctrl(maskall>0);
        currcrlb=currcrlb(maskall>0);
        currzasterr=currzasterr(maskall>0);
        currbg=currbg(maskall>0);
        currpr=currpr(maskall>0);
        
        driftstr=[currentfolder savename '_drift.mat'];
        interpflag=handles.Interpolation.Value;
        [xout{ss},yout{ss},zout{ss},shifts{ss}]=iPALM_driftcorrection_RedunLSv6(currx.*128,curry.*128,currzresult,currt,num_images,reverseflag,interpflag,driftstr);
%         [xout{ss},yout{ss},zout{ss},shifts{ss}]=iPALM_driftcorrection_RedunLSv11(currx.*128,curry.*128,currzresult,currt,num_images,reverseflag,interpflag);
%         [zout{ss},shifts_z]=iPALM_driftcorrection_RedunLSv9(currx.*128,currzresult,currt,num_images,0);
%         [xout{ss},yout{ss},zout{ss},shifts{ss}]=iPALM_driftcorrection_RedunLSv10(currx.*128,curry.*128,currzresult,currt,num_images,interpflag,driftstr);
               
        %         [flag]=iPALM2vutarav2(currentresultpath,[savename 'v13'],1,{xout},{yout},{zout},{round(tout./10)},{3*llrout});
        %         [flag]=iPALM2vutarav2(currentresultpath,[savename 'v13_driftcorrect'],1,{xout2},{yout2},{zout2},{round(tout./10)},{3*llrout});
        %
        %
%         [xout{ss} yout{ss} zout{ss} shifts{ss}]=iPALM_driftcorrection_RedunLSv7_3d(currx.*128,curry.*128,currzresult,currt,tmpld.num_images,reverseflag);
        
        z_err_out{ss}=z_err;
        zf_out{ss}=currzfresult;
        ll_out{ss}=currll;
        t_out{ss}=currt;
        I_out{ss}=currI;
        zangctrl_out{ss}=currzangctrl;
        crlb_out{ss}=currcrlb;
        Z_ast_err_out{ss}=currzasterr;
        bg_out{ss}=currbg;
        pr_out{ss}=currpr;
    end
    %      time = clock;
    %     stringTime = ['_' num2str(time(1)) '_' num2str(time(2)) '_' num2str(time(3)) '_' num2str(time(4)) '_' num2str(time(5)) '_' num2str(time(6)) '.mat'];
    %     save([currentresultpath savename 'DCresult' stringTime],'shifts','xout','yout','zout','z_err_out','zf_out','ll_out','t_out');
    
    save([currentfolder savename '_DCresult'],'shifts','xout','yout','zout','z_err_out','zf_out','ll_out','t_out','pr_out');
    programStatus_string =[programStatus_string;['Save DCresult to:' currentfolder savename '_DCresult']];
    set(handles.programStatus,'string',programStatus_string)
    drawnow update
    %% scaling based on fringe frequency
    
    % [zout]=iPALM_zscaling(zout,z_center,phaseslope,lambda,Nmedia);
    
    %% plot single layer image
    %     zflow=-1000;
    %     zfhigh=1000;
    %     zerrcut=100;
    %     llcut=600;
    %     for ss=1:1:numel(xout)
    %         maskzerr=z_err_out{ss}<zerrcut;
    %         maskzf=zout{ss}<zfhigh&zout{ss}>zflow;
    %         maskll=ll_out{ss}<llcut;
    %         mask2=maskzerr&maskzf&maskll;
    %
    %         tmpx=xout{ss}(mask2);
    %         tmpy=yout{ss}(mask2);
    %         tmpz=zout{ss}(mask2);
    %         tmpt=t_out{ss}(mask2);
    %
    %         [flag]=iPALM2vutarav2(currentresultpath,[savename '_3d_test_stack_60' num2str(ss-1) '_v18_driftcorrect'],1,{tmpx},{tmpy},{tmpz},{round(tmpt./10)});
    %     end
    
    %% align stack
    
    if numel(xout)==1
        shiftx_st=0;
        shifty_st=0;
        shiftz_st=0;
    else
        pixelsz=25;
        errorthresh=15;
        cutmeth='nocut';
        iniguess=[0 0 475];
        maskflag=0;
        stksort=[];
        [shiftx_st,shifty_st,shiftz_st,rankf]=iPALM_stack_RedunLSv9z(stksort,xout,yout,zout,pixelsz,errorthresh,cutmeth,iniguess,maskflag,[]);
    end
    % [shifty_st shiftz_st2 rankf]=iPALM_stack_RedunLSv6(stksort,ysort,zsort,pixelsz,errorthresh,cutmeth,iniguess,maskflag,[]);
    % [shiftx_st shifty_st rankf]=iPALM_stack_RedunLSv6(stksort,xsort,ysort,pixelsz,errorthresh,cutmeth,[0 0],maskflag,[]);
    [xcof]=shiftcoords_stack(xout,shiftx_st);
    [ycof]=shiftcoords_stack(yout,shifty_st);
    [zcof]=shiftcoords_stack(zout,shiftz_st);
    
    %% 2nd fileter result again
    zflow=-700;
    zfhigh=700;
    zerrcut=100;
    xoutf=[];
    youtf=[];
    zoutf=[];
    toutf=[];
    lloutf=[];
    Ioutf=[];
    zconf=[];
    crlbf=[];
    zasterf=[];
    bgf=[];
    prf=[];
    zerrf=[];
    
    for jj=1:1:numel(xcof)
        maskzerr=z_err_out{jj}<zerrcut;
        maskzf=zout{jj}<zfhigh&zout{jj}>zflow;
        mask2=maskzerr&maskzf;
        zerrf{jj}=z_err_out{jj}(mask2);
        xoutf{jj}=xcof{jj}(mask2);
        youtf{jj}=ycof{jj}(mask2);
        zoutf{jj}=zcof{jj}(mask2);
        toutf{jj}=t_out{jj}(mask2);
        Ioutf{jj}=I_out{jj}(mask2);
        lloutf{jj}=ll_out{jj}(mask2);
        zconf{jj}=zangctrl_out{jj}(mask2);
        crlbf{jj}=crlb_out{jj}(mask2);
        zasterf{jj}=Z_ast_err_out{jj}(mask2);
        bgf{jj}=bg_out{jj}(mask2);
        prf{jj}=pr_out{jj}(mask2);
    end
    
    
    %% output to Vutara    
    vutarax=cat(1,xoutf{:});
    vutaray=cat(1,youtf{:});
    vutaraz=cat(1,zoutf{:});
    vutarat=cat(1,toutf{:});
    vutaraI=cat(1,Ioutf{:});
    vutarall=cat(1,lloutf{:});
    vutarabg=cat(1,bgf{:});
    vutarazcon=cat(1,zconf{:});
    vutaracrlb=cat(1,crlbf{:});
    vutarazaster=cat(1,zasterf{:});
    vutarazerr=cat(1,zerrf{:});
    vutarapr=cat(1,prf{:});
    
    %     time = clock;
    %     stringTime = ['_' num2str(time(1)) '_' num2str(time(2)) '_' num2str(time(3)) '_' num2str(time(4)) '_' num2str(time(5)) '_' num2str(time(6)) '.mat'];
    %     save([currentresultpath savename '_642v20_60' stringTime],'vutarax','vutaray','vutaraz','vutarat');
    [flag]=iPALM2vutarav2(currentfolder,[savename '_ll'],1,{vutarax},{vutaray},{vutaraz},{ceil(vutarat/100)},{vutaraI},{vutaracrlb},{vutarall},{vutarabg},{vutarazcon},{vutarazerr});
    save([currentfolder savename '_' channel 'v20_60'],'vutarax','vutaray','vutaraz','vutarat','vutarall','vutaraI','vutarabg','vutarazcon','vutaracrlb','vutarazaster','vutarazerr','vutarapr');
    programStatus_string =[programStatus_string;['Save Vutarax File to:' currentfolder savename '_ll\' 'particles.csv']];
    set(handles.programStatus,'string',programStatus_string)
    drawnow update
    
    %% output image
    coords=[];
    coords(:,1)=vutarax/10;
    coords(:,2)=vutaray/10;
    zm=12.8;
    szx=168*zm;
    szy=168*zm;
    im=cHistRecon(szx,szy,single(coords(:,2)),single(coords(:,1)),0);
    gaussim=gaussf(im,[1 1]);
    str2=([currentfolder savename '_gauss_1.tif']);
    writeim(gaussim,str2,'tiff',1);
end