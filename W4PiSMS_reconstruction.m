function W4PiSMS_reconstruction(handles)

mainFolder = get(handles.pathMainfolder,'string');
index_selected = get(handles.pathSubfolder,'value');
filelist = get(handles.pathSubfolder,'string');
file_selected = filelist(index_selected);

photons=str2double(get(handles.photons,'string'));
llrthresh=str2double(get(handles.llrthresh2,'string'));
crlb_thresh = str2double(get(handles.crlb,'string'));

angctr=[0.4 1.2];
st_frame=0;             %startframe
frmthresh=3000*1000;    %stopframe
stopval=str2double(get(handles.stopValue,'string'));     %stop value for Dmap
zthresh=[0 1200];
zasterr=[0 0.05];
reverseflag=0; %642 ch does need to be reversed. from beginning to the end.

%%
freqld=load('frequency_simulate4Pi_594_690.mat');
for ff=1:numel(file_selected)

    currentfile=file_selected{ff};
    I=find(currentfile=='\',1,'last');
    currentfolder=currentfile(1:I);
    fileName=currentfile(I+1:end);
    if contains(fileName,'_642_')
        channel = '642';
        reverseflag=0; 
        freqstr='f_690';       
        eval(['freq=freqld.' freqstr ';']);
    elseif contains(fileName,'_561_')
        channel = '561';
        reverseflag=0;
        freqstr='f_594';
         eval(['freq=freqld.' freqstr ';']);
    else
        channel = 'not 561 or 642'
    end
    savename=fileName(1:end-27);%name format foldername_642_tempresult_yymmdd.mat
    datestring = datestr(now,'yyyymmdd');
    set(handles.programStatus,'string',{'Current Analized File:',currentfile, ['Current Channel:' channel]})
    drawnow update
    tmpld=load(file_selected{ff});
    imagesz=tmpld.para.imagesz;
    pixelsz=tmpld.para.pixelsz;  
    num_images=tmpld.para.num_images;
    
    %% reject beads
    Vin=[];
    Vin(:,1)=tmpld.xresult;
    Vin(:,2)=tmpld.yresult;
    Vin(:,3)=tmpld.tresult+1;
    maskbead=reject_beads(single(Vin),tmpld.para.num_images,imagesz);
    
    %% rejection   
    mtc=(tmpld.sxtot.^3./tmpld.sytot-tmpld.sytot.^3./tmpld.sxtot)./40*2*pi;   
    % mtc mask
    maskmtc=mtc<pi&mtc>-pi &(tmpld.sxtot+tmpld.sytot)<6;
    % ll threshold
    maskll=tmpld.llresult<llrthresh;
    % crlb
    maskcrlb=tmpld.CRLBresult(:,1)<crlb_thresh & tmpld.CRLBresult(:,2)<crlb_thresh;
    % frame number threshold
    maskt=tmpld.tresult>=st_frame & tmpld.tresult<frmthresh;
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
    % apply
    mask=maskI&maskll&maskzctr&maskmtc&maskcrlb&maskt&maskbead;%&maskxy;%&maskroi;
    
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

    %% dmap test
    segnum=ceil((max(tresult))/num_images);
    zest=[];
    zerr=[];
    zmask=[];
    Dmap0=[];
    Dmap=[];
    for ii=1:1:segnum
        close all
        st=(ii-1)*num_images;
        if ii==segnum
            ed=max(tresult);
        else
            ed=(ii)*num_images-1;
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
    
    str=([currentfolder savename '_dmap.tif']);
    tiffwrite(Dmap*10000,str);
    programStatus_string = get(handles.programStatus,'string');
    programStatus_string = [programStatus_string;['Save Dmap Image to:' str]];
    set(handles.programStatus,'string',programStatus_string)
    
    %% phase estimate for every optical section separately
    unistack=unique(stacktot);
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
        
        if numel(unistack)>1
            currt=currt-(ss-1)*num_images;
            cycle=floor(currt/(num_images*numel(unistack)));
            remn=rem(currt,num_images*numel(unistack));
            currt=cycle*num_images+remn;
        end
        display(['processing segment:' num2str(ss)]);
        
        centermtc=[];
        stopval1=stopval;
        [currzresult,z_err,mephimask]=Mephi_z_4PiSMS(currzang,currmtc,currt,num_images,1,stopval1,centermtc,freq);

        maskzerr=z_err<100;        
        maskall=maskzerr&(mephimask>0)&abs(currzresult)<700;
        currzresult=currzresult(maskall>0);        
        currx=currx(maskall>0);
        curry=curry(maskall>0);
        currt=currt(maskall>0);
            
        z_err=z_err(maskall>0);
        currzfresult=currzfresult(maskall>0);
        currll=currll(maskall>0);
        currI=currI(maskall>0);
        currzangctrl=currzangctrl(maskall>0);
        currcrlb=currcrlb(maskall>0);
        currzasterr=currzasterr(maskall>0);
        currbg=currbg(maskall>0);
        currpr=currpr(maskall>0);
        
        driftstr=[currentfolder,savename,'_',num2str(ss),'_drift.mat'];
        if length(unistack)>1
            interpflag=0;
        else
            interpflag=1;
        end
        
        [xout{ss},yout{ss},zout{ss},shifts{ss}]=W4PiSMS_driftcorrection_RedunLSv10(currx.*pixelsz,curry.*pixelsz,currzresult,currt,num_images,reverseflag,interpflag,driftstr);
        
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
    
    save([currentfolder savename '_DCresult'],'shifts','xout','yout','zout','z_err_out','zf_out','ll_out','t_out','pr_out');
    programStatus_string =[programStatus_string;['Save DCresult to:' currentfolder savename '_DCresult']];
    set(handles.programStatus,'string',programStatus_string)
    drawnow update   
    
    %% align stack
    if numel(xout)==1
        shiftx_st=0;
        shifty_st=0;
        shiftz_st=0;
    else
        pixelsz=25;
        errorthresh=15;
        cutmeth='nocut';
        iniguess=[0 0 500];
        maskflag=0;
        stksort=[];
        [shiftx_st,shifty_st,shiftz_st,rankf]=W4PiSMS_stack_RedunLSv9z(stksort,xout,yout,zout,pixelsz,errorthresh,cutmeth,iniguess,maskflag,[]);
    end
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
    
    [flag]=W4PiSMS2vutarav2(currentfolder,[savename '_ll'],1,{vutarax},{vutaray},{vutaraz},{vutarat},{vutaraI},{vutaracrlb},{vutarall},{vutarabg},{vutarazcon},{vutarazerr});
    save([currentfolder savename '_' channel 'v20_60'],'vutarax','vutaray','vutaraz','vutarat','vutarall','vutaraI','vutarabg','vutarazcon','vutaracrlb','vutarazaster','vutarazerr','vutarapr');
    programStatus_string =[programStatus_string;['Save Vutarax File to:' currentfolder savename '_ll\' 'particles.csv']];
    set(handles.programStatus,'string',programStatus_string)
    drawnow update
    
    %% output image
    coords=[];
    pixel_SR=10;    
    coords(:,1)=vutarax/pixel_SR;
    coords(:,2)=vutaray/pixel_SR;
    sz=imagesz*pixelsz/pixel_SR;
    im=cHistRecon(sz,sz,single(coords(:,2)),single(coords(:,1)),0);
    gaussim=gaussf(im,[1 1]);
    str2=([currentfolder savename '_gauss_1.tif']);
    writeim(gaussim,str2,'tiff',1);
end