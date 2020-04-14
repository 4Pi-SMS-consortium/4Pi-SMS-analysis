function W4PiSMS_getPosition(handles)

mainFolder = get(handles.pathMainfolder,'string');
index_selected = get(handles.pathSubfolder,'value');
folderlist = get(handles.pathSubfolder,'string');
folder_selected = folderlist(index_selected);

%parameters
det_thresh=str2num(get(handles.det_thresh,'string'));
CRLB_thresh=0;
sigma_thresh=0;
llrthresh=0;
center1 = str2num(get(handles.center1,'string'));
center2 = str2num(get(handles.center2,'string'));
center3 = str2num(get(handles.center3,'string'));
center4 = str2num(get(handles.center4,'string'));

centers = [center1 center2 center3 center4];
para.det_thresh = det_thresh;
para.CRLB_thresh = CRLB_thresh;
para.sigma_thresh = sigma_thresh;
para.llrthresh = llrthresh;
para.centers = centers;
para.imagesz = 168;   % image size, nm
para.pixel = 128;     % pixel size, nm

for ii=1:numel(folder_selected)
    close all;
    
    currentfolder=folder_selected{ii};
    I=find(currentfolder=='\',1,'last');
    childfoldername=currentfolder(I+1:end);
    
    I=find(currentfolder=='\',1,'last');
    parentFolder = currentfolder(1:I);
    
    files=dir([currentfolder '\*.dcimg']);
    imageName = files(1).name;
    if ~isempty(strfind(imageName,'_642_'))
        channel = '642';
    elseif ~isempty(strfind(imageName,'_561_'))
        channel = '561';
    else
        x = inputdlg({'Channel'},'Please specify the channel',1,{'488'});
        channel = x{1};        
    end
    para.channel = channel;
    
    %% 4PiSMS sanalysis
    foldername=[currentfolder '\'];
    nameroot=['*.dcimg'];
    savename=childfoldername;    
    
    tmpf=dir([mainFolder '*_' channel '_*Ast*.mat']);
    if numel(tmpf)>1
        set(handles.programStatus,'string','More than one *Ast*.mat calibration file detected!');
        return;
    elseif numel(tmpf)==0
        set(handles.programStatus,'string','No *Ast*.mat calibration file detected!');
        return;
    end    
    astfile=tmpf.name;
    astfile=[mainFolder astfile];
    
    tmpf=dir([mainFolder '*_' channel '_*dphi*.mat']);    
    if numel(tmpf)>1
        set(handles.programStatus,'string','More than one *dphi*.mat calibration file detected!');
        return;
    elseif numel(tmpf)==0
        set(handles.programStatus,'string','No *dphi*.mat calibration file detected!');
        return;
    end    
    anglefile=tmpf.name;
    anglefile=[mainFolder anglefile];
    
    tmpf=dir([mainFolder '*_' channel '_*FMTtransform*.mat']);    
    if numel(tmpf)>1
        set(handles.programStatus,'string','More than one *FMTtransform*.mat calibration file detected!');
        return;
    elseif numel(tmpf)==0
        set(handles.programStatus,'string','No *FMTtransform*.mat calibration file detected!');
        return;
    end    
    fmfile=tmpf.name;
    fmfile=[mainFolder fmfile];
       
    datestring = datestr(now,'yyyymmdd');
    set(handles.programStatus,'string',{'Current Analized Folder:',currentfolder, ['Current Channel:' channel],'Calibration Files Used:',astfile,anglefile,fmfile,'Output:',[mainFolder savename '_' channel '_tmpresult_' datestring '.mat']})
    drawnow update
    
    [xresult yresult tresult zfresult zangresult llresult CRLBresult Iresult ...
        zast_err_result stacktot num_images zangctrresult bgresult subimstot  sxtot sytot...
        ]=W4PiSMS_RM_analysisv8stack(foldername,nameroot,fmfile,astfile,anglefile,llrthresh,handles.scmos_cali_file,det_thresh,CRLB_thresh,sigma_thresh,centers,handles);
    
    para.num_images=num_images;
    save([parentFolder savename '_' channel '_tmpresult_' datestring],'foldername','sxtot','sytot','tresult','anglefile','astfile','xresult','yresult','zfresult',...
    'zangctrresult','bgresult','zangresult','llresult','CRLBresult','Iresult','zast_err_result','stacktot','para');
    
    %% output image
    coords=[];
    pixel=para.pixel;  
    pixel_SR=10;          % nm      
    coords(:,1)=xresult*pixel/pixel_SR;
    coords(:,2)=yresult*pixel/pixel_SR;
    sz=para.imagesz*pixel/pixel_SR;
    im=cHistRecon(sz,sz,single(coords(:,2)),single(coords(:,1)),0);
    gaussim=gaussf(im,[1 1]);
    str2=([parentFolder savename '_gauss_1.tif']);
    writeim(gaussim,str2,'tiff',1);   
end