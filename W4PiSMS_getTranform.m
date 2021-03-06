function W4PiSMS_getTranform(handles)

close all
posfile = get(handles.pathMainfolder, 'String');
center1 = str2num(get(handles.center1,'string'));
center2 = str2num(get(handles.center2,'string'));
center3 = str2num(get(handles.center3,'string'));
center4 = str2num(get(handles.center4,'string'));
centers = [center1 center2 center3 center4];

c = posfile(end); %get last character of path to image file
while c ~= '\' %get path only
    posfile(end) = []; %delete filename character
    if isempty(posfile) %break if no path present
        break;
    end
    c = posfile(end); %get last charakter of path to image file
end

[imageName,dataFolder] = uigetfile([posfile '*.dcimg'],'Open Calibration Images');
if imageName ~= 0
    if ~isempty(strfind(imageName,'_642_'))        
        resultpath = [pwd '\calibration_files\'];
        namestr = ['align' '_642_'];
        mkdir(resultpath);
        channel = '642';        
    elseif ~isempty(strfind(imageName,'_561_'))
        resultpath = [pwd '\calibration_files\'];
        namestr = ['align' '_561_'];
        mkdir(resultpath);
        channel = '561';
    else
        x = inputdlg({'Result Path','channel'},'File Directory',2,{[pwd '\calibration_files\'],''},options);
        resultpath = x{1};
        channel= x{2};
        namestr = ['align' '_' channel '_'];
        mkdir(resultpath);
    end
else
    return
end

[~,qds]=W4PiSMS_readdcimg([dataFolder imageName],centers);
meanqds1 = mean(qds(:,:,:,1),3);
meanqds2 = mean(qds(:,:,:,2),3);
meanqds3 = mean(qds(:,:,:,3),3);
meanqds4 = mean(qds(:,:,:,4),3);
qd1=qds(:,:,:,1)-median(meanqds1(:));
qd2=qds(:,:,:,2)-median(meanqds2(:));
qd3=qds(:,:,:,3)-median(meanqds3(:));
qd4=qds(:,:,:,4)-median(meanqds4(:));

%% find Fourier-Mellin transform
qdall=cat(4,qd1,qd2,qd3,qd4);
qdall(qdall<=0)=1e-37;
zm_all=[];
trans_all=[];
ang_all=[];
R=[];
invR=[];
for ii=1:1:4
    im1=mean(qdall(:,:,:,1),3);
    im2=mean(qdall(:,:,:,ii),3);
     
    [zm,trans,ang] = fmmatch(im2,im1);
    [out,R(:,:,ii)] = find_affine_trans(im2, im1, [[zm zm],trans,ang]);
    zm_fin=out(1:2);
    trans_fin=out(3:4);
    ang_fin=out(5);
    
    [imout]=affine_trans(im2,zm_fin,trans_fin,ang_fin);
    zm_all(ii,:)=zm_fin;
    trans_all(ii,:)=trans_fin;
    ang_all(ii,:)=ang_fin;
    joinchannels('RBG',im1,imout)
    
    [zm,trans,ang] = fmmatch(im1,im2);
    [~,invR(:,:,ii)] = find_affine_trans(im1, im2, [[zm zm],trans,ang]);
    set(handles.programStatus, 'String', ['Apply Affine Transform to Channel ' num2str(ii) ' of 4']); %show progress
    drawnow update;
end
para.centers = centers;
para.file = [dataFolder imageName];
datestring = datestr(now,'yyyymmdd');
save([resultpath namestr 'FMTtransform_' datestring],'zm_all','trans_all','ang_all','R','invR','para');
I=find(dataFolder=='\',2,'last');
parentFolder = dataFolder(1:I);
save([parentFolder namestr 'FMTtransform_' datestring],'zm_all','trans_all','ang_all','R','invR','para');
set(handles.programStatus, 'String', {'Channel Alignment', ['Color Channel:' channel],['Input:' dataFolder imageName],['Output:' parentFolder namestr 'FMTtransform' datestring '.mat']}); %show progress
drawnow update; %show progress
