function [ims,qds,cocrops]=W4PiSMS_readdcimg(filename,center,co_cropims)

files=dir(filename);
if numel(files)>1
    error('Multiple files with the same assigned name are detected');
elseif numel(files)==0
    error('No file detected');
end

[tmp,totalframes]= dcimgmatlab(0, filename);
ims=single(zeros(size(tmp,2),size(tmp,1),totalframes));
cutflag=0;
cocropflag=0;

qds=[];
cocrops=[];
if nargin>1
    cutflag=1;
end

if nargin>2
    cocropflag=1;
end

if cutflag==0
    for ii=1:1:totalframes
        ims(:,:,ii)=single(dcimgmatlab(ii-1, filename))';
    end
else
    for ii=1:1:totalframes
        ims(:,:,ii)=single(dcimgmatlab(ii-1, filename))';
    end
    
    flipsigns=[0 0 1 1];
    [qds]=single(W4PiSMS_scmos_makeqds(ims,center,flipsigns));
    ims=[];
    
    if cocropflag==1
        [cocrops]=single(W4PiSMS_scmos_makeqds(co_cropims,center,flipsigns));
    end   
end

clear mex

