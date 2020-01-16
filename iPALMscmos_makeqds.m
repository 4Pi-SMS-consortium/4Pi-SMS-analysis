% for v14 iPALMast_analysisv14scmos.m

function [qds]=iPALMscmos_makeqds(ims,center,fsign)
vsz=size(ims,1);
% vsz=150;
qds=[];
for ii=1:1:numel(center)
    tmp=ims(:,(center(ii)-vsz/2):(center(ii)+vsz/2-1),:);
%     tmp=ims(54:203,(center(ii)-vsz/2):(center(ii)+vsz/2-1),:);
    if fsign(ii)==1
        tmp=flipdim(tmp,2);
    end
    qds(:,:,:,ii)=tmp;
end
    