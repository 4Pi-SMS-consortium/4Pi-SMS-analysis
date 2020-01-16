% for v14 iPALMast_analysisv14scmos.m

% Version numbers:
% July30th 2014: updated with recent changes in FMTtransform format

function [q1 q2 q3 q4]=iPALMast_RotAlign_FMT(qd1,qd2,qd3,qd4,fmtname)
tmp=load(fmtname);
zm=tmp.zm_all;
trans=tmp.trans_all;
ang=tmp.ang_all;

for ii=1:1:size(qd1,3)
    
q1(:,:,ii)=qd1(:,:,ii);

q2(:,:,ii)=double(affine_trans(qd2(:,:,ii),zm(2,:),trans(2,:),ang(2)));

q3(:,:,ii)=double(affine_trans(qd3(:,:,ii),zm(3,:),trans(3,:),ang(3)));

q4(:,:,ii)=double(affine_trans(qd4(:,:,ii),zm(4,:),trans(4,:),ang(4)));

if mod(ii,1000)==0
    display([num2str(ii) ' out of ' num2str(size(qd1,3)) ' is done...']);
end
end