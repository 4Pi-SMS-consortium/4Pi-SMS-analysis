% for v14 iPALMast_analysisv14scmos.m

% patch notes:
% Dec23rd 2014: change to include only 4 segements for drift correction.
% Reason: large drift observed and eventually out of the boundary of fft
% max finding routine.
function [R_shift1 R_shift2]=correct_drift_LS(x1in,x2in,toutin,frmnum,pixelsz,errothresh,cutmeth)

if nargin<6
    pixelsz=15;
    errothresh=7.5; % in nm
end

rndsz=2;
thresh=errothresh/pixelsz;
% parse the data
x1=floor((x1in-min(x1in))./pixelsz);
x2=floor((x2in-min(x2in))./pixelsz);

x1sz=max(24,ceil(max(x1)/rndsz)*rndsz+1);
x2sz=max(24,ceil(max(x2)/rndsz)*rndsz+1);
display(['drift_correction_reconstruction size: ' num2str(x1sz) ', ' num2str(x2sz)]);

if nargin<=6
    bomask1=(x1>=max(x1)*0.2)&(x1<=max(x1)*0.8);
    bomask2=(x2>=max(x2)*0.2)&(x2<=max(x2)*0.8);
elseif strcmp(cutmeth,'percentile')
    p1_20 = prctile(x1,20);
    p1_80 = prctile(x1,80);
    p2_20 = prctile(x2,20);
    p2_80 = prctile(x2,80);
    bomask1=(x1>=p1_20)&(x1<=p1_80);
    bomask2=(x2>=p2_20)&(x2<=p2_80);
elseif strcmp(cutmeth,'nocut')
    bomask1=x1>=1e-37;
    bomask2=x2>=1e-37;
end


% determine segments
segnum=ceil((max(toutin)+1)/frmnum);
% AA=zeros(segnum*(segnum-1)/2,segnum-1);
kk=1;
for ii=1:1:segnum-1
    st1=(ii-1)*frmnum;
    ed1=(ii)*frmnum-1;
    maskt=toutin>=st1&toutin<=ed1;
    
    x1_1=x1(maskt&bomask1&bomask2);
    x2_1=x2(maskt&bomask1&bomask2);
    im1=cHistRecon(x1sz,x2sz,x1_1,x2_1,0);
    if sum(im1(:))==0
        continue
    end
%     max(x1_1)
%     max(x2_1)
    for jj=(ii+1):1:min(segnum,(ii+1+10))
        % calculate permuted pairwise distances
        st2=(jj-1)*frmnum;
        ed2=(jj)*frmnum-1;
        maskt=toutin>=st2&toutin<=ed2;
        
        x1_2=x1(maskt&bomask1&bomask2);
        x2_2=x2(maskt&bomask1&bomask2);
        
        im2=cHistRecon(x1sz,x2sz,x1_2,x2_2,0);
        
        [shift1(kk) shift2(kk)]=drift_correction_core(double(gaussf(im1,1)),double(gaussf(im2,1)),[0 0]);
        AA(kk,ii:jj-1)=1;
        kk=kk+1;
    end
end

R_shift1=pinv(AA)*shift1';
R_shift2=pinv(AA)*shift2';

error1=sqrt(((AA*R_shift1)-shift1').^2);
error2=sqrt(((AA*R_shift2)-shift2').^2);

newAA=AA;
oldAA=newAA;
oldshift=shift1;
while (rank(newAA)>=min(size(AA)))&&(max(error1)>thresh)
    oldAA=newAA;
    oldshift=shift1;
    [maxval ind]=max(error1);
    newAA(ind,:)=[];
    error1(ind)=[];
    shift1(ind)=[];
end

if rank(newAA)<min(size(AA))
    newAA=oldAA;
    shift1=oldshift;
end

R_shift1=pinv(newAA)*shift1';


newAA=AA;
oldAA=newAA;
oldshift=shift2;
while (rank(newAA)>=min(size(AA)))&&(max(error2)>thresh)
    oldAA=newAA;
    oldshift=shift2;
    [maxval ind]=max(error2);
    newAA(ind,:)=[];
    error2(ind)=[];
    shift2(ind)=[];
end

if rank(newAA)<min(size(AA))
    newAA=oldAA;
    shift2=oldshift;
end

R_shift2=pinv(newAA)*shift2';

R_shift1=R_shift1*pixelsz;
R_shift2=R_shift2*pixelsz;