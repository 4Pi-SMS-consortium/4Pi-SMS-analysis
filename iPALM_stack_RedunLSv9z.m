function [R_shift1 R_shift2 R_shift3 flag error1_final error2_final error3_final]=iPALM_stack_RedunLSv9z(stksort,x1in,x2in,x3in,pixelsz,errothresh,cutmeth,iniguess,maskflag,cropsz)

flag=[0 0];
if nargin<6
    pixelsz=15;
    errothresh=15; % in nm
end

if nargin<7
    inipx=[0 0];
else
    inipx=round(iniguess./pixelsz);
end

rndsz=2;
thresh=errothresh/pixelsz;
% parse the data

minx1=min(cat(1,x1in{:}));
minx2=min(cat(1,x2in{:}));
minx3=min(cat(1,x3in{:}));

x1=cellfun(@(x) (x-minx1)./pixelsz, x1in,'UniformOutput',false);
x2=cellfun(@(x) (x-minx2)./pixelsz, x2in,'UniformOutput',false);
x3=cellfun(@(x) (x-minx3)./pixelsz, x3in,'UniformOutput',false);

% x1=(x1in-min(cat(1,x1in{:})))./pixelsz;
% x2=(x2in-min(cat(1,x2in{:})))./pixelsz;
% x3=(x3in-min(cat(1,x3in{:})))./pixelsz;
% x1sz=ceil(max(x1)/rndsz)*rndsz;
% x2sz=ceil(max(x2)/rndsz)*rndsz;

x1sz=ceil(max((cat(1,x1{:})))./rndsz)*rndsz+1;
x2sz=ceil(max((cat(1,x2{:})))./rndsz)*rndsz+1;
x3sz=ceil(max((cat(1,x3{:})))./rndsz)*rndsz+1;

% x3sz=min(59,max(25,x3sz));

if nargin<=5
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
%     bomask1=x1>=1e-37;
%     bomask2=x2>=1e-37;
end


% determine segments
% seg=unique(stksort);
segnum=numel(x1in);
% AA=zeros(segnum*(segnum-1)/2,segnum-1);
kk=1;
for ii=1:1:segnum
%     maskt=(stksort==seg(ii));
    
    x1_1=x1{ii};
    x2_1=x2{ii};
    x3_1=x3{ii};
%     im1=cHistRecon(x1sz,x2sz,single(x1_1),single(x2_1),1);
    im1=cHistRecon3D(x1sz,x2sz,x3sz,single(x1_1),single(x2_1),single(x3_1),0);
%     if maskflag==1
%         %         imref=cHistRecon(x1sz,x2sz,single(x1_1),single(x2_1),1);
%         finishflag=0;
%         while finishflag==0
%             [pickco]=pickspot(gaussf(im1,1));
%             if size(pickco,1)<2
%                 cropsz=512;
%             end
%             
%             if isempty(cropsz)
%                 cropsz=2*max(abs(pickco(1,:)-pickco(2,:)));
%             end
%             
%             try
%                 finishflag=1;
%                 im1=cropim(im13d,[[pickco(1,:)]-round(cropsz/2) 0],[cropsz cropsz x3sz]);
%                 if sum(im1(:))==0||isempty(im1)
%                     finishflag=0;
%                 end
%             catch
%                 finishflag=0;
%             end
%         end
%         
%     end
    
    for jj=ii+1:1:min(ii+1,segnum)
%         maskt=(stksort==seg(jj));
        % calculate permuted pairwise distances
        %         st2=(jj-1)*frmnum+min(toutin);
        %         ed2=(jj)*frmnum+min(toutin)-1;
        %         maskt=toutin>=st2&toutin<=ed2;
        
        x1_2=x1{jj};
        x2_2=x2{jj};
        x3_2=x3{jj};
        
        im2=cHistRecon3D(x1sz,x2sz,x3sz,single(x1_2),single(x2_2),single(x3_2),0);
        
%         if maskflag==1
%             im2= cropim(im23d,[[pickco(1,:)]-round(cropsz/2) 0],[cropsz cropsz x3sz]);
%         end
%         corrim3d=normxcorr3(double(gaussf(im1,1)),double(gaussf(im2,1)),'full');
%         [shifts]=findshift(gaussf(im1,4),gaussf(im2,4),'iter');
        [shift1(kk) shift2(kk) shift3(kk)]=drift_correction_core3d(double(im1),double(im2),inipx*(jj-ii));
        AA(kk,ii:jj-1)=1;
        kk=kk+1;
    end
end

R_shift1=pinv(AA)*shift1';
R_shift2=pinv(AA)*shift2';
R_shift3=pinv(AA)*shift3';

error1=sqrt(((AA*R_shift1)-shift1').^2);
error2=sqrt(((AA*R_shift2)-shift2').^2);
error3=sqrt(((AA*R_shift3)-shift3').^2);

newAA=AA;
oldAA=newAA;
oldshift=shift1;
while rank(newAA)>=min(size(AA))&&max(error1)>thresh
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
    flag(1)=1;
end

R_shift1=pinv(newAA)*shift1';
error1_final=sqrt(((newAA*R_shift1)-shift1').^2);

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
    flag(2)=1;
end

error2_final=sqrt(((newAA*R_shift2)-shift2').^2);
R_shift2=pinv(newAA)*shift2';

newAA=AA;
oldAA=newAA;
oldshift=shift3;
while rank(newAA)>=min(size(AA))&&max(error3)>thresh
    oldAA=newAA;
    oldshift=shift3;
    [maxval ind]=max(error3);
    newAA(ind,:)=[];
    error3(ind)=[];
    shift3(ind)=[];
end

if rank(newAA)<min(size(AA))
    newAA=oldAA;
    shift3=oldshift;
    flag(1)=1;
end

R_shift3=pinv(newAA)*shift3';
error3_final=sqrt(((newAA*R_shift3)-shift3').^2);

R_shift1=R_shift1*pixelsz;
R_shift2=R_shift2*pixelsz;
R_shift3=R_shift3*pixelsz;