function iPALMgetPhaseShift(handles)

close all
center1 = str2num(get(handles.center1,'string'));
center2 = str2num(get(handles.center2,'string'));
center3 = str2num(get(handles.center3,'string'));
center4 = str2num(get(handles.center4,'string'));
centers = [center1 center2 center3 center4];

posfile = get(handles.pathMainfolder, 'String');
c = posfile(end); %get last character of path to image file
while c ~= '/' %get path only
    posfile(end) = []; %delete filename character
    if isempty(posfile) %break if no path present
        break;
    end
    c = posfile(end); %get last charakter of path to image file
end

[imageName,dataFolder] = uigetfile([posfile '*.dcimg'],'Open Calibration Images','MultiSelect','on');
I=find(dataFolder=='/',2,'last');
parentFolder = dataFolder(1:I);
% clear; clc;close all
% addpath '../calibration_files/'
% addpath '../'
% %% directory setup
% datafolder=['G:\4PISCMOS\Mar23rd-2016\642\Cell03\'];
% nameroot=['Cell03_642*.dcimg'];

if ~isequal(imageName,0)
    if ~isempty(strfind(imageName{1},'_642_'))
        channel = '642';
    elseif ~isempty(strfind(imageName{1},'_561_'))
        channel = '561';
    else
        x = inputdlg({'Channel'},'Please specify the channel',1,{'488'});
        channel = x{1};
    end
else
    return
end

tmpf=dir([parentFolder '*_' channel '_*FMTtransform*.mat']);
if numel(tmpf)~=1
    [FMTtransformName,FMTtransformFolder] = uigetfile([posfile '*.mat'],'Open FMTtransform');
    if isempty(FMTtransformName)
        return
    else
        fmtname = [FMTtransformFolder FMTtransformName];
        
    end
else
    fmtname = [parentFolder tmpf.name];
end
resultpath = [pwd '/calibration_files/'];
namestr= ['bead_' channel];
mkdir(resultpath);

%%
tmpld=load(handles.scmos_cali_file);
offsetim=tmpld.offsetim; % for scmos, yes for now
varim=tmpld.varim;
gainim=tmpld.g3;

% offsetim=ones(164,1108)*100;
% varim=ones(164,1108)*13;
% gainim=ones(164,1108)*2.2;

caliims=cat(3,offsetim,varim,gainim);
% offsetim=ones(256,2048)*110;
% varim=ones(256,2048);
% gainim=ones(256,2048)*2;
% caliims=cat(3,offsetim,varim,gainim);
% offset=110;
%     gain=1;
% resultpath='..\calibration_files\';
%     mkdir(resultpath);
% totfiles=dir([datafolder nameroot]);

aveqd1=[];
aveqd2=[];
aveqd3=[];
aveqd4=[];
set(handles.programStatus,'string','Read Beads Images')
drawnow update
for ii=1:1:numel(imageName)
    [~,qds,calicrops]=iPALM_readdcimg([dataFolder imageName{ii}],centers,caliims);
    tmpim=calicrops(:,:,3,:);%gain
    tmpim(tmpim<1.3|tmpim>3.5)=mean(mean(mean(calicrops(:,:,3,:))));
    calicrops(:,:,3,:)=tmpim;
    offsetim=squeeze(calicrops(:,:,1,:));
    varim=squeeze(calicrops(:,:,2,:));
    gainim=squeeze(calicrops(:,:,3,:));
    aveqd1(:,:,ii)=(mean(qds(:,:,:,1),3)-offsetim(:,:,1))./gainim(:,:,1);
    aveqd2(:,:,ii)=(mean(qds(:,:,:,2),3)-offsetim(:,:,2))./gainim(:,:,2);
    aveqd3(:,:,ii)=(mean(qds(:,:,:,3),3)-offsetim(:,:,3))./gainim(:,:,3);
    aveqd4(:,:,ii)=(mean(qds(:,:,:,4),3)-offsetim(:,:,4))./gainim(:,:,4);
%     aveqd1(:,:,ii)=mean(qds(:,:,:,1),3);
%     aveqd2(:,:,ii)=mean(qds(:,:,:,2),3);
%     aveqd3(:,:,ii)=mean(qds(:,:,:,3),3);
%     aveqd4(:,:,ii)=mean(qds(:,:,:,4),3);
end

%% rotate and align
[q1 q2 q3 q4]=iPALMast_RotAlign_FMT(aveqd1,aveqd2,aveqd3,aveqd4,fmtname);
[v1 v2 v3 v4]=iPALMast_RotAlign_FMT(varim(:,:,1),varim(:,:,2),varim(:,:,3),varim(:,:,4),fmtname);
[g1 g2 g3 g4]=iPALMast_RotAlign_FMT(gainim(:,:,1),gainim(:,:,2),gainim(:,:,3),gainim(:,:,4),fmtname);
maskg=g1<=1|g2<=1|g3<=1|g4<=1;
v1(maskg)=1e7;
v2(maskg)=1e7;
v3(maskg)=1e7;
v4(maskg)=1e7;

g1(maskg)=1e-7;
g2(maskg)=1e-7;
g3(maskg)=1e-7;
g4(maskg)=1e-7;

colorim=joinchannels('RGB',q3,q1);
colorim2=joinchannels('RGB',q4,q2);
%% find local centers
for kk=1:1
    sumim1=q1+q4+q2+q3;
%     eval(['sumim1=q',num2str(kk),';']);
    sumv=v1+v2+v3+v4;
    sumvg=v1./g1./g1+v2./g2./g2+v3./g3./g3+v4./g4./g4;
    sumim1(sumim1<=1e-6)=1e-6;
    subsz=21;
    cutim=[];
    subvar_g=[];
    cut1=[];
    cut2=[];
    cut3=[];
    cut4=[];
    N=size(q1,3);
    % [coords]=hand_pick(sumim1(:,:,50));
    h=dipshow(sumim1(:,:,round(N/2)),'lin');
    set(handles.programStatus,'string','Select Bead')
    drawnow update
    if kk==1
        coords = dipgetcoords;
        coords = coords(1:2);
    end
    tempim = cut(sumim1(:,:,round(N/2)),[subsz subsz],coords-10);
    dipshow(tempim,'lin');
    diptruesize('off')
    pause(0.5)
    set(handles.programStatus,'string','Confirm Selection')
    drawnow update
    choice = questdlg('Confirm Selection?', ...
        '', ...
        'Yes','No','No');
    while strcmp(choice,'No')
        dipshow(h,sumim1(:,:,round(N/2)),'lin');
        coords = dipgetcoords;
        coords = coords(1:2);
        tempim = cut(sumim1(:,:,round(N/2)),[subsz subsz],coords-10);
        dipshow(tempim,'lin');
        diptruesize('off')
        pause(0.5)
        choice = questdlg('Confirm Selection?', ...
            '', ...
            'Yes','No','No');
    end
    %     coords=[84,84];
    % meanim=mean(sumim1,3);
    % [cx,cy]=find(meanim==max(meanim(:)));
    % coords=[cy,cx];
    
    for ii=1:1:size(sumim1,3)
        tmp=sumim1(:,:,ii);
        cutim(:,:,ii)=cut(tmp,[subsz subsz],coords-10);
        subvar_g(:,:,ii) = cut(sumvg,[subsz subsz],coords-10);
        cut1(:,:,ii)=cut(q1(:,:,ii),[subsz subsz],coords-10);
        cut2(:,:,ii)=cut(q2(:,:,ii),[subsz subsz],coords-10);
        cut3(:,:,ii)=cut(q3(:,:,ii),[subsz subsz],coords-10);
        cut4(:,:,ii)=cut(q4(:,:,ii),[subsz subsz],coords-10);
    end
    set(handles.programStatus,'string','Fit the Beads')
    drawnow update
    %% fit sigmax sigmay
    
    % [x1 y1]=iPALM_radial_symmetry(cut1);
    % [P CRLB LL]=gaussmlev2(single(cutim),2,200,4);
    [P CRLB LL]=sCMOS_MLE_CUDA42(single(cutim),2,200,4,single(subvar_g),single(subvar_g.*0+1));
    
    set(handles.programStatus,'string','Plot and fit Sigmax and SigmaY')
    drawnow update
    
    xstep=(1:1:length(P(:,1)))'*20;
    xco=P(:,2);
    yco=P(:,1);
    I=P(:,3);
    llr=-2*LL;
    sigmax=P(:,5);
    sigmay=P(:,6);
    
    % plot(xstep,sigmax,'r*');
    % hold on
    % plot(xstep,sigmay,'b*');
    
    sx=P(:,5);
    sy=P(:,6);
    A=(sx.^3./sy-sy.^3./sx)./40*2*pi;
    
    S=abs(sx-sy);
    id=find(S==min(S))
end
%% generate spline fit

stacksz=1;
sigxmean=mean(reshape(sigmax,[stacksz numel(sigmax)/stacksz]),1);
sigymean=mean(reshape(sigmay,[stacksz numel(sigmay)/stacksz]),1);
zstep=(1:1:length(sigxmean))*20;
sigxp=spline(zstep,sigxmean);
sigyp=spline(zstep,sigymean);

%%
figure
[estx model]=fit_ast(xstep,sigmax);
[~,fit_curve]=model(estx);
plot(xstep,fit_curve,'b--','linewidth',4);
hold on
plot(xstep,sigmax,'bo');

[esty model]=fit_ast(xstep,sigmay);
[~,fit_curve]=model(esty);
plot(xstep,fit_curve,'r--','linewidth',4);
plot(xstep,sigmay,'ro');

xlabel('relative z (nm)');
ylabel('sigma');
legend('Obs. \sigma_x','Fit. \sigma_x','Fit. \sigma_y','Location','North');
xlim([0 1200]);
% saveas(gcf,[resultpath namestr 'sigmaxy_fit_z_top.png'],'png');
%% save
dat = now;
datestring=datestr(now,'yyyymmdd');
save([resultpath namestr '_Astfit_' datestring],'model','estx','zstep','esty','sigmax','sigmay','dat');
save([parentFolder namestr '_Astfit_' datestring],'model','estx','zstep','esty','sigmax','sigmay','dat');

%% the following part of code is incopperated in dphi file with name '_dev_find_phaseshift'
set(handles.programStatus,'string','Find phase shift')
drawnow update
%% get rms rmp
xf=xco;
yf=yco;
subims=cat(4,cut1,cut2,cut3,cut4);
[rms,rmp]=iPALMast_findmom_givenC(subims,xf,yf);

%% calculate reduced moment
rm1=rms./sqrt(rms.^2+rmp.^2);
rm2=rmp./sqrt(rms.^2+rmp.^2);
% rm1=rms;
% rm2=rmp;

figure;
plot(rm1,'r');
hold on;
plot(rm2,'b');
plot(-rm1,'--r');
plot(-rm2,'--b');
legend('rm1','rm2','inv_rm1','inv_rm2');

%% estimate the phase delay
%
A1=0.9;
w=0.43;
phi1=1.4;
b1=0.1;
A2=0.9;
phi2=phi1-1.65;
b2=b1;
% id=27;
rmsin=rm1(id-20:id+20);
rmpin=rm2(id-20:id+20);
% rmsin=rm1;
% rmpin=rm2;
x=1:1:numel(rmsin);

% fit indivual curve first
[est1 mode1l]=iPALM_find_delay_single(x,rmsin,[A1 w phi1 b1]);
[est2 model2]=iPALM_find_delay_single(x,rmpin,[A2 w phi2 b2]);

figure
plot(x,rmsin,'r');
hold on
plot(x,rmpin,'b');
start_point=[est1(1,1) est1(1,2) est1(1,3) est1(1,4) est2(1,1) est2(1,3) est2(1,4)];
[est model]=iPALM_find_delay(x,rmsin,rmpin,start_point);

[~,fitrms,fitrmp]=feval(model,est);

plot(x,fitrms,'--r');
hold on
plot(x,fitrmp,'--b');

% phase difference is phirms=est(3), phirmp=est(6)

dphi=est(3)-est(6)
phi_s=est(3);
phi_p=est(6);
% phi_p=phi_s-1.57;

%% get angle from rms and rmp
c=rm1+1i*rm2;
ang=angle(c);   % change to MLE
ang2=[];
parfor ii=1:1:numel(rms)
    [ang2(ii,:)]=iPALM_est_angle(rms(ii),rmp(ii),phi_s,phi_p,ang(ii)-phi_s);
end

figure
plot(ang2(:,2))
hold on
plot(ang-phi_s,'r');
% figure
% plot(unwrap(ang2(:,2)),'--b');
% hold on
% plot(unwrap(ang-phi_s),'--r');

%%
[ang_coeff angp]=polyfit(xstep(1:N),unwrap(ang2(1:N,2)),1);
figure
plot(xstep(1:N),unwrap(ang2(1:N,2)),'b');
hold on
plot(xstep(1:N),unwrap(ang(1:N)-phi_s),'r');
plot(xstep(1:N),xstep(1:N)*ang_coeff(1)+ang_coeff(2),'--r');
% plot(ang-2*pi,'*b')

para.centers = centers;
para.file = [dataFolder imageName];

save([resultpath namestr '_dphi_cali_' datestring],'phi_s','phi_p','ang_coeff','para');
save([parentFolder namestr '_dphi_cali_' datestring],'phi_s','phi_p','ang_coeff','para');
set(handles.programStatus, 'String', {'Find Phase Shift', ['Color Channel:' channel],['Input:' dataFolder],'Calibration File used',fmtname,['Output:' parentFolder namestr 'dphi_cali' datestring '.mat']}); %show progress
drawnow update; %show progress
