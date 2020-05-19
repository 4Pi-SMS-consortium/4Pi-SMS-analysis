function W4PiSMS_getPhaseShift(handles)

close all
center1 = str2num(get(handles.center1,'string'));
center2 = str2num(get(handles.center2,'string'));
center3 = str2num(get(handles.center3,'string'));
center4 = str2num(get(handles.center4,'string'));
centers = [center1 center2 center3 center4];

posfile = get(handles.pathMainfolder, 'String');
c = posfile(end); %get last character of path to image file
while c ~= '\' %get path only
    posfile(end) = []; %delete filename character
    if isempty(posfile) %break if no path present
        break;
    end
    c = posfile(end); %get last charakter of path to image file
end

[imageName,dataFolder] = uigetfile([posfile '*.dcimg'],'Open Calibration Images','MultiSelect','on');
I=find(dataFolder=='\',2,'last');
parentFolder = dataFolder(1:I);

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
resultpath = [pwd '\calibration_files\'];
namestr= ['bead_' channel];
mkdir(resultpath);

%%
load(handles.scmos_cali_file);
caliims=cat(3,offsetim,varim,gainim);

aveqd1=[];
aveqd2=[];
aveqd3=[];
aveqd4=[];
set(handles.programStatus,'string','Read Beads Images')
drawnow update
for ii=1:1:numel(imageName)
    [~,qds,calicrops]=W4PiSMS_readdcimg([dataFolder imageName{ii}],centers,caliims);
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
end

%% rotate and align
[q1 q2 q3 q4]=W4PiSMS_RotAlign_FMT(aveqd1,aveqd2,aveqd3,aveqd4,fmtname);
colorim1=joinchannels('RGB',q3,q1);
colorim2=joinchannels('RGB',q4,q2);

%% find local centers
for kk=1:1
    sumim1=q1+q4+q2+q3;
    sumim1(sumim1<=1e-6)=1e-6;
    subsz=21;
    cutim=[];
    cut1=[];
    cut2=[];
    cut3=[];
    cut4=[];
    N=size(q1,3);
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
    
    for ii=1:1:size(sumim1,3)
        tmp=sumim1(:,:,ii);
        cutim(:,:,ii)=cut(tmp,[subsz subsz],coords-10);
        cut1(:,:,ii)=cut(q1(:,:,ii),[subsz subsz],coords-10);
        cut2(:,:,ii)=cut(q2(:,:,ii),[subsz subsz],coords-10);
        cut3(:,:,ii)=cut(q3(:,:,ii),[subsz subsz],coords-10);
        cut4(:,:,ii)=cut(q4(:,:,ii),[subsz subsz],coords-10);
    end
    set(handles.programStatus,'string','Fit the Beads')
    drawnow update
    
    % fit sigmax sigmay
    set(handles.programStatus,'string','Plot and fit Sigmax and SigmaY')
    drawnow update
    [P,CRLB,LL]=mleFit_LM(single(cutim),4,50,1.4,0,0,0);    
    zstep=(1:1:length(P(:,1)))'*20;
    xco=P(:,2);
    yco=P(:,1);
    sigmax=P(:,5);
    sigmay=P(:,6);
    S=abs(sigmax-sigmay);
    id=find(S==min(S));
    disp(['Focus position (step number): ',num2str(id)]);
end

%%
figure
[estx,model]=fit_ast(zstep,sigmax);
[~,fit_curve]=model(estx);
plot(zstep,sigmax,'bo');
hold on
plot(zstep,fit_curve,'b--','linewidth',4);

[esty,model]=fit_ast(zstep,sigmay);
[~,fit_curve]=model(esty);
plot(zstep,sigmay,'ro');
plot(zstep,fit_curve,'r--','linewidth',4);

xlabel('relative z (nm)');
ylabel('sigma');
legend('Obs. \sigma_x','Fit. \sigma_x','Obs. \sigma_x','Fit. \sigma_y','Location','North');
xlim([0 1200]);

%% save
datestring=datestr(now,'yyyymmdd');
save([resultpath namestr '_Astfit_' datestring],'model','estx','zstep','esty','sigmax','sigmay');
save([parentFolder namestr '_Astfit_' datestring],'model','estx','zstep','esty','sigmax','sigmay');
set(handles.programStatus,'string','Find phase shift')
drawnow update

%% get rms rmp
xf=xco;
yf=yco;
subims=cat(4,cut1,cut2,cut3,cut4);
[rms,rmp]=W4PiSMS_findmom_givenC(subims,xf,yf);

%% calculate reduced moment
rm1=rms./moment_normalization(rms);
rm2=rmp./moment_normalization(rmp);

figure;
plot(rm1,'r');
hold on;
plot(rm2,'b');
plot(-rm1,'--r');
plot(-rm2,'--b');
legend('rm1','rm2','inv_rm1','inv_rm2');

%% estimate the phase delay
A1=0.9;
w=0.43;
phi1=1.4;
b1=0.1;
A2=0.9;
phi2=phi1-1.57;
b2=b1;

rmsin=rm1(id-20:id+20);
rmpin=rm2(id-20:id+20);
x=1:1:numel(rmsin);

% fit indivual curve first
[est1 mode1l]=W4PiSMS_find_delay_single(x,rmsin,[A1 w phi1 b1]);
[est2 model2]=W4PiSMS_find_delay_single(x,rmpin,[A2 w phi2 b2]);

figure
plot(x,rmsin,'r');
hold on
plot(x,rmpin,'b');
start_point=[est1(1,1) est1(1,2) est1(1,3) est1(1,4) est2(1,1) est2(1,3) est2(1,4)];
[est model]=W4PiSMS_find_delay(x,rmsin,rmpin,start_point);

[~,fitrms,fitrmp]=feval(model,est);

plot(x,fitrms,'--r');
hold on
plot(x,fitrmp,'--b');

phi_s=est(3);
phi_p=est(6);
dphi=phi_s-phi_p;
disp(['Phase shift (rad): ',num2str(dphi)]);

%% get angle from rms and rmp
c=rm1+1i*rm2;
ang=angle(c);   % change to MLE
ang2=[];
for ii=1:1:numel(rms)
    [ang2(ii,:)]=W4PiSMS_est_angle(rms(ii),rmp(ii),phi_s,phi_p,ang(ii)-phi_s);
end

figure
plot(ang2(:,2))
hold on
plot(ang-phi_s,'r');

%%
[ang_coeff,angp]=polyfit(zstep(1:N),unwrap(ang2(1:N,2)),1);

figure
plot(zstep(1:N),unwrap(ang2(1:N,2)),'b');
hold on
plot(zstep(1:N),unwrap(ang(1:N)-phi_s),'r');
plot(zstep(1:N),zstep(1:N)*ang_coeff(1)+ang_coeff(2),'--r');

para.centers = centers;
para.file = [dataFolder imageName];
save([resultpath namestr '_dphi_cali_' datestring],'phi_s','phi_p','ang_coeff','para');
save([parentFolder namestr '_dphi_cali_' datestring],'phi_s','phi_p','ang_coeff','para');
set(handles.programStatus, 'String', {'Find Phase Shift', ['Color Channel:' channel],['Input:' dataFolder],'Calibration File used',fmtname,['Output:' parentFolder namestr 'dphi_cali' datestring '.mat']}); %show progress
drawnow update; %show progress