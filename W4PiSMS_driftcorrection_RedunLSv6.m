function [xout2,yout2,zout2,shifts]=W4PiSMS_driftcorrection_RedunLSv6(xout,yout,zout,tout,frmnum,reverseflag,interpflag,str)

cutmeth='nocut';
pixelsz=16; % nm
thresh=8; % nm

if exist(str,'file')
    load(str);
else   
    [shiftx,shifty]=correct_drift_LS(single(xout),single(yout),tout,frmnum,pixelsz,thresh,cutmeth);
    if sum(shiftx.^2)<sum(shifty.^2)
        [shiftx2,shiftz]=correct_drift_LS(single(xout),single(zout),tout,frmnum,pixelsz,thresh,cutmeth);
    else
        [shifty2,shiftz]=correct_drift_LS(single(yout),single(zout),tout,frmnum,pixelsz,thresh,cutmeth);
    end
    save(str,'shiftx','shifty','shiftz');
end

if ~interpflag
    [xout2]=shiftcoords_LS(xout,shiftx,tout,frmnum,reverseflag);
    [yout2]=shiftcoords_LS(yout,shifty,tout,frmnum,reverseflag);
    [zout2]=shiftcoords_LS(zout,shiftz,tout,frmnum,reverseflag);
    shifts=[cumsum(shiftx(:)),cumsum(shifty(:)),cumsum(shiftz(:))];
else
    ntotalframe=max(tout)+1;
    nbinframe=length(shiftx)+1;
    indexinterp=zeros(nbinframe+2,1);
    indexinterp(1)=1;
    indexinterp(nbinframe+2)=ntotalframe;
    indexinterp(2:nbinframe+1)=round(frmnum/2):frmnum:frmnum*nbinframe-1;
    drift=cumsum(shiftx);
    finaldrift(:,1)=interp1(indexinterp,[0 0 drift' drift(end,1)],1:ntotalframe,'linear')';
    drift=cumsum(shifty);
    finaldrift(:,2)=interp1(indexinterp,[0 0 drift' drift(end,1)],1:ntotalframe,'linear')';
    drift=cumsum(shiftz);
    finaldrift(:,3)=interp1(indexinterp,[0 0 drift' drift(end,1)],1:ntotalframe,'linear')';
    if reverseflag
        finaldrift=finaldrift-ones(length(finaldrift),1)*finaldrift(end,:);
    end
    shift=finaldrift(tout+1,1);
    xout2=xout-shift;
    shift=finaldrift(tout+1,2);
    yout2=yout-shift;
    shift=finaldrift(tout+1,3);
    zout2=zout-shift;
    shifts=finaldrift;
end