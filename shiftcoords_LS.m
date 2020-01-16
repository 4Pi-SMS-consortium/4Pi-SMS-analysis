% for v14 iPALMast_analysisv14scmos.m

function [out]=shiftcoords_LS(xout,shiftx,toutin,frmnum,reverseflag)

if nargin<=4
    segnum=ceil((max(toutin)+1)/frmnum);
    for ii=2:1:segnum
        st1=(ii-1)*frmnum;
        ed1=(ii)*frmnum-1;
        maskt=toutin>=st1&toutin<=ed1;
        xout(maskt)=xout(maskt)-sum(shiftx(1:ii-1,1));
    end
    out=xout;
elseif reverseflag==0
    segnum=ceil((max(toutin)+1)/frmnum);
    for ii=2:1:segnum
        st1=(ii-1)*frmnum;
        ed1=(ii)*frmnum-1;
        maskt=toutin>=st1&toutin<=ed1;
        xout(maskt)=xout(maskt)-sum(shiftx(1:ii-1,1));
    end
    out=xout;
else
    segnum=ceil((max(toutin)+1)/frmnum);
    for ii=1:1:segnum-1
        st1=(ii-1)*frmnum;
        ed1=(ii)*frmnum-1;
        maskt=toutin>=st1&toutin<=ed1;
        xout(maskt)=xout(maskt)+sum(shiftx(ii:segnum-1,1));
    end
    out=xout;
end


    