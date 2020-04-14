function [flag]=W4PiSMS2vutarav2(resultfolder,filename,colornum,xcell,ycell,zcell,tcell,Icell,crlbcell,LLRcell,conf1,conf2,conf3)

Mtot=[];
flag=0;

rfolder=[resultfolder filename '\'];
mkdir(rfolder);

c =['image-ID,cycle,z-step,frame,accum,probe,photon-count,photon-count11,photon-count12,photon-count21,photon-count22,psfx,psfy,psfz,psf-photon-count,x,y,z,stdev,amp,background11,background12,background21,background22,chisq,log-likelihood,custom-confidence1,custom-confidence2,custom-confidence3,valid'];
fid = fopen([rfolder 'particles.csv'],'wt');
fprintf(fid,'%s\n',c);
fclose(fid);
for cc=1:1:colornum
    Mtottmp=[];
    x=xcell{cc};
    y=ycell{cc};
    z=zcell{cc};
    col1=tcell{cc};
    col2=ones(size(x));
    col3=ones(size(x));
    col4=tcell{cc};
    col5=ones(size(x));
    col6=ones(size(x)).*0+cc-1;
    col7=Icell{cc};
    col8=ones(size(x));
    col9=ones(size(x));
    col10=ones(size(x));
    col11=ones(size(x));
    col12=ones(size(x));
    col13=ones(size(x));
    col14=ones(size(x));
    col15=ones(size(x));   
    col16=x;
    col17=y;
    col18=z;
    col19=ones(size(x));
    col20=ones(size(x));
    col21=ones(size(x));
    col22=ones(size(x));
    col23=ones(size(x));
    col24=ones(size(x));
    col25=crlbcell{cc};
    col26=LLRcell{cc};
    col27=conf1{cc};
    col28=conf2{cc};
    col29=conf3{cc};
    col30=ones(size(x));
    kk=30;
    
    %% cat them together
    for ii=1:1:kk
        eval(['Mtottmp=cat(2,Mtottmp,' ['col' num2str(ii)] ');']);
    end
    
    if cc==1
        Mtot=Mtottmp;
    else
        Mtot=cat(1,Mtot,Mtottmp);
    end
end

try
    dlmwrite([rfolder 'particles.csv'],Mtot,'-append')
    flag=1;
catch
    error('Error in writing csv file, please check write permission in current folder');
end