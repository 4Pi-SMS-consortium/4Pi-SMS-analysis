function [q1,q2,q3,q4]=W4PiSMS_RotAlign_FMT(qd1,qd2,qd3,qd4,fmtname)
tmp=load(fmtname);
zm=tmp.zm_all;
trans=tmp.trans_all;
ang=tmp.ang_all;

q1=single(qd1);
q2=single(zeros(size(qd2)));
q3=single(zeros(size(qd3)));
q4=single(zeros(size(qd4)));

for ii=1:1:size(qd1,3)
    q2(:,:,ii)=single(affine_trans(qd2(:,:,ii),zm(2,:),trans(2,:),ang(2)));
    q3(:,:,ii)=single(affine_trans(qd3(:,:,ii),zm(3,:),trans(3,:),ang(3)));
    q4(:,:,ii)=single(affine_trans(qd4(:,:,ii),zm(4,:),trans(4,:),ang(4)));
    if mod(ii,1000)==0
        disp([num2str(ii) ' out of ' num2str(size(qd1,3)) ' is done...']);
    end
end