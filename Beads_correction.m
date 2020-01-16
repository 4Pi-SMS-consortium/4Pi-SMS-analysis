function [Vout, drift, flag]=Beads_correction(Vin,reverseflag)
V=Vin;
m=sum((V(:,3)==1));
n=0;
N=max(Vin(:,3))-min(Vin(:,3))+1;
for i=1:m
    id=abs(V(:,1)-V(i,1))<2 & abs(V(:,2)-V(i,2))<2;
    V1=V(id,:);
    V1(:,5)=find(id>0);    
    U=unique(V1(:,3));
    if length(U)>=N-1000                 %% exist for every frame 
        n=n+1;
        V(i,4)=n;
        V(id,4)=n;  
        if length(V1(:,1))>length(U)%size(V1,1)>=N-3000 && size(V1,1)<=N+3000      %% more than frame number
            %% find frame contains more than one molecule
            M=hist(V1(:,3),U);
            FI=find(M>1);           
            for k=1:length(FI)
                f=FI(k);
                id=V1(:,3)==U(f);
                Vnew=V1(id,:);
                id=V1(:,3)==U(f)-1;
                cx=mean(V1(id,1));
                cy=mean(V1(id,2));
                if isempty(cx)
                    continue
                end
                dist=(Vnew(:,1)-cx).^2+(Vnew(:,2)-cy).^2;
                FI1=find(dist==min(dist));             %% keep the closest molecule
                for j=1:size(Vnew,1)
                    if j~=FI1
                        n1=Vnew(j,5);
                        V(n1,4)=0;
                    end
                end
            end
        end
    end
end

m=max(V(:,4));
drift=zeros(N,4);
Vout=V(:,1:2);
Vout(:,3:4)=V(:,6:7);
n=0;
for i=1:m
    IX=V(:,4)==i;
    I=V(IX,5);
    if sum(IX)<N
        Vtemp=NaN(N,7);
        Vnew=V(IX,:); 
        Vtemp(Vnew(:,3),:)=Vnew;
        Vtemp=fillnan(Vtemp);
        drift(:,1:2)=drift(:,1:2)+Vtemp(:,1:2);
        drift(:,3)=drift(:,3)+unwrap(Vtemp(:,6));
%         drift(:,3)=drift(:,3)+Vtemp(:,6)-Vtemp(1,6);
%         drift(:,3)=drift(:,3)+wrapToPi(Vtemp(:,6)-Vtemp(1,6)+pi-0.2);
        drift(:,4)=drift(:,4)+Vtemp(:,7);
        n=n+1;
    end
    if sum(IX)==N
        Vnew=V(IX,:);
        drift(:,1:2)=drift(:,1:2)+Vnew(:,1:2);
%         drift(:,3)=drift(:,3)+Vnew(:,6)-Vnew(1,6);
        drift(:,3)=drift(:,3)+unwrap(Vnew(:,6))-Vnew(1,6);
%         drift(:,3)=drift(:,3)+wrapToPi(Vnew(:,6)-Vnew(1,6))+2;
        drift(:,4)=drift(:,4)+Vnew(:,7);
        n=n+1;
    end     
end

disp([num2str(n),' beads are used for drift correction.']);
if n>0
    drift(:,1)=smooth(drift(:,1),100);
    drift(:,2)=smooth(drift(:,2),100);
    drift(:,3)=smooth(drift(:,3),100);
    drift(:,4)=smooth(drift(:,4),100);
    if reverseflag
        drift=drift-ones(N,1)*drift(end,:);
    else
        drift=drift-ones(N,1)*drift(1,:);
    end
    drift=drift/n;
    flag=1;
    D=drift(V(:,3),1:4);
    Vout=Vout-D;
else
    flag=0;
end

% disp([num2str(n),' beads are used for drift correction.']);
% if n==0
%     load('G:\4PISCMOS\2016-12-22\Cell07_drift.mat');
%     flag=1;
%     D=drift(V(:,3),1:4);
%     Vout=Vout-D;
% else
%     flag=0;
% end