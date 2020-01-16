function result = MyCalPahse_new(tint1,tint2,phase1,phase2)
%     tint1 = (int1a-int1b)./(int1a+int1b);
%     tint2 = (int2a-int2b)./(int2a+int2b);
    
    A=cos(phase1);
    B=sin(phase1);
    C=cos(phase2);
    D=sin(phase2);
    E=tint1;
    F=tint2;
    
    tcos = (D.*E-B.*F)./(A.*D-B.*C);
    tsin = (C.*E-A.*F)./(A.*D-B.*C);
        
    alpha = asincos(tsin, tcos);
    
    G=cos(alpha+phase1);
    H=cos(alpha+phase2);
    I=(E-F)./(G-H);
    result(:,1)=I;
    result(:,2)=alpha;    
end

function result = asincos(sind,cosd)
    temp = sind.^2+cosd.^2;
    sind = sind./sqrt(temp);
    cosd = cosd./sqrt(temp);
    
%     result = zeros(size(sind));
    result = acos(cosd);
    result(sind<0) = -result(sind<0);%+2*pi;
%     if(sind(1)>=0) %0 - pi
%         result = acos(cosd);
%     else%-pi - 0
%         result = -acos(cosd);
%     end
    
end