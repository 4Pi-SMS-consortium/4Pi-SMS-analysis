function [x]=W4PiSMS_est_angle(rms,rmp,phi1,phi2,ang0)
a0=0.5;
options = optimset('Display','off');
x=fsolve(@fun,[a0 ang0],options);
    function f=fun(x)
        f(1)=abs(x(1))*cos(x(2)+phi1)-rms;
        f(2)=abs(x(1))*cos(x(2)+phi2)-rmp;
    end
end