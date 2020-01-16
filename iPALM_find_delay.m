function [est model]=iPALM_find_delay(x,rms,rmp,start_point)
% assuming rm1=cos(theta)
% sin(theta)=sqrt(1-rm1.^2);
% startpoint=[-1.5];
model = @phifun;
est = fminsearch(model, start_point);
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, fitrms,fitrmp] = phifun(params)
        A1 = params(1);
        w=params(2);
        phi1=params(3);
        b1=params(4);
        A2=params(5);
        phi2=params(6);
        b2=params(7);
        fitrms=abs(A1).*cos(w.*x+phi1)+b1;
        fitrmp=abs(A2).*cos(w.*x+phi2)+b2;
%         FittedCurve = A .* exp(-lambda * xdata);
        ErrorVector = (fitrms(:) - rms).^2+(fitrmp(:) - rmp).^2;
        sse = sum(ErrorVector(:));
    end
end