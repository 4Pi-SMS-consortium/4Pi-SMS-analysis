function [est model]=W4PiSMS_find_delay_single(x,rm,start_point)

model = @phifun;
est = fminsearch(model, start_point);
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, fitrm] = phifun(params)
        A1 = params(1);
        w=params(2);
        phi1=params(3);
        b1=params(4);
        fitrm=abs(A1).*cos(w.*x+phi1)+b1;
%         FittedCurve = A .* exp(-lambda * xdata);
        ErrorVector = (fitrm(:) - rm).^2;
        sse = sum(ErrorVector(:));
    end
end