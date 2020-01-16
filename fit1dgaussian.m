function [c]=fit1dgaussian(y,x,inis)

start_point = [inis min(y(:))];
model = @Gaussian1D;
c = fminsearch(model, start_point);
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = Gaussian1D(params)
        A=params(1);
        center=params(2);
        sig=params(3);
        bg=params(4);
        
        FittedCurve = A .* exp(-(x-center).^2./2./sig./sig)+bg;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end
end