function [estimates, model] = fit_ast(zdata, sigma)
% Call fminsearch with a random starting point.
start_point = [1.1 600 1200 0.1 0.1];
model = @astfun;
estimates = fminsearch(model, start_point);
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = astfun(params)
        w = params(1);
        c = params(2);
        d = params(3);
        A = params(4);
        B = params(5);
        FittedCurve = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
        ErrorVector = FittedCurve - sigma;
        sse = sum(ErrorVector .^ 2);
    end
end