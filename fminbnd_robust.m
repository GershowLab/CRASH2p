function [x,fval, exitflag] = fminbnd_robust(fun,x1, x2, ndivisions)
%function [x,fval, exitflag] = fminbnd_robust(fun,x1, x2, ndivisions)
%divides the interval into ndivisions subintervals (default 10)
%finds the coarse minimum value in each interval, then 
%then finds fine minimum value and location using location of minimum
%over all intervals
%
%attempts to avoid converging to local rather than global minimum

existsAndDefault('ndivisions', 10);

xx = linspace(x1, x2, ndivisions+1);
delta = diff(xx(1:2));
f = zeros([1 ndivisions]);
xm = f;
op = optimset('fminbnd');
op.TolX = delta/10;
op.Display = 'off';
for j = 1:ndivisions
    [xm(j), f(j)] = fminbnd(fun, xx(j), xx(j+1), op);
end
[~,I] = min(f);
[x,fval,exitflag] = fminbnd(fun, max(x1,xm(I)-delta), min(x2,xm(I)+delta));


end

