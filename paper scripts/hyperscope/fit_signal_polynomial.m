function [p,ffit] = fit_signal_polynomial(deg, taxis, signal, xdata, amplitude, f, weight)
%function [p,ffit] = fit_signal_polynomial(deg, taxis, signal, xdata, amplitude, f, weight)
%
% taxis time (or phase) axis, matching the size of signal
% xdata is size of f (e.g. f = X x Y x T)
% last dimension of f should match size of taxis, signal
% amplitude = amplitude of signal at each point in f
% deltat = polyval(pfit, xdata)
% minimizes weight*abs(f-signal(t - deltat))
% vel_x = axes of magnitudes of velocity (0 to vmax) in units of taxis / xaxis 
% theta_x = axes of directions (e.g. 0 to 2 pi) 

existsAndDefault('amplitude', ones(size(f)));
existsAndDefault('weight',  ones(size(f)));
existsAndDefault('taxis', 1:length(signal));

t = taxis;
T = max(t)-min(t);
t(end+1) = 2*t(end) - t(end-1);
tnorm = (t - t(1))/T;
F = griddedInterpolant(tnorm, signal([1:end 1]));
sz = size(f);
rr = sz; rr(1:(end-1)) = 1;
tt = repmat(reshape((taxis-taxis(1))/T, rr), sz(1:(end-1))); 


xr = max(xdata,[],'all','omitnan') - min(xdata,[],'all','omitnan');
xmid = 0.5 * (max(xdata,[],'all','omitnan') - min('xdata',[],'all','omitnan'));
xdatanorm = (xdata-xmid)/xr;

p = [];
fun = @(p) mean(weight.*abs(f-projectSignal(xdatanorm, tt, amplitude, 1, F, p)),'all','omitnan');

%find linear fit by exploring all parameter space first
[xx,yy] = meshgrid(linspace(-1,1,20), linspace(-1,1,20));
xx = xx(:);
yy = yy(:);
ff = zeros(size(xx));
for j = 1:numel(xx)
    ff(j) = fun([xx(j), yy(j)]);
end
p0 = [xx(argmin(ff)) yy(argmin(ff))];
lb = -ones(size(p0));
ub = ones(size(p0));
p = fmincon(fun, p0, [], [], [], [], lb, ub, [], optimset('display','off'));

%if deg > 1, then add one at a time 
for d = 2:deg
    q = fminbnd(@(q) fun([q p]), -1,1,optimset('display','off'));
    p0 = [q p];
    lb = -ones(size(p0));
    ub = ones(size(p0));
    p = fmincon(fun, p0, [], [], [], [], lb, ub, [], optimset('display','off'));
end

ffit = projectSignal(xdatanorm, tt, amplitude, 1, F, p);

%dt/T = p(1)*(x-xmid)/xr + p(2)
%dt = T*p(1)/xr * x + p(2)-p(1)*T/xr*xmid
p(2) = p(2)-p(1)*T/xr*xmid;
p(1) = p(1)*T/xr;



end

function ps = projectSignal(xdata, tt, amplitude, T, F, p)


    tt = tt - polyval(p,xdata);
    ps = amplitude.*F(mod(tt,T));
    

end
