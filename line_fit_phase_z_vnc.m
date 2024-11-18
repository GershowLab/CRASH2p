function [k,phi0,zfit] = line_fit_phase_z_vnc (x,y,z)
% minimizes sum |(z-|z| exp(i (kx x + ky y  + phi0))|^2
% maximizes (sum(Re |z| conj(z) * exp(i (kx x + ky y  + phi0)))
%p0 = [kx ky phi0]
% vnc specific assumptions: x is the long axis
% phase changes at most 2pi over x range
% phase changes at most pi/2 over y range

if (numel(x) ~= numel(z))
    [y,x] = meshgrid(x,y);
end

valid = isfinite(x)&isfinite(y)&isfinite(z);
x = x(valid);
y = y(valid);
z = z(valid);

z0 = mean(z);
z = z./z0;

prefac = abs(z).*conj(z);
x0 = mean(x,'all','omitnan');
y0 = mean(y,'all','omitnan');

x = x - x0;
y = y - y0;

xgrid = linspace(min(x),max(x),20);
ygrid = linspace(min(y),max(y),10);

[xa,mz_x] = meanyvsx(x,z,xgrid);
valid = isfinite(mz_x) & isfinite(xa);
xa = xa(valid);
mz_x = mz_x(valid); 
mz_x = lowpass1D(mz_x,2,'padType','circular');
phi = unwrap(angle(mz_x));

kx0 = median(deriv(phi,1))./median(deriv(xa,1));

[ya,mz_y] = meanyvsx(y,z,ygrid);
valid = isfinite(mz_y) & isfinite(ya);
ya = ya(valid);
mz_y = mz_y(valid); 
phi = unwrap(angle(mz_y));
ky0 = median(deriv(phi,1))./median(deriv(ya,1));

phi0 = angle(mean(z.*exp(-1i*(kx0*x + ky0*y))));
p0 = [kx0 ky0 phi0];


function [v,g,h] = objfun(p)
    phi = p(1)*x + p(2)*y + p(3);
    q = prefac.*exp(1i*phi);
    v = -sum(real(q),'all');
    g = [sum(imag(q).*x,'all'), sum(imag(q).*y,'all'), sum(imag(q),'all')];
    hxx = sum(real(q).*x.^2,'all');
    hyy = sum(real(q).*y.^2,'all');
    hxy = sum(real(q).*x.*y,'all');
    hx0 = sum(real(q).*x,'all');
    hy0 = sum(real(q).*y,'all');
    h00 = -v;
    h = [hxx hxy hx0; hxy hyy hy0; hx0 hy0 h00];
end

dx = max(x)-min(x);
dy = max(y)-min(y);
ub = [2*pi/dx 0.5*pi/dy pi];
lb = -ub;
op = optimoptions("fmincon","SpecifyObjectiveGradient",true,"MaxIterations",1e4,"Display","off");

p0(p0 >= ub) = 0.95*ub(p0 >= ub);
p0(p0 <= lb) = 0.95*lb(p0 <= lb);


objfun(p0);
p = fmincon(@objfun, p0, [],[],[],[],lb,ub,[],op);

k = p(1:2);
phi0 = p(3) - k(1)*x0 - k(2)*y0 + angle(z0);
phifit = p(1)*x + p(2)*y + p(3) + angle(z0);
zfit = abs(z)*abs(z0).*exp(1i*phifit);


end
