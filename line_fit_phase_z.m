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

sz = size(z);

prefac = abs(z).*conj(z);
x0 = mean(x,'all','omitnan');
y0 = mean(y,'all','omitnan');

x = x(:) - x0;
y = y(:) - y0;
z = z(:);
xgrid = linspace(min(x),max(x),50);
ygrid = linspace(min(y),max(y),25);



valid = isfinite(x)&isfinite(y)&isfinite(z);


% 
% 
% if (nargin < 4)
%     z2 = imfilter(imfilter(z,gaussKernel(4),'same','circular'),gaussKernel(4)','same','circular');
%     phi2 = unwrap(angle(z2),[],1); phi2 = phi2 - mean(phi2); 
%     kx = sum(x.*phi2,'all')./sum(x.*x,'all');
% %     plot(mean(x,2), mean(phi2,2), mean(x,2), kx*mean(x,2))
% 
%     phi2 = unwrap(angle(z2),[],1); phi2 = phi2 - mean(phi2); 
%     ky = sum(y.*phi2,'all')./sum(y.*y,'all');
% 
% 
%     p0 = [kx ky angle(mean(z2.*exp(-1i*kx*x - 1i*ky*y),'all'))];
%     
%     [k,phi0,zfit0] = line_fit_phase_z(x,y,z2,p0);
%     p0 = [k phi0];
% %     subplot(1,2,1)
% %     pcolor(y,x,angle(z2)); shading flat; axis equal
% %     subplot(1,2,2)
% %     pcolor(y,x,angle(zfit0)); shading flat; axis equal
% end



valid = isfinite(x)&isfinite(y)&isfinite(z);
x = x(valid);
y = y(valid);
z = z(valid);
prefac = prefac(valid);

if (nargin < 4)
    p0 = [0 0 mod(angle(mean(z)) + pi, 2*pi)-pi];
end
function [v,g] = initfun(p)

    [v,g] = objfun([p(1) 0 p(2)]);
    g = g([1 3]);


end

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
    lscale = max(max(x,[],'all')-min(x,[],'all'), max(y,[],'all')-min(y,[],'all'));%2*(sqrt(max(x.^2,'all') + max(y.^2,'all')));
    op = optimoptions("fmincon","SpecifyObjectiveGradient",true,"MaxIterations",1e4,"Display","off");
    %ub = [2*pi./(max(x,[],'all')-min(x,[],'all')) 2*pi./(max(y,[],'all')-min(y,[],'all')) pi];

%         objfun(p0)

%     p0 = fmincon(@initfun, [p0(1) p0(3)], [],[],[],[],[-2*pi/lscale -pi],[2*pi/lscale pi],[],op);
% 
%     p0 = [p0(1) 0 p0(2)];

    ub = [8*pi/lscale 8*pi/lscale 2*pi];

    

    lb = -ub;
%     objfun(p0)
    p = fmincon(@objfun, p0, [],[],[],[],lb,ub,[],op);
%     objfun(p)
    k = p(1:2);
    phi0 = p(3) - k(1)*x0 - k(2)*y0;
    phifit = p(1)*x + p(2)*y + p(3);
    zfit = reshape(abs(z).*exp(1i*phifit),sz);


end
