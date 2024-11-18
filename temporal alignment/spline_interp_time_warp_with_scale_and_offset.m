function [taligned, deltat, cpts, scale, offset] = spline_interp_time_warp_with_scale_and_offset (t_target, y_target, t, y, n_ctrl_pts, max_fractional_shift)
%function taligned = spline_interp_time_warp  (t_target, y_target, t, y, n_ctrl_pts, normalize=false, max_fractional_shift=1/8)
%
% t_target is 1 x T0
% y_target is 1 x T0
% t is 1 x T1
% y is 1 x T1
% 
% taligned is 1 x T1
% deltat is 1 x n_ctrl_pts
% cpts is 1 x n_ctrl_pts
%
% scale rescales y_target to match y
% offset adds to scale*y_target
%
% taligned = t + the b-spline interpolation of deltat at cpts
%   where cpts are n_ctrl_pts spread evenly from min(t) to max(t)
%   and deltat(1...n_ctrl_pts) minimizes sum (y - scale*interp(t_target,y_target,taligned)-offset)^2
%    subject to the constraint |deltat(i)| < (max(t)-min(t))*max_fractional_shift for all i
%
%
% it is suggested that t_target extend beyond t by at least max_fractional
%   shift in each direction
% otherwise nearest neighbor interpolation will be used for the endpoints


    

    existsAndDefault('max_fractional_shift', .125);

    valid = isfinite(t_target)&isfinite(y_target);
    t_target = t_target(valid);
    y_target = y_target(valid);

    f = griddedInterpolant(t_target, y_target, 'linear', 'nearest');
    df = griddedInterpolant(t_target, gradient(y_target)./gradient(t_target), 'linear', 'nearest');

    b = spline_component_polynomials(t, n_ctrl_pts);

    function [v,g] = obj_sqerr(x)
        tnew = t + x(3:end)*b;
        ftn = f(tnew);
        yf = ftn*x(1) + x(2);
        v = 0.5*sum((yf-y).^2);
        dtdx = b;
        dvdt = (yf-y).*df(tnew);
        g = 0*x;
        g(3:end) = sum(dtdx.*dvdt,2);
        g(1) = sum((yf-y).*ftn);
        g(2) = sum((yf-y));
    end

%     function [v,g] = obj_prod(x)
%         tnew = t + x*b;
%         yf = f(tnew);
%         v = -sum(yf.*y);
%         dtdx = b;
%         dvdt = -y.*df(tnew);
%         g = sum(dtdx.*dvdt,2);
%     end

    x0 = zeros(1,n_ctrl_pts+2);
    x0(1) = 1;

    obj_sqerr(x0);

    problem.objective = @obj_sqerr;
    problem.x0 = x0;
    problem.ub = max_fractional_shift*(max(t)-min(t))*ones(size(x0));
    problem.ub(1) = 1000;
    problem.ub(2) = max(abs(y))+max(abs(y_target)); 
    problem.lb = -problem.ub;
    problem.solver = 'fmincon';
    problem.options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','off');

    x = fmincon(problem);
    deltat = x(3:end);
    scale = x(1);
    offset = x(2);

    taligned = t + deltat*b;
    cpts = linspace(min(t),max(t),n_ctrl_pts);

    

end