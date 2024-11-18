function [x,yfit] = fminArbCostFunction (fun,x0, xdata, ydata, costfun, op)
% function fminArbCostFunction (f, xdata, ydata, costfun, op)
%
% minimizes sum_xdata (costfun(ydata - f(x0,xdata))
% if (op.GradObj is 'on') then both costfun and f should return a gradient
% x0 is 1 x D
% xdata, ydata are N x 1
% fun(xdata) is N x 1
% costfun has a N x 1 derivative, f has a N x D  derivative 

    existsAndDefault('op', optimoptions('fminunc'));
    hasderiv = strcmpi(op.GradObj, 'on');
    if (size(xdata,1) == 1)
        xdata = xdata.';
    end
    if (size(ydata, 1) == 1)
        ydata = ydata.';
    end
    
    x = fminunc(@objfun, x0, op);
    yfit = fun(x, xdata);
    

    function [v,g] = objfun(x)
        if(hasderiv)
            [f,dfdx] = fun(x,xdata);
            [vv,dg] = costfun(ydata - f);
            v = sum(vv);
            g = -sum(repmat(dg, [1 size(dfdx, 2)]).*dfdx);
        else
            v = sum(costfun(ydata - fun(x,xdata)));
        end
    end
end