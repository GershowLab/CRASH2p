function [y,ty] = resample_to_t (tx, x, tr)
%function [y,ty] = resample_to_t (tx, x, tr)
%uses the matlab resample command to map x onto tr
%takes care of end points; tr must be evenly spaced
if isempty(tx) || isempty(x) || isempty(tr)
    y = [];
    ty = [];
    return;
end
tx = tx(:).';
tr = tr(:).';
trout = false;
if (size(x,2) ~= size(tx, 2))
    x = x.';
    trout = true;
end
if (size(x,2) ~= size(tx, 2))
    error ('length of x does not match length of tx');
end
dtr = mean(diff(tr));
txf = (tr(1)-3*dtr):median(diff(tx)):(tr(end)+3*dtr);
xf = interp1(tx, x, txf, 'linear', NaN);
xn = interp1(tx, x, txf, 'nearest', 'extrap');
xf(~isfinite(xf)) = xn(~isfinite(xf));
% 
% 
% xe = interp1(tx, x, tr([1 end]), 'nearest', 'extrap');
% ii = tx > tr(1) & tx < tr(end);
% txx = [tr(1) tx(ii) tr(end)];
% xx = [xe(:,1) x(:,ii) xe(:,end)];


[y,ty] = resample(xf, txf, 1/mean(diff(tr)));
% if (length(y) ~= length(tr))
    y = interp1(ty, y, tr, 'linear', 'extrap');
    ty = tr;
% end

if (trout)
    y = y.';
end


end

