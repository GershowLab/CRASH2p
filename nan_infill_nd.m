function x = nan_infill_nd(x)


for j = 1:max(size(x)) %should terminate before reaching this point, but prevents infinite loop in case all of x is nan
    nhood = conndef(ndims(x), 'minimal');
    valid = isfinite(x);
    x(~valid) = 0;
    num = convn(x, nhood, 'same');
    den = convn(double(valid), nhood, 'same');
    x(~valid) = num(~valid)./den(~valid);
    if(all(isfinite(x),'all'))
        break;
    end

end
