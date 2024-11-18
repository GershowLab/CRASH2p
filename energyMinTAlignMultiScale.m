function [tfit,mag] = energyMinTAlignMultiScale (txknown, vknown, vunknown, tends, mindelta,  alpha_spring)
%function [tfit] = energyMinTAlignMultiScale (txknown, vknown, vunknown, tends, mindelta,  alpha_spring)
%sigma for derivatives
%existsAndDefault('sigma', 0.01);


if (size(vknown,1) < size(vknown,2))
    vknown = vknown';
end
if (size(vunknown,1) < size(vunknown,2))
    vunknown = vunknown';
end

initFitLength = 128;
if (size(vunknown,1) > initFitLength)
    scale = 2;
    % = downsamp(txknown, scale, 'linear');
    [vknown_res, txknown_res] = downsamp(vknown, scale, 'ends', txknown);
    [vuknown_res, inds] = downsamp(vunknown, scale, 'ends', (1:size(vunknown, 1))');
%    inds = downsamp(1:size(vunknown, 1), scale, 'linear');
    
    [tfit,mag] = energyMinTAlignMultiScale(txknown_res, vknown_res, vuknown_res, tends, mindelta, alpha_spring);
    tinitial = interp1([0; inds; (size(vunknown, 1)+1)], [tends(1); tfit; tends(2)], 1:size(vunknown, 1), 'linear');
    maginitial = mag;
else
    scale = 1;
    tinitial = linspace(tends(1), tends(2), size(vunknown,1)+2);
    tinitial = tinitial(2:(end-1));
    maginitial = [];
end

[tfit,mag] = energyMinTAlign(txknown, vknown, vunknown, tends, mindelta*scale, alpha_spring, tinitial, maginitial);
figure();
if (isempty(maginitial))
    maginitial = 1;
end
plot (txknown, vknown(:,1), tinitial, vunknown(:,1)/maginitial, tfit, vunknown(:,1)/mag);

function [xdown, t] = downsamp(x, scale, padType, t)
    %t is not lowpassed
    if (size(x,1) < size(x,2))
        x = x';
    end
    x = lowpass1D(x, scale/2, 'padType', padType);
    ii = 1:size(x,1);
    xdown = interp1(ii, x, downsample(ii,scale, round(scale/2)));
    if (nargin > 3 && nargout > 1)
        t = downsample(t, scale, round(scale/2));
    end
end
    
    
end
