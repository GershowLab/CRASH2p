function zus = upsample_circular(z,us,dim)
%zus = upsample_circular(z,us,dim)

    if nargin < 3
        dim = ndims(z);
    end

    sz = size(z);    


    z = permute(z, circshift(1:ndims(z),(1-dim))); 
    args = cell(1,ndims(z));
    for j = 1:length(args)
        args{j} = ':';
    end
    t = 0:size(z,1);
    args{1} = 1;
    z = cat(1,z,z(args{:}));
    tus = (0:(us*t(end)-1))/us;
    zus = interp1(t,z,tus);
    zus = permute(zus,circshift(1:ndims(z),(dim-1)));
    szus = size(zus);
    sz(dim) = length(tus);
    if (~all(sz == szus))
        try
            zus = zus';
        catch
        end
    end

end