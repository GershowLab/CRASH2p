function im = imblur(im, sigma)
%function im = imblur(im, sigma)
%im can be any dimension
%sigma is same ordering as im, if sigma is 1 element same used for all


for j = 1:ndims(im)
    sz = num2cell(ones([ndims(im) 1]));
    sz{j} = [];
    gk = reshape(gaussKernel(sigma(min(j,length(sigma)))), sz{:});
    im = imfilter(im, gk, 'same', 'replicate');
end
