classdef FastAssembler
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xedges;
        yedges;
        zedges;
        tedges;
        dxdz = 0;
        dydz = 0; 
        maxDwellExpansion = 10;
        zOfP = @(p) cos(p);
        depthCorrection = @(p) abs(sin(p))/0.7846; % .7846 instead of 2/pi = .6366 reflects normalization between 20&160 over 22 bins
        
        %ordering is im(x,y,z) NOT im(row,col,z)
    end
    
    properties (Access = public)
        x0 %x0,y0,z0 is center of origin. 
        y0
        z0
        t0
        
        nx
        ny
        nz
        nt
        
        dx
        dy
        dz
        dt
        
    end
    
    methods
        
        function fa = FastAssembler(xedges,yedges,zedges,tedges, makeEdgesFromCenters, z_scale)
             if (nargin < 5 || isempty(makeEdgesFromCenters))
                 makeEdgesFromCenters = false;
             end
             if (nargin >= 6 && ~isempty(z_scale))
                 fa.zOfP = @(p) z_scale*cos(p);
             end
             [fa.xedges, fa.x0, fa.dx, fa.nx] = FastAssembler.regularizeGrid(xedges, makeEdgesFromCenters);
             [fa.yedges, fa.y0, fa.dy, fa.ny] = FastAssembler.regularizeGrid(yedges,makeEdgesFromCenters);
             if (nargin > 2 && length(zedges) > 1)
                 [fa.zedges, fa.z0, fa.dz, fa.nz] = FastAssembler.regularizeGrid(zedges,makeEdgesFromCenters);
             end  
             if (nargin > 3 && length(tedges) > 1)
                 [fa.tedges, fa.t0, fa.dt, fa.nt] = FastAssembler.regularizeGrid(tedges,makeEdgesFromCenters);
             end  
        end
        
        function [xc,yc,zc,tc] = getBinCenters(fa)
            xc = [];
            yc = [];
            zc = [];
            tc = [];
            if (numel(fa.xedges) > 1)
                xc = fa.xedges(1:end-1) + fa.dx/2;
               % xc = 0.5*(fa.xedges(1:(end-1)) + fa.xedges(2:end));
            end
            if (numel(fa.yedges) > 1)
                yc = fa.yedges(1:end-1) + fa.dy/2;
                %yc = 0.5*(fa.yedges(1:(end-1)) + fa.yedges(2:end));
            end
            if (numel(fa.zedges) > 1)
                zc = fa.zedges(1:end-1) + fa.dz/2;
                %zc = 0.5*(fa.zedges(1:(end-1)) + fa.zedges(2:end));
            end
            if (numel(fa.tedges) > 1)
                tc = fa.tedges(1:end-1) + fa.dt/2;
              %  tc = 0.5*(fa.tedges(1:(end-1)) + fa.tedges(2:end));
            end
            
        end
        
        function [bin, subz] = accumVals(fa, val, x, y, z)
            %The value in each row of the m-by-n matrix, subs, specifies an n-dimensional index into the output, A. For example, if subs is a 3-by-2 matrix, it contains three 2-D subscripts. subs also can be a column vector of indices, in which case the output, A, is also a column vector.
            sz = fa.getImSize();
            threed = length(sz) > 2 && sz(3) > 1 && nargin > 4;
            [xi, inrange] = FastAssembler.valToInd(x, fa.x0, fa.dx, fa.nx);
            [yi, iry] = FastAssembler.valToInd(y, fa.y0, fa.dy, fa.ny);
            inrange = inrange & iry;
            if (threed)
                [zi, irz] = FastAssembler.valToInd(z, fa.z0, fa.dz, fa.nz);
                inrange = inrange & irz;
            end
            xi = xi(inrange);
            yi = yi(inrange);
            val = val(inrange);
            if (threed)
                zi = zi(inrange);
                sz = sz(1:3);
                subs = [xi(:) yi(:) zi(:)];
            else
                sz = sz(1:2);
                subs = [xi(:) yi(:)];
            end
            bin = accumarray(subs, val, sz, @sum, 0);
            if (nargout > 1)
                subz = NaN(length(x), size(subs,2));
                subz(inrange,:) = subs;
            end
        end
        
        function [bin,subz] = accumVals4(fa, val, x, y, z, t)
            %The value in each row of the m-by-n matrix, subs, specifies an n-dimensional index into the output, A. For example, if subs is a 3-by-2 matrix, it contains three 2-D subscripts. subs also can be a column vector of indices, in which case the output, A, is also a column vector.
            sz = fa.getImSize();
            [xi, inrange] = FastAssembler.valToInd(x, fa.x0, fa.dx, fa.nx);
            [yi, iry] = FastAssembler.valToInd(y, fa.y0, fa.dy, fa.ny);
            inrange = inrange & iry;
            [zi, irz] = FastAssembler.valToInd(z, fa.z0, fa.dz, fa.nz);
            inrange = inrange & irz;
            [ti, irt] =  FastAssembler.valToInd(t, fa.t0, fa.dt, fa.nt);
            inrange = inrange & irt;
            xi = xi(inrange);
            yi = yi(inrange);
            val = val(inrange);
            zi = zi(inrange);
            ti = ti(inrange);
            subs = [xi(:) yi(:) zi(:) ti(:)];
            bin = accumarray(subs, val, sz, @sum, 0);
            
            if (nargout > 1)
                subz = NaN(length(x), size(subs,2));
                subz(inrange,:) = subs;
            end
        end
        
        function [bin,subs] = binWithWeightAndAffine (fa, x, y, p, weight, at)
            x = x(:);
            y = y(:);
            p = p(:);           
            z = fa.zOfP(p);
            x = x + fa.dxdz*z;
            y = y + fa.dydz*z;
            if (nargin < 5 || isempty(weight))
                weight = ones(size(x));
            end
            
            if (nargin >= 6 && ~isempty(at))
                ploc = at*[x';y';z';ones(size(x'))];
                x = ploc(1,:)';
                y = ploc(2,:)';
                z = ploc(3,:)';            
            end
            [bin,subs] = fa.accumVals(weight, x,y, z);  
            
        end
            
        function [bin,subs] = binPhotons (fa, x, y, p, t, at)  
            %        function [bin,subs] = binPhotons (fa, x, y, p, t, at)  

            x = x(:);
            y = y(:);
            p = p(:);           
            z = fa.zOfP(p);
            weight = fa.depthCorrection(p);
            x = x + fa.dxdz*z;
            y = y + fa.dydz*z;
            if (nargin >= 6 && ~isempty(at))
                ploc = at*[x';y';z';ones(size(x'))];
                x = ploc(1,:)';
                y = ploc(2,:)';
                z = ploc(3,:)';
                
            end
            
            if (nargin < 5 || isempty(t))
                [bin,subs] = fa.accumVals(weight, x,y, z);               
            else
                t = t(:);
                [bin,subs] = fa.accumVals4(weight, x, y, z, t);                
            end
%             if (nargin >= 6 && ~isempty(at))
%                 for j = 1:size(bin,4)
%                     bin(:,:,:,j) = fa.applyAffine(bin(:,:,:,j), at(:,:,min(j,size(at,3))));
%                 end
%             end
        end
        
        
        
        
        function bin = binDwell (fa, x, y, tau, t, at)
            %function bin = binDwell (fa, x, y, tau, t, at)
            %x,y are the locations when z = 0
            if (nargin < 4 || isempty(tau))
                tau = ones(size(x));
            end
            if (length(tau) == 1)
                tau = ones(size(x))*tau;
            end
            
            
            
            if (nargin > 4 && ~isempty(t))
                [ti, inrange] =  FastAssembler.valToInd(t, fa.t0, fa.dt, fa.nt);
                ti = ti(inrange);
                x = x(inrange);
                y = y(inrange);
                tau = tau(inrange);
                bin = zeros(fa.nx, fa.ny, fa.nz, fa.nt);
                for j = 1:fa.nt
                    ii = ti == j;
                    bin(:,:,:,j) = fa.binDwell(x(ii), y(ii), tau(ii));
                end
                if (nargin >= 6 && ~isempty(at))
                    for j = 1:size(bin,4)
                        bin(:,:,:,j) = fa.applyAffine(bin(:,:,:,j), at(:,:,min(j,size(at,3))));
                    end
                end
                return;
            end
 
            maxDeltaX = max(abs(diff(x)))/fa.dx;
            maxDeltaY = max(abs(diff(y)))/fa.dy;
            expansion = ceil(min(fa.maxDwellExpansion, max(maxDeltaX, maxDeltaY)));
            if (expansion > 1)
                x = interpft(x, expansion*length(x));
                y = interpft(y, expansion*length(y));
                tau = interpft(tau/expansion, expansion*length(tau));
            end
            bin = fa.stretchDwellInZ(fa.accumVals(tau, x, y)); 
            if (nargin >= 6 && ~isempty(at))
                bin = fa.applyAffine(bin, at);
            end
        end
        
        function bin3 = stretchDwellInZ (fa, bin2)
            if (fa.nz < 1)
                bin3 = bin2;
                return;
            end
            if (fa.dxdz == 0 && fa.dydz == 0)
                zmin = min(fa.zOfP(linspace(-pi,pi, 361))); %just in case there's something weird and it's not just cos(p)
                zmax = max(fa.zOfP(linspace(-pi,pi, 361))); %just in case there's something weird and it's not just cos(p)
                valid = (fa.zedges(1:end-1) >= zmin) & (fa.zedges(2:end) < zmax);
                bin3 = repmat(bin2, [1 1 fa.nz])/sum(valid);
                bin3(:,:,~valid) = 0;
                return;
            end
            zx = fa.z0 + (0:(fa.nz-1))*fa.dz;
            bin3 = FastAssembler.fftTranslate(bin2, fa.dxdz*zx, fa.dydz*zx)/fa.nz;
        end
        
        function sz = getImSize(fa)
                sz = [fa.nx, fa.ny, fa.nz, fa.nt];
        end
        
        function [xl,yl,zl] = getImBoundary(fa)
            xl = fa.xedges([1 end]);
            yl = fa.yedges([1 end]);
            zl = fa.zedges([1 end]);
        end
        
        function cm = getImCenter(fa)
            [xl,yl,zl] = fa.getImBoundary;
            cm = [mean(xl) mean(yl) mean(zl)];
        end
        
        function ndim_out = applyAffine(fa, ndim, at)
            if (all(size(at) == [3,4]))
                at = [at; 0 0 0 1]';
            end
            if (all(size(at) == [4,3]))
                at = [at [0; 0; 0; 1]];
            end
            if (any(at(1:3,4) ~= 0) && all(at(4,:) == [0 0 0 1])) %this is the form if y = at*x instead of y = x*at
                at = at';
            end
            
            tform = affine3d(at);
            sz = fa.getImSize;
            [xc,yc,zc] = fa.getBinCenters();
            ra = imref3d(sz, xc([1 end]), yc([1 end]), zc([1 end]));
            ndim_out = pagetranspose(imwarp(pagetranspose(ndim), ra, tform, 'OutputView', ra));
        end
        
      
    end
    
    methods(Static)
        function [ind, inrange] = valToInd(x, x0, dx, nx)
            if (nx < 1)
                ind = ones(size(x));
                inrange = true(size(x));
                return;
            end
            ind = floor((x - x0)/dx + 1); %changed from round
            inrange = ind >= 1 & ind <= nx; 
        end
        
        function transims = fftTranslate(im, deltax, deltay)
            %adapted from example accompanying dftregistration.m on matlab
            %file exchange
            %assumes image is ndgrid oriented: ie im(x,y) not im(row, col)
            
            %zero-pad
            
            im_nx = size(im,1);
            im_ny = size(im,2);
            nx = 2^nextpow2(im_nx + max(abs(deltax)));
            ny = 2^nextpow2(im_ny + max(abs(deltay)));
            x0 = floor((nx-im_nx)/2);
            y0 = floor((ny-im_ny)/2);
            f = zeros(nx,ny);
            f(x0 + (1:im_nx), y0 + (1:im_ny)) = im;            
            Nx = ifftshift((-fix(nx/2):ceil(nx/2)-1));
            Ny = ifftshift((-fix(ny/2):ceil(ny/2)-1));
            [Nx,Ny] = ndgrid(Nx,Ny);
            transims = zeros(im_nx, im_ny, length(deltax));
            
            f = fft2(f);
            for j = 1:length(deltax)             
                g = ifft2(f.*exp(-1i*2*pi*(deltax(j)*Nx/nx+deltay(j)*Ny/ny)));
                transims(:,:,j) = abs(g(x0 + (1:im_nx), y0 + (1:im_ny)));
            end
           
        end
        
        function [xedges, x0, dx, nx] = regularizeGrid(xedges,makeEdgesFromCenters)
            xedges = unique(xedges);
            x0 = xedges(1);
            dx = xedges(2)-xedges(1);
            nx = round((xedges(end)-x0)/dx);
            if (makeEdgesFromCenters)
                x0 = x0 - dx/2;
                nx = nx+1;
            end
            xedges = x0 + (0:(nx))*dx;
        end
        
        function r = getBinRank(bin, subs)
            r = NaN(size(subs,1), 1);
            valid = true(size(r));
            for j = 1:size(subs,2)
                valid = valid & subs(:,j) >= 1 & subs(:,j) <= size(bin,j);
            end
            [~,I] = sort(bin(:), 'descend');
%             ii = sub2ind(size(bin), subs(valid,1), subs(valid,2));
              for j = 1:size(subs,2)
                 args{j} = subs(valid,j);
             end
             ii = sub2ind(size(bin), args{:});
            r(valid) = interp1(I, 1:length(I), ii, 'nearest', NaN);
            
        end
        
    end
    
end

