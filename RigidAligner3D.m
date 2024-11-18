classdef RigidAligner3D
    % Calculates image values, derivatives, and Hessians 
    % at specified locations where photon was observed
    
    properties (Constant)
        
    end
    
    properties
        invertEval = true; %allows minimizing F using eval directly        
        fitParameters; %1 x nparams 
        fitAffine = [1 0 0 0;0 1 0 0;0 0 1 0]; %3D affine trans = (A|delta), A is 3x3, delta is 3x1;
        startTime;
        endTime;
        funPerPhoton;
        ptime; %1 x N - time
        ploc; %3 x N - [x;y;z] - for now keep z = 0
        pphase; %1 x N - tag lens phase for each photon 


        dtime; %1 x NN - time coordinate (not duration) of dwell
        dloc; %3 x NN - dwell coordinate - for now keep z = 0
        dphase; %1 x NN - tag lens phase for each dwell coordinate
        dweight; %1 x NN - dwell duration * rate_scaling for template image (e.g. if tau = 25 ns and rate_scale = 10^7 Hz, dweight = .25)

        z_scale = 25; %z-axis in microns
        
        xrange;
        yrange;
        czrange;
        rsqrange = Inf;
        
        binnedRaw;
        
        normByPhotons = true;
        
        mifit;
        mifit2D;      
        
        lb = [-pi/4,       -pi/6, -pi/6, -5, -5, -10]; 
        ub = [9*pi/4, pi/6, pi/6, 5, 5, 10];

        nparams = 6; %theta_z, theta_y, theta_x, dx, dy, dz
        
        hasdwell = false;
    end
    
    properties (Access = public)
        plocvec; %4 x N, [x;y;z;1] (last row all 1s)
        pweightvec; %1 X N, weight applied to each photon
        
        dlocvec; %as for ploc
        dweightvec; %as for pweight

        %image function and derivatives, all gridded interpolant
        immodel;
        
        initAffine; %for debugging
       
    end
    
    
    methods
       
        function ra = RigidAligner3D(immodel, photonloc, photonphase, photontime,  z_scale)
            %rfc = RasterFitChunk(photonloc, immodel)
            %photonloc = [xphoton;yphoton;time] (time is optional but used
            %to split rfc into smaller components)
            %immodel has fields u, v, F - function of (u,v) -- u is
            %horizontal axis (x-axis) -- ie F is functional, not matrix form
            %immodel can have Fextrap 
            if (nargin < 3)
                return;
            end
            if (nargin < 4)
                photontime = [];
            end
            if (nargin >= 5)
               ra.z_scale = z_scale;
            end
            
            ra = ra.setPhotonLoc(photonloc, photonphase, photontime);
            ra = ra.createImModel(immodel);
        end
        
        function ra = affineTransformPoints(ra, at)
            pl = at*[ra.ploc;ones([1 size(ra.ploc,2)])];
            ra.ploc = pl(1:3,:);
        end
        
        function fa = getFastAssembler(ra, nbins, at)

            %transform bounding box by at if provided (xrange, yrange,
            %czrange)
            if (numel(nbins) < 3)
                nbins = nbins(1)*[1 1 1];
            end
            if (nargin < 3 || isempty(at))
                fa = FastAssembler(linspace(min(ra.xrange), max(ra.xrange), nbins(1)+1), linspace(min(ra.yrange), max(ra.yrange), nbins(2)+1), ra.z_scale*linspace(min(ra.czrange), max(ra.czrange), nbins(3)+1));
            else
                [xx,yy,zz] = ndgrid(ra.xrange, ra.yrange, ra.czrange*ra.z_scale);
                bb = at*([xx(:) yy(:) zz(:) ones(size(xx(:)))]');
                fa = FastAssembler(linspace(min(bb(1,:)), max(bb(1,:)), nbins(1)+1), linspace(min(bb(2,:)), max(bb(2,:)), nbins(2)+1), linspace(min(bb(3,:)), max(bb(3,:)), nbins(3)+1));
            end
                
            fa.zOfP = @(p) ra.z_scale*cos(p);
            
        end

        function fa = getFastAssemblerBinsize(ra, binsize, at)
            if (nargin < 3 || isempty(at))              
                nbins = [max(10, ceil(abs(diff(ra.xrange)/binsize(1)))), ceil(max(10, abs(diff(ra.yrange)/binsize(2)))), ceil(max(10, abs(ra.z_scale*diff(ra.czrange)/binsize(3))))];
                fa = ra.getFastAssembler(nbins);
            else
                [xx,yy,zz] = ndgrid(ra.xrange, ra.yrange, ra.czrange*ra.z_scale);
                bb = at*([xx(:) yy(:) zz(:) ones(size(xx(:)))]');
                xr = [min(bb(1,:)), max(bb(1,:))];
                yr = [min(bb(2,:)), max(bb(2,:))];
                zr = [min(bb(3,:)), max(bb(3,:))];
                nbins = [max(10, ceil(abs(diff(xr)/binsize(1)))), ceil(max(10, abs(diff(yr)/binsize(2)))), ceil(max(10, abs(diff(zr)/binsize(3))))];
                fa = FastAssembler(linspace(xr(1), xr(2), nbins(1)+1), linspace(yr(1), yr(2), nbins(2)+1), linspace(zr(1), zr(2), nbins(3)+1));
            end
        end
        
        function ra = setPhotonLoc(ra, photonloc, photonphase, photontime)
            ra.ploc = photonloc(1:3,:);  
            ra.pphase = photonphase;
            if (nargin < 4)
                photontime = [];
            end
            ra.ptime = photontime;
            
            ra.startTime = min(ra.ptime);
            ra.endTime = max(ra.ptime);
                       
        end
        
        function ra = setDwellLoc(ra, dwellloc, dwellphase, dwellweight, dwelltime)
            if (isempty(dwellloc))
                return;
            end
            ra.dloc = dwellloc(1:3,:);
            ra.dphase = dwellphase;
            if (nargin < 5)
                dwelltime = [];
            end
            ra.dtime = dwelltime;
            ra.dweight = dwellweight;
            ra.hasdwell = true;

        end

        function fitMontage (ra, autozoom)
            if (nargin < 2 || isempty(autozoom))
                autozoom = true;
            end


            if (ra.hasdwell)
                [counts,dwell] = ra.binCountsAndDwell();
                [counts2,dwell2] = ra.binCountsAndDwell([],false);
                im = imblur(counts,2)./imblur(dwell,2);
                im2 = imblur(counts2,2)./imblur(dwell2,2);
            else
                im = ra.makeImage();
                im2 = ra.makeImage([], false);
            end

            ul = ra.immodel.u([find(any(counts > 0, [2 3]), 1, 'first') find(any(counts > 0, [2 3]), 1, 'last')]);
            vl = ra.immodel.v([find(any(counts > 0, [1 3]), 1, 'first') find(any(counts > 0, [1 3]), 1, 'last')]);
            wl = ra.immodel.w([find(any(counts > 0, [1 2]), 1, 'first') find(any(counts > 0, [1 2]), 1, 'last')]);

            axdim = [2 1; 3 2; 3 1];
            projdim = [3; 1; 2];
            axvec = {ra.immodel.u, ra.immodel.v, ra.immodel.w};
            lims = {ul, vl, wl};
            axlabel = {'x','y','z'};
            maxpts = 40;
            s = 0.5;
            for qq = 1:3
                figure(); clf();
                refim = squeeze(mean(ra.immodel.F_im,projdim(qq),'omitnan'));
                p = FastPeakFind(refim, percentile(refim(:), .05));

                if (ra.hasdwell)
                    ims = {squeeze(imblur(sum(counts,projdim(qq),'omitnan'),s)./imblur(sum(dwell,projdim(qq),'omitnan'),s)),...
                        refim,...
                        squeeze(imblur(sum(counts2,projdim(qq),'omitnan'),s)./imblur(sum(dwell2,projdim(qq),'omitnan'),s))};
                else
                    ims = {squeeze(sum(im,projdim(qq),'omitnan')),squeeze(sum(ra.immodel.F_im,projdim(qq),'omitnan')), squeeze(sum(im2,projdim(qq),'omitnan'))};
                end
    
                if (size(p,2) > maxpts)
                    [~,I] = sort(refim(sub2ind(size(refim), p(2,:),p(1,:))), 'descend');
                    p = p(:,I(1:maxpts));
                end
                for j = 1:3
                    ax(j) = subplot(1,3,j); %#ok<AGROW> 
                    pcolor(axvec{axdim(qq,1)}, axvec{axdim(qq,2)}, ims{j}); 
                    shading flat; axis equal
                    xlim(lims{axdim(qq,1)});
                    ylim(lims{axdim(qq,2)});
                    if (ra.hasdwell)
                        v = ims{j};
%                         clim(percentile(v(isfinite(v)),[.05 .99]));
                        set(gca,'CLim',percentile(v(isfinite(v)),[.05 .99])); % MATLAB renamed caxis to clim in R2022a and it's really annoying
                    end

                    axis equal;
                    hold on;
                    x = axvec{axdim(qq,1)}(p(1,:));
                    y = axvec{axdim(qq,2)}(p(2,:));
                    valid = x > ax(j).XLim(1) & x < ax(j).XLim(2) & y > ax(j).YLim(1) & y < ax(j).YLim(2);
                    x = x(valid); y = y(valid);
                    for k = 1:length(x)
                        text(x(k), y(k), num2str(k), 'Color','m');
                    end
                end
                xlabel(ax(2), axlabel(axdim(qq,1)));
                ylabel(ax(1), axlabel(axdim(qq,2)));

            end
           
           return
            

            
            figure(2);clf(2);

            refim = squeeze(mean(ra.immodel.F_im,3,'omitnan'));
            p = FastPeakFind(refim, percentile(refim(:), .05));

            if (size(p,2) > maxpts)
                [~,I] = sort(refim(sub2ind(size(refim), p(2,:),p(1,:))), 'descend');
                p = p(:,I(1:maxpts));
            end

            if (ra.hasdwell)
                ims = {squeeze(imblur(sum(counts,3,'omitnan'),s)./imblur(sum(dwell,3,'omitnan'),s)),...
                    squeeze(mean(ra.immodel.F_im,3)),...
                    squeeze(imblur(sum(counts2,3,'omitnan'),s)./imblur(sum(dwell2,3,'omitnan'),s))};
            else
                ims = {squeeze(sum(im,3,'omitnan')),squeeze(sum(ra.immodel.F_im,3,'omitnan')), squeeze(sum(im2,3,'omitnan'))};
            end

            for j = 1:3
                subplot(1,3,j); pcolor(ra.immodel.v, ra.immodel.u, ims{j}); shading flat; axis equal
                xlim(vl);
                ylim(ul);
                if (ra.hasdwell)
                    v = ims{j};
%                     clim(percentile(v(isfinite(v)),[.05 .99]));
                    set(gca,'CLim',percentile(v(isfinite(v)),[.05 .99]));
                end
                %clim([min(sum(ra.immodel.F_im,2,'omitnan'),[],'all') max(sum(ra.immodel.F_im,2,'omitnan'),[],'all')]);
        
                axis equal;
                hold on;
                for k = 1:length(p)
                    text(ra.immodel.v(p(1,k)), ra.immodel.u(p(2,k)), num2str(k),'Color','m');
                end
            end

        end

        function [im, ax] = makeImage(ra, ax, usefit)
            % function [im, ax] = makeImage(ra, ax, usefit)
            % if ax not provided, uses immodel ax
            if (nargin < 2 || isempty(ax))
                ax.x = ra.immodel.u;
                ax.y = ra.immodel.v;
                ax.z = ra.immodel.w;
            end
            if (nargin < 3 || isempty(usefit))
                usefit = true;
            end
            if (usefit)
                at = ra.fitAffine;
            else
                at = [1 0 0 0; 0 1 0 0; 0 0 1 0];
            end
            fa = FastAssembler(ax.x, ax.y, ax.z,[],true);
            fa.zOfP = @(p) ra.z_scale*cos(p);
            valid = ra.ploc(1,:) >= min(ra.xrange) & ra.ploc(1,:) <= max(ra.xrange) & ra.ploc(2,:) >= min(ra.yrange) & ra.ploc(2,:) <= max(ra.yrange) & cos(ra.pphase) >= min(ra.czrange) & cos(ra.pphase) <= max(ra.czrange);
            valid = valid&(ra.ploc(1,:).^2 + ra.ploc(2,:).^2 < ra.rsqrange);

            im = fa.binPhotons(ra.ploc(1,valid), ra.ploc(2,valid), ra.pphase(valid),[],at);

        end
        function [counts,dwell] = binCountsAndDwell (ra, fa, usefit)
            if (nargin < 2 || isempty(fa))               
                fa = FastAssembler(ra.immodel.u, ra.immodel.v, ra.immodel.w,[],true);
            end
            
            if (nargin < 3 || isempty(usefit))
                usefit = true;
            end
            if (usefit)
                at = ra.fitAffine;
            else
                at = [1 0 0 0; 0 1 0 0; 0 0 1 0];
            end
            valid = ra.ploc(1,:) >= min(ra.xrange) & ra.ploc(1,:) <= max(ra.xrange) & ra.ploc(2,:) >= min(ra.yrange) & ra.ploc(2,:) <= max(ra.yrange) & cos(ra.pphase) >= min(ra.czrange) & cos(ra.pphase) <= max(ra.czrange);
            valid = valid&(ra.ploc(1,:).^2 + ra.ploc(2,:).^2 < ra.rsqrange);
            counts = fa.binPhotons(ra.ploc(1,valid), ra.ploc(2,valid), ra.pphase(valid),[],at);
            if (ra.hasdwell)
                valid = ra.dloc(1,:) >= min(ra.xrange) & ra.dloc(1,:) <= max(ra.xrange) & ra.dloc(2,:) >= min(ra.yrange) & ra.dloc(2,:) <= max(ra.yrange) & cos(ra.dphase) >= min(ra.czrange) & cos(ra.dphase) <= max(ra.czrange);
                valid = valid&(ra.dloc(1,:).^2 + ra.dloc(2,:).^2 < ra.rsqrange);
                dwell = fa.binWithWeightAndAffine(ra.dloc(1,valid), ra.dloc(2,valid), ra.dphase(valid),ra.dweight,at);
            else
                dwell =0 * counts;
            end

        end

        function [g,gd] = gradientTest(ra, pvec, delta)
            [v,g] = ra.eval(pvec);
            if (nargin <3)
                delta = 1e-3;
            end
            for j = 1:length(pvec)
                d = 0*pvec;
                d(j) = delta;
                v2 = ra.eval(pvec + d);
                gd(j) = (v2-v)/delta;

            end
        end
        
        function ra = binPhotonsByBinsize(ra, binsize)
             if (nargin < 2 || isempty(binsize))
                binsize = [median(diff(ra.immodel.u)) median(diff(ra.immodel.v)) median(diff(ra.immodel.w))];            
             end
             if (length(binsize) < 3)
                 binsize = binsize(1) * [1 1 1];
             end
             nbins = [max(10, abs(diff(ra.xrange)/binsize(1))), max(10, abs(diff(ra.yrange)/binsize(2))), max(10, abs(ra.z_scale*diff(ra.czrange)/binsize(3)))];
             ra = ra.binPhotons(nbins);
        end
        
        function ra = binPhotons (ra, nbins)
            % function rfc = binPhotons (rfc)
            % function rfc = binPhotons (rfc, nbins)
            % function rfc = binPhotons (rfc, binx, biny, binz)

            fa = ra.getFastAssembler(nbins);
            
            valid = ra.ploc(1,:) >= min(ra.xrange) & ra.ploc(1,:) <= max(ra.xrange) & ra.ploc(2,:) >= min(ra.yrange) & ra.ploc(2,:) <= max(ra.yrange) & cos(ra.pphase) >= min(ra.czrange) & cos(ra.pphase) <= max(ra.czrange);
            valid = valid&(ra.ploc(1,:).^2 + ra.ploc(2,:).^2 < ra.rsqrange);
            
            ra.binnedRaw.im = fa.binPhotons(ra.ploc(1,valid), ra.ploc(2,valid), ra.pphase(valid));
            [ra.binnedRaw.x, ra.binnedRaw.y, ra.binnedRaw.z] = fa.getBinCenters();
            
            w = ra.binnedRaw.im;
            [x,y,z] = ndgrid(ra.binnedRaw.x, ra.binnedRaw.y, ra.binnedRaw.z);
            
            w = w(:);
            ra.plocvec = [x(w>0)';y(w>0)';z(w>0)']; ra.plocvec(4,:) = 1;
            ra.pweightvec = w(w>0)';
            
            if (ra.hasdwell)
                valid = ra.dloc(1,:) >= min(ra.xrange) & ra.dloc(1,:) <= max(ra.xrange) & ra.dloc(2,:) >= min(ra.yrange) & ra.dloc(2,:) <= max(ra.yrange) & cos(ra.dphase) >= min(ra.czrange) & cos(ra.dphase) <= max(ra.czrange);
                valid = valid&(ra.dloc(1,:).^2 + ra.dloc(2,:).^2 < ra.rsqrange);
            

                ra.binnedRaw.dim = fa.binWithWeightAndAffine (ra.dloc(1,valid), ra.dloc(2,valid), ra.dphase(valid), ra.dweight);

            
                w = ra.binnedRaw.dim(:);
                ra.dlocvec = [x(w>0)';y(w>0)';z(w>0)']; ra.dlocvec(4,:) = 1;
                ra.dweightvec = w(w>0)';
            end

            
        end
        
        function ra = estimateMutualInformation(ra)
            q = @(x) round(255*(x(:)-min(x(:)))./(max(x(:))-min(x(:))));
            im = ra.makeImage();
            if (isempty(im))
                ra.mifit = 0;
                ra.mifit2D = 0;
                return;
            end
            ra.mifit = MI_GG(q(im), q(ra.immodel.F_im));
            ra.mifit2D = MI_GG(q(sum(im,3)), q(sum(ra.immodel.F_im,3)));
        end
        
        function [val, grad] = eval(ra, pvec)
            %[val, grad] = eval(rfc, pvec)
            %
            % part of objective function for fitter
            % calculates image energy and gradients 
            % does not calculate continuity terms
            
            % if rfc is a single element,
            %   pvec = [theta_z theta_y theta_x dx dy dz] (1x6)
            %
           
            eval2D = false;
            
            if (isempty(ra.pweightvec) || sum(ra.pweightvec) == 0 || ~isfinite(sum(ra.pweightvec)))
                val = 0;
                grad = 0*pvec;
                return;
            end
            
            if (length(pvec) == 3)
                eval2D = true;
                grad = 0 * pvec;
                pvec = [pvec(1) 0 0 pvec(2) pvec(3) 0];
            end
            
            if (ra.invertEval)
                m = -1;
            else
                m = 1;
            end
            if (ra.normByPhotons)
                m = m/sum(ra.pweightvec);
            end
            
            [tmat, gradmat] = RigidAligner3D.rigidTransform(pvec);
            
            
            newpts = tmat*ra.plocvec; %3x4 * 4xN --> 3xN
            x = newpts(1,:);
            y = newpts(2,:);
            z = newpts(3,:);

            if (ra.hasdwell)
                dpts = tmat*ra.dlocvec; %3x4 * 4xN --> 3xN
            end

            if (ra.hasdwell)
                 if (eval2D)
                    f = ra.immodel.F_2D(x,y);
                    val = m * (sum(ra.pweightvec.*log(f)) - sum(ra.dweightvec.*ra.immodel.F_2D(dpts(1,:),dpts(2,:))));
                 else
                     f = ra.immodel.F(x,y,z);
                    val = m * (sum(ra.pweightvec.*log(f),'all','omitnan') - sum(ra.dweightvec.*ra.immodel.F(dpts(1,:),dpts(2,:),dpts(3,:)),'all','omitnan'));
                 end
                 if (nargout < 2)
                    return;
                end
                if (eval2D)
                    gu = m*(ra.pweightvec.*ra.immodel.dFdu_2D(x,y)./f);
                    gv = m*(ra.pweightvec.*ra.immodel.dFdv_2D(x,y)./f);
                    gw = 0*gv;
                else
                    gu = m*(ra.pweightvec.*ra.immodel.dFdu(x,y,z)./f);
                    gv = m*(ra.pweightvec.*ra.immodel.dFdv(x,y,z)./f);
                    gw = m*(ra.pweightvec.*ra.immodel.dFdw(x,y,z)./f);
                end
                pgrad = [gu(:)';gv(:)';gw(:)'];
                if (eval2D)
                    gu = m*(ra.dweightvec.*ra.immodel.dFdu_2D(dpts(1,:),dpts(2,:)));
                    gv = m*(ra.dweightvec.*ra.immodel.dFdv_2D(dpts(1,:),dpts(2,:)));
                    gw = 0*gv;
                else
                    gu = m*(ra.dweightvec.*ra.immodel.dFdu(dpts(1,:),dpts(2,:), dpts(3,:)));
                    gv = m*(ra.dweightvec.*ra.immodel.dFdv(dpts(1,:),dpts(2,:), dpts(3,:)));
                    gw = m*(ra.dweightvec.*ra.immodel.dFdw(dpts(1,:),dpts(2,:), dpts(3,:)));
                end
                dgrad = [gu(:)';gv(:)';gw(:)'];
                if (eval2D)
                    grad(1) = sum(pgrad .* (gradmat(:,:,1)*ra.plocvec), 'all','omitnan') - sum(dgrad .* (gradmat(:,:,1)*ra.dlocvec), 'all','omitnan');
                    grad(2) = sum(pgrad .* (gradmat(:,:,4)*ra.plocvec), 'all','omitnan') - sum(dgrad .* (gradmat(:,:,4)*ra.dlocvec), 'all','omitnan');
                    grad(3) = sum(pgrad .* (gradmat(:,:,5)*ra.plocvec), 'all','omitnan') - sum(dgrad .* (gradmat(:,:,5)*ra.dlocvec), 'all','omitnan');
                else
                    grad = pvec*0;
                    for j = 1:6
                        grad(j) = sum(pgrad .* (gradmat(:,:,j)*ra.plocvec), 'all','omitnan') - sum(dgrad .* (gradmat(:,:,j)*ra.dlocvec), 'all','omitnan');
                    end
                end
            else
                if (eval2D)
                    val = m * sum(ra.pweightvec.*ra.immodel.F_2D(x,y));
                else
                    val = m * sum(ra.pweightvec.*ra.immodel.F(x,y,z));
                end
                if (nargout < 2)
                    return;
                end
                if (eval2D)
                    gu = m*(ra.pweightvec.*ra.immodel.dFdu_2D(x,y));
                    gv = m*(ra.pweightvec.*ra.immodel.dFdv_2D(x,y));
                    gw = 0*gv;
                else
                    gu = m*(ra.pweightvec.*ra.immodel.dFdu(x,y,z));
                    gv = m*(ra.pweightvec.*ra.immodel.dFdv(x,y,z));
                    gw = m*(ra.pweightvec.*ra.immodel.dFdw(x,y,z));
                end
                imgrad = [gu(:)';gv(:)';gw(:)'];

                if (eval2D)
                    grad(1) = sum(imgrad .* (gradmat(:,:,1)*ra.plocvec), 'all');
                    grad(2) = sum(imgrad .* (gradmat(:,:,4)*ra.plocvec), 'all');
                    grad(3) = sum(imgrad .* (gradmat(:,:,5)*ra.plocvec), 'all');
                else
                    grad = pvec*0;
                    for j = 1:6
                        grad(j) = sum(imgrad .* (gradmat(:,:,j)*ra.plocvec), 'all');
                    end
                end
            end
        end
        
        function [theta_z,fval] = orientRotationBounded(ra, theta_bounds, ndivisions)
            if (nargin < 3 || isempty(ndivisions))
                ndivisions = 3;
            end
            ra = ra.binPhotonsByBinsize();%uses immodel binning
            
            funr = @(t) ra.eval([t 0 0]); 
            
%             tx = linspace(0,2*pi,100);
%             for j = 1:length(tx)
%                 f(j) = funr(tx(j));
%             end
%             plot(180/pi * tx, f); 
            
            [theta_z,fval] = fminbnd_robust(funr, theta_bounds(1), theta_bounds(2), ndivisions);
        end
        
       
        function [theta_z, dx,dy] = findInitial2DOrientation (ra, binsize) 
            if (nargin < 2 || isempty(binsize))
                binsize = 2; %microns
            end
            ra = ra.binPhotonsByBinsize(binsize);
            
            
            %find best initial rotation by minimizing just with regard to
            %theta around z-axis
            funr = @(t) ra.eval([t 0 0]);
%             op = optimset('fminbnd');
%             op.Display = 'off';
%             op.MaxIter = 100;
%             op.TolX = 0.1;
%             tx = linspace(-pi/4,2*pi+pi/4, 11);
%             for j = 1:(length(tx) -1)
%                 [~,f(j)] = fminbnd(funr, tx(j), tx(j+1), op); %#ok<AGROW>
%             end
%             [~,I] = min(f);
%             op.TolX = 0.01;
%             op.MaxIter = 500;
%            theta_init = mod(fminbnd(funr, tx(I), tx(I+1),op), 2*pi);
            theta_init = mod(fminbnd_robust(funr, -pi/8, 2*pi), 2*pi);
            
            lb = ra.lb([1 4 5]);
            ub = ra.ub([1 4 5]);
            x0 = [theta_init, 0, 0];
            %fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB)
            op = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Display', 'notify-detailed', 'MaxFunctionEvaluations', 3e4);

            %here to throw an error if there's a problem with initial
            %conditions
            [v0,g0] = ra.eval(x0);
            
            %find best 2D rotation and translation starting at initial
            %rotation guess
            [x1,v1] = fmincon(@(x) ra.eval(x), x0, [], [], [], [], lb, ub, [], op);
            
            x0_pi = [mod(theta_init + pi, 2*pi), 0, 0];
            
            %find best 2D rotation and tranlsation starting 180 degrees
            %away from best guess
            [x2,v2] = fmincon(@(x) ra.eval(x), x0_pi, [], [], [], [], lb, ub, [], op);

            %choose best fit
            if (v1 < v2)
                x = x1;
            else 
                x = x2;
            end
            theta_z = x(1);
            dx = x(2);
            dy = x(3);        
        end
        
        function [pvecinit, pvecinit180] = initialMomentsHorizontal(ra, percentile_cut, binsize)
            if (nargin < 2 || isempty(percentile_cut))
                percentile_cut = 0.95;
            end
            if (percentile_cut > 1)
                percentile_cut = percentile_cut / 100;
            end
            if (nargin < 3 || isempty(binsize))
                binsize = [median(diff(ra.immodel.u)) median(diff(ra.immodel.v)) median(diff(ra.immodel.w))];
            end
            ra = ra.binPhotonsByBinsize(binsize);
            inds = ra.pweightvec >= percentile(ra.pweightvec, percentile_cut);
            [~,pvecinit,~,pvecinit180] = orientPoints3D_no_delta(ra.plocvec(1,inds), ra.plocvec(2,inds), ra.plocvec(3,inds), ra.pweightvec(inds));
            
%             ra.fitParameters = pvecinit;
%             ra.fitAffine = at;
%             [im, ax] = ra.makeImage([],true);
%             subplot(2,2,1); pcolor(ax.x, ax.y, sum(im,3)'); shading flat; axis equal;
%             subplot(2,2,2); pcolor(ax.z, ax.y, squeeze(sum(im,1))); shading flat; axis equal;
%             subplot(2,2,3); pcolor(ra.immodel.u, ra.immodel.v, ra.immodel.F_2D_im'); shading flat; axis equal;
%             subplot(2,2,4); pcolor(ra.immodel.w, ra.immodel.v, squeeze(sum(ra.immodel.F_im,1))); shading flat; axis equal;
                 
        end
         
        function [ra,pvec] = fitRigidFromScratch(ra, percentile_cut, binsize)
            %recommend using a coarse grid and template for this
             if (nargin < 2 || isempty(percentile_cut))
                percentile_cut = 0.95;
            end
            if (percentile_cut > 1)
                percentile_cut = percentile_cut / 100;
            end
            if (nargin < 3 || isempty(binsize))
                binsize = [median(diff(ra.immodel.u)) median(diff(ra.immodel.v)) median(diff(ra.immodel.w))];
            end
            [pvecinit, pvecinit180] = ra.initialMomentsHorizontal(percentile_cut, binsize);
            
          
            [ra1,pvec1] = ra.fitRigid(pvecinit, binsize);
            [ra,pvec] = ra.fitRigid(pvecinit180, binsize);
%             
%             figure(1);
%             [im,ax] = ra1.makeImage();
%             subplot(2,2,1); pcolor(ax.x, ax.y, sum(im,3)'); shading flat; axis equal;
%              subplot(2,2,2); pcolor(ax.z, ax.y, squeeze(sum(im,1))); shading flat; axis equal;
%              subplot(2,2,3); pcolor(ra.immodel.u, ra.immodel.v, ra.immodel.F_2D_im'); shading flat; axis equal;
%              subplot(2,2,4); pcolor(ra.immodel.w, ra.immodel.v, squeeze(sum(ra.immodel.F_im,1))); shading flat; axis equal;
%             
%               figure(2);
%             [im,ax] = ra.makeImage();
%             subplot(2,2,1); pcolor(ax.x, ax.y, sum(im,3)'); shading flat; axis equal;
%              subplot(2,2,2); pcolor(ax.z, ax.y, squeeze(sum(im,1))); shading flat; axis equal;
%              subplot(2,2,3); pcolor(ra.immodel.u, ra.immodel.v, ra.immodel.F_2D_im'); shading flat; axis equal;
%              subplot(2,2,4); pcolor(ra.immodel.w, ra.immodel.v, squeeze(sum(ra.immodel.F_im,1))); shading flat; axis equal;
%             
%              ra1.funPerPhoton
%              ra.funPerPhoton
%              
            if (ra.funPerPhoton > ra1.funPerPhoton)
                ra = ra1;
                pvec = pvec1;
            end
            

        end
        
        function [ra,pvec] = oldfitRigidFromScratch(ra, bsinit, bsfinal)
            
            if (nargin < 2 || isempty(bsinit))
                bsinit = 2;
            end
            if (nargin < 3 || isempty(bsfinal))
                bsfinal = min(median(diff(ra.immodel.u)), median(diff(ra.immodel.v)));
            end
            [theta_z, dx,dy] = findInitial2DOrientation (ra, bsinit);
            pvec_init = [theta_z, 0, 0, dx, dy, 0];
            [ra, pvec] = ra.fitRigid(pvec_init, bsinit);
            if (any(bsfinal < bsinit))
                [ra, pvec] = ra.fitRigid(pvec, bsfinal);
            end
            
        end
        
        function [ra, pvec] = fitRigid (ra, pvec_init, binsize)
            if (nargin < 3)
                binsize = []; 
            end
            if (nargin < 2 || isempty(pvec_init))
                 [theta_z, dx,dy] = findInitial2DOrientation (ra, binsize);
                 pvec_init = [theta_z, 0, 0, dx, dy, 0];
            end
            if (isempty(binsize))
                binsize = [median(diff(ra.immodel.u)) median(diff(ra.immodel.v)) median(diff(ra.immodel.w))];            
            end
            ra = ra.binPhotonsByBinsize(binsize);
%             nbins = [max(10, abs(diff(ra.xrange)/binsize)), max(10, abs(diff(ra.yrange)/binsize)), max(10, abs(ra.z_scale*diff(ra.czrange)/binsize))];
%             ra = ra.binPhotons(nbins);
            op = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Display', 'notify', 'MaxFunctionEvaluations', 3e4, 'StepTolerance', 1e-8);


            [pvec, val] = fmincon(@(x) ra.eval(x), pvec_init, [], [], [], [], ra.lb, ra.ub, [], op);
            pvec(1) = mod(pvec(1), 2*pi);
            ra.fitParameters = pvec;
            if (ra.normByPhotons)
                ra.funPerPhoton = val;
            else
                ra.funPerPhoton = val/(sum(ra.pweightvec));
            end
            ra.fitAffine = RigidAligner3D.rigidTransform(pvec);
        end

        function ra = createImModel(ra, template)
            if (isa(template, 'Immodel'))
                [ra.immodel] = deal(template);
            else
                [ra.immodel] = deal(Immodel(template));
            end        
        end
        
        
        
        
    end
    
    methods(Static)
        
     
        
        function pvec2 = pvec1to2 (pvec1)
            %reshapes 1D parameter list to be nrfcXnparams sized
            Np = RasterFitterChunk3D.nparams;
            nr = numel(pvec1)/Np;
            pvec2 = reshape(pvec1, nr, Np);          
        end
        function pvec1 = pvec2to1 (pvec2)
            pvec1 = reshape(pvec2, 1, []);
        end
        

       
        
        
    
        function [tmat,gradmat] = rigidTransform(theta_z, theta_y, theta_x, dx, dy, dz)
            %[tmat,gradmat] = rigidTransform(theta_z, theta_y, theta_x, dx, dy, dz)
            % rotate about z axis, then rotate about x axis, then rotate
            % about y axis
            %   tmat = transform matrix [3x4 affine]
            %   gradmat --> if y = tmat*[x;1], then gradmat.theta_z * (x_1 x_2 x_3 1)' is (dy_1/dtheta_z,
            %   dy_2/dtheta_z, dy_3/dtheta_z)' etc.
            %   or 
            %   gradmat(:,:,1) * (x_1 x_2 x_3 1)' is (dy_1/dpvec(1),
            %   dy_2/dpvec(1), dy_3/dpvec(1))' etc.
            
            if (nargin == 1)
                p = theta_z;
                theta_z = p(1);
                theta_y = p(2);
                theta_x = p(3);
                dx = p(4);
                dy = p(5); 
                dz = p(6);            
            end
            
            
            cz = cos(theta_z);
            sz = sin(theta_z);
            rz = [cz -sz 0; sz cz 0; 0 0 1];
            drz = [-sz -cz 0; cz -sz 0; 0 0 0];
            
            cy = cos(theta_y);
            sy = sin(theta_y);
            ry = [cy 0 sy;0 1 0; -sy 0 cy];
            dry = [-sy 0 cy; 0 0 0; -cy 0 -sy];
            
            cx = cos(theta_x);
            sx = sin(theta_x);
            rx = [1 0 0;0 cx -sx;0 sx cx];
            drx = [0 0 0; 0 -sx -cx; 0 cx -sx];
            
            tmat = ry * rx *  rz;
            tmat(:,4) = [dx;dy;dz];
            
            
           
            gradmat.theta_z = ry*rx*drz;
            gradmat.theta_y = dry*rx*rz;
            gradmat.theta_x = ry*drx*rz;
            
             gradmat.theta_z(:,4) = 0;
             gradmat.theta_y(:,4) = 0;
             gradmat.theta_x(:,4) = 0;
             
             gradmat.dx = zeros(3,4); gradmat.dx(1,4) = 1;
             gradmat.dy = zeros(3,4); gradmat.dy(2,4) = 1;
             gradmat.dz = zeros(3,4); gradmat.dz(3,4) = 1;
             
             
            if (nargin == 1)
                gradmat = cat(3, gradmat.theta_z, gradmat.theta_y, gradmat.theta_x, gradmat.dx, gradmat.dy, gradmat.dz);
            end
             
        end
        
        function pvec = rigidToPvec(at)

             if (all(size(at) == [3,4]))
                at = [at; 0 0 0 1];
            end
            if (all(size(at) == [4,3]))
                at = [at [0; 0; 0; 1]]';
            end
            if (any(at(4,1:3) ~= 0) && all(at(:,4) == [0 0 0 1])) %this is the form if y = x*at instead of y = x*at
                at = at';
            end
            attarget = at(1:3,:);

            
            %at =  delta*ry * rx *  rz;
            theta_z = mod(atan2(at(2,1), at(1,1)), 2*pi);
            cz = cos(theta_z);
            sz = sin(theta_z);
            rz = [cz -sz 0 0; sz cz 0 0; 0 0 1 0; 0 0 0 1];
            
            at = at/rz;
            theta_x = atan2(at(3,2), at(2,2));
            if (cos(theta_x) < 0) %greater than 90 degree rotation - so don't use as initial guess
                theta_x = 0;
            end
            cx = cos(theta_x);
            sx = sin(theta_x);
            rx = [1 0 0 0;0 cx -sx 0;0 sx cx 0; 0 0 0 1];

            at = at/rx;
            theta_y = atan2(at(1,3), at(1,1));
            if (cos(theta_y) < 0) %greater than 90 degree rotation - so don't use as initial guess
                theta_y = 0;
            end
            cy = cos(theta_y);
            sy = sin(theta_y);
            ry = [cy 0 sy 0;0 1 0 0; -sy 0 cy 0; 0 0 0 1];
            at = at/ry;
            
            pvec_init = [theta_z, theta_y, theta_x at(1:3,4)'];
            function l2n = atdist(p)
                a = RigidAligner3D.rigidTransform(p)-attarget;
                l2n = sum(a(:).^2);
            end
            lb = [-pi/4 -pi/2 -pi/2 -Inf -Inf -Inf];
            ub = [9*pi/4 pi/2 pi/2 +Inf +Inf +Inf];
            op = optimoptions('fmincon', 'display', 'notify-detailed', 'MaxFunctionEvaluations', 3e4, 'StepTolerance', 1e-8);

            
            pvec = fmincon(@atdist, pvec_init, [],[],[],[],lb,ub,[],op);
            pvec(1) = mod(pvec(1),2*pi);
            
          
            
            
        end
        
        
    end
    
end

