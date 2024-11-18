classdef RasterFitter
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
       
        czrange = [-2 2]; %cos(phase) - zo cutoff range for import; default -2 to 2 includes all possible phases and offsets 
        
        k_a_times_second = 10; %2x2 matrix spring constant per photon times second
        k_delta_times_second = 1 %translation spring constant per photon times second        
        
        display = 'none'
        prefix = ''
        verbose = 0
        stripTpAxes = [1 2]; %axes on which to remove turning points (edges)

        %must be provided; template to match
        %should contain fields: x, y, rate (= rate(x,y) ), rate_background
        %(optional)
        template
        
        %extracted from VolumeReconstructor
        photons = zeros([3 0]); %[x;y;t]       
        
        nostretch = true;
        noshear = false;

        rfc;
        
        tx;
        avec;
    end
    
    methods
        
        function rf = 
        
        function [rf, rffit, success] = registerVR(rf, vr, timerange, temporalchunksize, beamnum, color, czrange) 
            existsAndDefault('beamnum', 2);
            existsAndDefault('color', 'g');
            existsAndDefault('czrange', [-0.5 0.5]);
            alignToFrames = false;
            numchunks = ceil(diff(timerange)/temporalchunksize);
            timeedges = linspace(timerange(1), timerange(end), numchunks+1);
            
            if (alignToFrames)
                timeedges = interp1(vr.timer.frameEdges, vr.timer.frameEdges, timeedges, 'nearest');
            end
            
            if (isempty(rf.template))
                [t, ax] = vr.getXYTemplate(timerange, [], [], [], beamnum, color);
                rf.template.rate = t;
                rf.template.x = ax.x;
                rf.template.y = ax.y;
                rf.template.xe = ax.xe;
                rf.template.ye = ax.ye;
            end
            rf2 = rf.stripMemoryHogs();
            rf2.template = rf.template;
            rf2.verbose = min(rf2.verbose, 1);
            rffit = repmat(rf2, [1 numchunks]);
            t2 = timeedges(2:end);
            ts1 = tic;
            rf.display = 'off';
           
            success = false(size(rffit));
            
            if (~isempty(vr.beam)) %vr already has photon data loaded in, so can't be passed to parallel loops
                for j = 1:numchunks
                    rffit(j) = rffit(j).importVolumeReconstructor(vr, beamnum, color, [timeedges(j) t2(j)], czrange);
                    rffit(j).prefix = sprintf('loop %d: ', j);
                end
                disp('imported vr data');
                parfor j = 1:numchunks   
                    ts2 = tic;
                    fprintf(1, 'loop %d: starting, total time = %.0fs\n', j,  toc(ts1));
                    try
                        rf2 = rffit(j).setParamSteps();
                        rf2 = rf2.configureProblem();
                        rf2 = rf2.fitRasterCascade();
                        rffit(j) = rf2.stripMemoryHogs();
                        fprintf(1, 'loop %d: loop time = %.0fs, total time = %.0fs\n', j, toc(ts2), toc(ts1));
                        success(j) = true;
                    catch me
                        fprintf(2, 'loop %d failed with error\n%s', j, me.getReport());
                    end
                end
            else
                dname = sprintf('rasterFit %d %d %s', round(min(timerange)), round(max(timerange)), datestr(now,30));
                fdir = fullfile(fileparts(vr.filename), dname);
                mkdir(fdir);
                fprintf(1, 'total number of chunks = %d\n', numchunks);
                %todo implement waitbar for parfor by yun pu on file
                %exchange
                parfor j = 1:numchunks   
                    try
                        ts2 = tic;
                        fprintf(1, 'loop %d: starting, total time = %.0fs\n', j,  toc(ts1));
                        rf2 = rffit(j).importVolumeReconstructor(vr, beamnum, color, [timeedges(j) t2(j)], czrange);
                        rf2.prefix = sprintf('loop %d: ', j);
                        rf2 = rf2.setParamSteps();
                        rf2 = rf2.removeTurningPointsFromScan(0.1);
                        rf2 = rf2.configureProblem();
                        rf2 = rf2.fitRasterCascade();
                        rffit(j) = rf2.stripMemoryHogs();
                        fprintf(1, 'loop %d: loop time = %.0fs, total time = %.0fs\n', j, toc(ts2), toc(ts1));
                        saveMe(rffit(j), fullfile(fdir, ['rasterfitter_', num2str(j,3), '.mat']));
                        success(j) = true;
                    catch me
                        fprintf(2, 'loop %d failed with error\n%s', j, me.getReport());
                    end
                end
            end
            ot = cat(2,rffit.offsettime);
            fo = cat(2,rffit.fitoffset);
            [rf.offsettime,~,ic] = unique(ot);
            rf.fitoffset = zeros(2,length(rf.offsettime));
            rf.fitoffset(1,:) = accumarray(ic, fo(1,:), [], @mean)';
            rf.fitoffset(2,:) = accumarray(ic, fo(2,:), [], @mean)';
            
            it = cat(2,rffit.intensitytime);
            fi = cat(2,rffit.fitintensity);
            [rf.intensitytime,~,ic] = unique(it);
            rf.fitintensity = accumarray(ic, fi, [], @mean)';
            
            rt = cat(2,rffit.rotationtime);
            fr = cat(2,rffit.fitrotation);
            [rf.rotationtime,~,ic] = unique(rt);
            rf.fitrotation = accumarray(ic, fr, [], @polarmean)';
           
            
        end
        function saveMe(rf, fname)
            save(fname, 'rf');
        end
        
        
        
        function v = negLogP_photons(rf, rxy)
            %just the value of the log of the rate function for the
            %adjusted photon locations; if rxy is just 1 number, only
            %rotate, if 3 numbers, rotate and translate
            if(size(rxy,2) > 1)
                v = zeros([1 size(rxy,2)]);
                for j = 1:size(v,2)
                    v(j) = negLogP_photons(rf, rxy(:,j));
                end
                return;
            end
            
            if size(rxy,1) >= 3
                offset = rxy(2:3);
            else
                offset = [0;0];
            end
            %function [val, grads, hess] = rasterLineThroughImageValuesAndDerivatives(imdata, pts, offset, theta, scaling, islog,dwelltime)
            v = -RasterFitter.rasterLineThroughImageValuesAndDerivatives(rf.rateim, rf.photons(1:2,:), offset, rxy(1), 1, true);
            
        end
        
        function [v,grad,hess, checks] = negLogP(rf, params)
            [vdiff, gdiff] = rf.continuityTerm(params); %cheap, no need to check nargout
            if (nargout > 2)
                [vphoton, gphoton, hphoton] = rf.rasterLineThroughImageParameterizedValuesAndDerivatives(rf.logim, rf.photons(1:2,:), params, rf.weights_photon, true, [], rf.useintensity);
                if (rf.incldwell)
                    [vdwell, gdwell, hdwell] = rf.rasterLineThroughImageParameterizedValuesAndDerivatives(rf.rateim, rf.dwell(1:2,:), params, rf.weights_dwell, false, rf.dwelltime, rf.useintensity);
                else
                    vdwell = 0;
                    gdwell = 0;
                    hdwell = 0;
                end
                v = vdwell - vphoton - vdiff;
                grad = gdwell - gphoton - gdiff;
                hess = hdwell - hphoton - rf.hess_diffusion;
                if (nargout > 3)
                    checks.vp = isfinite(vphoton);
                    checks.gp = all(isfinite(gphoton));
                    checks.hp = all(isfinite(hphoton),[1 2]);
                    checks.vdw = isfinite(vdwell);
                    checks.gdw = all(isfinite(gdwell));
                    checks.hdw = all(isfinite(hdwell),[1 2]);
                    checks.vdf = isfinite(vdiff);
                    checks.gdf = all(isfinite(gdiff));
                    checks.hdf = all(isfinite(rf.hess_diffusion),[1 2]);
                    checks.allok = isfinite(v) && all(isfinite(grad)) && all(isfinite(hess),[1 2]);
                end
                return;
            end
            if (nargout > 1)
                [vphoton, gphoton] = rf.rasterLineThroughImageParameterizedValuesAndDerivatives(rf.logim, rf.photons(1:2,:), params, rf.weights_photon, true, [], rf.useintensity);
                 if (rf.incldwell)
                    [vdwell, gdwell] = rf.rasterLineThroughImageParameterizedValuesAndDerivatives(rf.rateim, rf.dwell(1:2,:), params, rf.weights_dwell, false, rf.dwelltime, rf.useintensity);
                else
                    vdwell = 0;
                    gdwell = 0;
                    
                end
                v = vdwell - vphoton - vdiff;
                grad = gdwell - gphoton - gdiff;
                return;
            end
            
            vphoton = rf.rasterLineThroughImageParameterizedValuesAndDerivatives(rf.logim, rf.photons(1:2,:), params, rf.weights_photon, true, [], rf.useintensity);
             if (rf.incldwell)
                    vdwell = rf.rasterLineThroughImageParameterizedValuesAndDerivatives(rf.rateim, rf.dwell(1:2,:), params, rf.weights_dwell, false, rf.dwelltime, rf.useintensity);
                else
                    vdwell = 0;
                    
             end
             v = vdwell - vphoton - vdiff;
            return;
        end
        
        
        function [v,grad] = continuityTerm(rf, params)
           
            grad = zeros([rf.Mo*2 + rf.Mr + rf.Mi, 1]);
            
            [vx,grad(rf.inds_xo)] = rf.logDiffProbability(params(rf.inds_xo), rf.v_offset);
            [vy,grad(rf.inds_yo)] = rf.logDiffProbability(params(rf.inds_yo), rf.v_offset);
            [vr,grad(rf.inds_r)] = rf.logDiffProbability(params(rf.inds_r), rf.v_rotation);
            if (rf.incldwell && rf.useintensity)
                [vi,grad(rf.inds_i)] = rf.logDiffProbability(params(rf.inds_i), rf.v_intensity);
            else
                vi = 0;
            end
            v = vx+vy+vr+vi;
            
        end
        
        function rf = importVolumeReconstructor(rf, vr, beamnum, color, timerange)
            photarr = vr.getPhotons(beamnum, color, timerange);
            rf = rf.importPhotArr(photarr);
            if (rf.incldwell)
                dwellarr = vr.getDwell(beamnum, timerange, color);
                rf = rf.importDwellArr(dwellarr);
            end
        end
        
        function rf = importDataFromFile(rf, photfid, dwellfid, append)
            existsAndDefault('append', false);
            if (fseek(photfid, 0, 'bof') < 0)
                warning ('invalid photon fid');
                return;
            end
            photarr = fread(photfid,[6 Inf], 'double');
            rf = rf.importPhotArr(photarr, append);
            if (nargin < 3 || isempty(dwellfid))
                return;
            end
            if (fseek(dwellfid, 0, 'bof') < 0)
                warning ('invalid dwell fid');
                return;
            end
            dwellarr = fread(dwellfid,[6 Inf], 'double');
            rf = rf.importDwellArr(dwellarr, append);
        end
        function rf = importPhotArr(rf, photarr, append)
             %photarr = [time;phase;pos;zo];
             existsAndDefault('append', false);
             z = cos(photarr(2,:)) - photarr(6,:);
             valid = z >= rf.czrange(1) & z <= rf.czrange(2);
             if append
                 rf.photons = [rf.photons photarr([3 4 1], valid)];
                 [~,I] = sort(rf.photons(3,:), 'ascend');
                 rf.photons = rf.photons(:,I);
             else
                 rf.photons = photarr([3 4 1], valid);
             end
        end
        function rf = importDwellArr(rf, dwellarr, append)
            % dwellarr = [time;period;pos;zo];   
            existsAndDefault('append', false);
            if append
                 rf.dwell = [rf.dwell dwellarr([3 4 1], valid)];
                 czr = min([1 1], max(rf.czrange, [-1 -1]));
                 rf.dwelltime = [rf.dwelltime dwellarr(2,:)/2 * (abs(diff(acos(czr)))/(2 * pi))];
                 [~,I] = sort(rf.dwell(3,:), 'ascend');
                 rf.dwell = rf.dwell(:,I);
                 rf.dwelltime = rf.dwelltime(I);
            else
                rf.dwell = dwellarr([3 4 1],:);
                czr = min([1 1], max(rf.czrange, [-1 -1]));
                rf.dwelltime = dwellarr(2,:)/2 * (abs(diff(acos(czr)))/(2 * pi));
            end
        end
        
        
        function rf = setParamSteps (rf, no, nr, ni)
            if (rf.incldwell && ~isempty(rf.dwell))
                trange = [min(rf.photons(3,1), rf.dwell(3,1)),max(rf.photons(3,end), rf.dwell(3,end))];
            else
                trange = rf.photons(3,[1 end]);
            end
            existsAndDefault('no', max(2, ceil(diff(trange)/rf.offset_dt)));
            existsAndDefault('nr', max(2, ceil(diff(trange)/rf.rotation_dt)));
            if (rf.incldwell && rf.useintensity)
                existsAndDefault('ni', max(2, ceil(diff(trange)/rf.intensity_dt)));
            else
                ni = 0;
            end
            rf.offsettime = linspace(min(trange), max(trange), max(2,no));
            rf.rotationtime = linspace(min(trange), max(trange), max(2,nr));
            rf.intensitytime = linspace(min(trange), max(trange), max(2,ni));

            rf.Mo = numel(rf.offsettime);
            rf.Mr = numel(rf.rotationtime);
            rf.Mi = numel(rf.intensitytime);
        end

        function rf = downsample(rf, frac)
            assert(frac > 0 && frac <= 1);
            if (rf.incldwell)
                numdwell = frac*size(rf.dwell,2);
                numdwell = round(min(size(rf.dwell,2), numdwell));
                dvalid = zeros([1,size(rf.dwell,2)]);
                di = randperm(size(rf.dwell,2), numdwell);
                dvalid(di) = 1;
                pvalid = logical(interp1(rf.dwell(3,:), dvalid, rf.photons(3,:), 'nearest', 0));
                dvalid = logical(dvalid);
                rf.dwell = rf.dwell(:,dvalid);
                rf.dwelltime = rf.dwelltime(dvalid);
                rf.photons = rf.photons(:,pvalid);
            else
                numphot = ceil(size(rf.photons, 2)*frac);
                pvalid = false([1, size(rf.photons,2)]);
                pvalid(randperm(size(rf.photons,2),numphot)) = true;
                rf.photons = rf.photons(:,pvalid);
                
            end
            rf = rf.configureProblem();
        end
%         
%         function rf = removeTurningPointsFromScan (rf, edgefrac)
%             existsAndDefault('edgefrac', 0.1);
%            
%             vv = true(size(rf.dwell(3,:)));
%             sigma = 5e-5/median(diff(rf.dwell(3,:))); %
%             for j = rf.stripTpAxes
%                 a0 = std(rf.dwell(j,:))/.5774;
%                 ci = find(diff(sign(deriv(rf.dwell(j,:), sigma))));
%                 f0 = 0.5*nnz(ci)./diff(rf.dwell(3,[1 end]));               
%                 ft = fftshift(fft(rf.dwell(j,:)));
%                 f = linspace(-0.5,0.5,length(ft))./median(diff(rf.dwell(3,:)));
%                 p0 = pi/2+angle(interp1(f,ft,f0));
%                 myfun = @(x,xdata) sawtooth(2*pi*(x(1)+f0)*xdata + x(2) , 0.5);
%                 
%                 x0 = [0 p0];
%                 if (rf.verbose > 1)
%                     fprintf(1, '%sinitial guess axis %d: f = %.5f Hz, a = %.2f microns, phi = %.2f rad\n',rf.prefix, j, f0 , a0, p0);
%                 end
%                 
%                 op = optimoptions(@lsqcurvefit);
%                 op.Display = 'off';
%                 op.Algorithm = 'levenberg-marquardt';
%                 lb = [];
%                 ub = [];
%              
%                 xfit = lsqcurvefit(myfun, [0 0], rf.dwell(3,1:10:end), (rf.dwell(j,1:10:end)-mean(rf.dwell(j,1:10:end)))/(a0), lb,ub,op);
%                
%                 xfit(2) = atan2(sin(xfit(2)), cos(xfit(2)));
%                 
%                 vv = vv & (abs(myfun(xfit, rf.dwell(3,:))) < (1-edgefrac));
%                 if (rf.verbose)
%                     fprintf(1, '%saxis %d: f = %.5f Hz, a = %.2f microns, phi = %.2f rad\n',rf.prefix, j, xfit(1) + f0,a0, xfit(2));
%                 end
%                 if (rf.verbose > 2)
%                     figure(1);
%                     subplot(2,2,j); plot ( rf.dwell(3,1:10:end),rf.dwell(j,1:10:end),  rf.dwell(3,1:10:end), a0*myfun(xfit, rf.dwell(3,1:10:end)), rf.dwell(3,1:10:end), a0*myfun(x0, rf.dwell(3,1:10:end)));
%                     legend('data', 'fit', 'initial guess');
%                     xlim(mean(rf.dwell(3,:)) + [-0.002 0.002]);
%                     subplot(2,2,j+2); plot ( rf.dwell(3,1:10:end),rf.dwell(j,1:10:end)-a0*myfun(x0, rf.dwell(3,1:10:end)),rf.dwell(3,1:10:end),rf.dwell(j,1:10:end)-a0*myfun(xfit, rf.dwell(3,1:10:end))); title ('scan - fit');
%                     legend('initial', 'fit');
%              
%                 end
%             end
%             pvalid = logical(interp1(rf.dwell(3,:), double(vv), rf.photons(3,:), 'nearest', 'extrap'));
%             rf.dwell = rf.dwell(:,vv);
%             rf.dwelltime = rf.dwelltime(:,vv);
%             rf.photons = rf.photons(:,pvalid);
%             
%             
%         end
        
        function rf = configureProblem (rf)
           
            rf.useintensity = rf.useintensity && rf.incldwell;
            rf.weights_photon = rf.getWeightMatrix(rf.photons(3,:), rf.offsettime, rf.rotationtime, rf.intensitytime, rf.useintensity);
            if (rf.incldwell)
                rf.weights_dwell =  rf.getWeightMatrix(rf.dwell(3,:), rf.offsettime, rf.rotationtime, rf.intensitytime, rf.useintensity);
            end
            if (~isfield(rf.template, 'x') && isfield(rf.template, 'u'))
                rf.template.x = rf.template.u;
            end
            if (~isfield(rf.template, 'y') && isfield(rf.template, 'v'))
                rf.template.y = rf.template.v;
            end
            if (~isfield(rf.template, 'xe'))
                x = rf.template.x;
                x = [2*x(1)-x(2) x 2*x(end)-x(end-1)];
                rf.template.xe = conv(x,[0.5 0.5],  'valid');
            end
             if (~isfield(rf.template, 'ye'))
                y = rf.template.y;
                y = [2*y(1)-y(2) y 2*y(end)-y(end-1)];
                rf.template.ye = conv(y,[0.5 0.5],  'valid');
            end
            rf.rateim = struct('u', rf.template.x, 'v', rf.template.y);
            rf.logim = rf.rateim;
            F = rf.template.rate;
            if (~isfield(rf.template, 'rate_background'))
                Fextrap = min(F(F>0));
            else
                Fextrap = rf.template.rate_background;
            end
            F(F < Fextrap) = Fextrap;
            rf.rateim.F = F;
            rf.rateim.Fextrap = Fextrap;
            rf.logim.F = log(F);
            rf.logim.Fextrap = log(Fextrap);
            
           % if (rf.incldwell)
                rf.rateim = rf.createImModel(rf.rateim);
            %end
            rf.logim = rf.createImModel(rf.logim);
            
            rf.v_offset = diff(rf.offsettime)*rf.Doffset;
            rf.v_rotation = diff(rf.rotationtime)*rf.Drotation;
            rf.v_intensity = diff(rf.intensitytime)*rf.Dintensity;
            
            rf.Mo = numel(rf.offsettime);
            rf.Mr = numel(rf.rotationtime);
            if (rf.useintensity)
                rf.Mi = numel(rf.intensitytime);
                rf.inds_i = ((2*rf.Mo+rf.Mr+1):(2*rf.Mo+rf.Mr+rf.Mi));
                [ii,ji,hi] = rf.logDiffProgbabilityHessian(rf.v_intensity);
            else
                rf.Mi = 0;
                rf.inds_i = [];
            end
            rf.inds_xo = 1:rf.Mo;
            rf.inds_yo = ((rf.Mo+1):(2*rf.Mo));
            rf.inds_r = ((2*rf.Mo+1):(2*rf.Mo+rf.Mr));
            
            [io,jo,ho] = rf.logDiffProgbabilityHessian(rf.v_offset);
            [ir,jr,hr] = rf.logDiffProgbabilityHessian(rf.v_rotation);

            if (rf.useintensity)
                i = [io io+rf.Mo ir+2*rf.Mo ii+2*rf.Mo+rf.Mr];
                j = [jo jo+rf.Mo jr+2*rf.Mo ji+2*rf.Mo+rf.Mr];
                h = [ho ho hr hi];
            else
                i = [io io+rf.Mo ir+2*rf.Mo];
                j = [jo jo+rf.Mo jr+2*rf.Mo];
                h = [ho ho hr];
            end
            rf.hess_diffusion = sparse(i,j,h,2*rf.Mo+rf.Mr+rf.Mi, 2*rf.Mo+rf.Mr+rf.Mi);
            
        end

        
        
        function [params_init, initrotation, initoffset, initintensity] = getInitialParameters(rf, initrotation, initoffset,initintensity)
            blurmicron = 5;
            edgebuffer = min(diff(rf.template.x([1 end])), diff(rf.template.y([1 end])))/8; %central 75%
            if (nargin == 1 || isempty(initrotation))                
                [initrotation, initoffset] = rf.initialRotationAndOffset(blurmicron,edgebuffer);
            end
            if (nargin == 2)
                rf2 = rf.initialAlignmentRasterFitter(blurmicron, edgebuffer).fitRaster(initrotation);
                initrotation = polarmean(rf2.fitrotation);
                initoffset = mean(rf2.fitoffset,2);
            end
            if (numel(initoffset) ~= 2*rf.Mo)        
                ix = repmat(initoffset(1,1),[rf.Mo 1]);
                iy = repmat(initoffset(2,1),[rf.Mo 1]);
            else
                ix = initoffset(1,:); 
                iy = initoffset(2,:);
            end
            if (numel(initrotation) ~= rf.Mr)
                ir = repmat(initrotation, [rf.Mr 1]);
            else
                ir = initrotation;
            end
            if (nargin < 4)
                initintensity = 1;
            end
            if (numel(initintensity) ~= rf.Mi)
                ii = repmat(initintensity, [rf.Mi 1]);
            else
                ii = initintensity;
            end
            if (isrow(ix)), ix = ix'; end
            if (isrow(iy)), iy = iy'; end
            if (isrow(ir)), ir = ir'; end
            if (isrow(ii)), ii = ii'; end
           % ii = ones([rf.Mi 1]);
            params_init = [ix;iy;ir;ii];
        end

        function [params, offset, rotation, intensity] = adaptParams(rfto, rffrom, params)
            %function [offset, rotation, intensity] = adaptParams(rfto, rffrom, params)
            if nargin < 3
                params = rffrom.fitparams;
            end
            o = params([rffrom.inds_xo;rffrom.inds_yo]);
            r = params(rffrom.inds_r);
            ii = params(rffrom.inds_i);
            offset = interp1(rffrom.offsettime, o', rfto.offsettime, 'linear', 'extrap')';
            rotation = interp1(rffrom.rotationtime, r, rfto.rotationtime, 'linear', 'extrap');
            if (~isempty(ii))
                intensity = interp1(rffrom.intensitytime, ii, rfto.intensitytime, 'linear', 'extrap');
            else
                intensity = [];
            end
            params = [offset(1,:)';offset(2,:)';rotation';intensity'];
        end

%         function rf = fitRasterCascade (rf, varargin)
%             % rf = fitRasterCascade(rf, initrotation, initoffset, initintensity)
%             
%             %just make sure we're not fitting outside the range by accident
%     %        trange = [min([rf.photons(3,:) rf.dwell(3,:)]) max([rf.photons(3,:) rf.dwell(3,:)])];
%             trange = rf.photons(3,[1 end]);
%             rf.offsettime = rf.offsettime(rf.offsettime >= trange(1) & rf.offsettime <= trange(2));
%             rf.intensitytime = rf.intensitytime(rf.intensitytime >= trange(1) & rf.intensitytime <= trange(2));
%             rf.rotationtime = rf.rotationtime(rf.rotationtime >= trange(1) & rf.rotationtime <= trange(2));
%             
%             if (length(rf.offsettime) < 2 || length(rf.rotationtime) < 2)
%                 error ('offset and/or rotation times don''t match data -- run setParamSteps()');
%             end
%             
%                 
%             
%             rf = rf.configureProblem();
%             params = getInitialParameters(rf, varargin{:});
%             maxedge = min(diff(rf.template.x([1 end])), diff(rf.template.y([1 end])))/6;
%             
%             ts1 = tic;
%             
%             
% %             [u,v] = rf.adjustedCoords(rf.dwell(1:2,:), params, rf.weights_dwell);
% %             edge0 = max([max(rf.template.x(1)-u), max(u-rf.template.x(end)), max(rf.template.y(1)-v), max(v-rf.template.y(end))]);
% %             
%             nsteps = ceil(max(log2(max([rf.Mo, rf.Mi, rf.Mr]))));
%             no = round(2.^(linspace(1, log2(rf.Mo), nsteps)));
%             nr = round(2.^(linspace(1, log2(rf.Mr), nsteps)));
%            % numdwell = round(2.^(linspace(log2(5e4), log2(size(rf.dwell,2)),nsteps)));
%             blur = 2.^(linspace(log2(5), -1, nsteps));
%             if (rf.Mi > 0)
%                 ni = round(2.^(linspace(1, log2(rf.Mi), nsteps)));
%             else
%                 ni = 2*ones(size(no));
%             end
%             if (rf.verbose > 0)
%                 fprintf(1, '%sincreasing parameter resolution while decreasing blur, total steps = %d\n',rf.prefix, nsteps);
%             end
%             for j = 1:nsteps
%                 if (rf.verbose > 1)
%                     fprintf(1, '%sincreasing parameter resolution while decreasing blur, step %d of %d\n',rf.prefix, j, nsteps);
%                 end
%                 %edge = min(edge0+blur(j), maxedge);
%                 rf2 = rf.blurAndCrop(blur(j), maxedge);
%              %   rf2 = rfbc.downsample(numdwell(j));
%                 rf2 = rf2.setParamSteps(no(j), nr(j), ni(j));
%                 rf2.Doffset = rf.Doffset * rf.Mo/rf2.Mo;
%                 rf2.Drotation = rf.Drotation * rf.Mr/rf2.Mr;
%                 if (rf2.Mi > 0)
%                     rf2.Dintensity = rf.Dintensity * rf.Mi/rf2.Mi;
%                 end
%                 rf2 = rf2.configureProblem();
%                 [~,offset, rotation, intensity] = rf2.adaptParams(rf, params);
%                 if (rf.verbose > 1)
%                     fprintf(1,'%sinitial mean offset = %g,%g ; initial mean rotation = %.1f\n',rf.prefix, mean(offset,2), polarmean(rotation));
%                 end
%                 rf2 = rf2.fitRaster(rotation, offset, intensity);
%                 if (rf.verbose > 1)
%                     fprintf(1,'%sfit mean offset = %g,%g ; fit mean rotation = %.1f\n',rf.prefix, mean(rf2.fitoffset,2), polarmean(rf2.fitrotation));
%                 end
% %                 if (rf.verbose > 2)
% %                     figure(1); clf(1);
% %                     subplot(2,2,1); plot(rf2.offsettime, offset(1,:),'b.-', rf2.offsettime, rf2.fitoffset(1,:),'r.-'); title ('x-offset'); legend('init', 'fit')
% %                     subplot(2,2,2); plot(rf2.offsettime, offset(2,:),'b.-', rf2.offsettime, rf2.fitoffset(2,:),'r.-'); title ('y-offset'); legend('init', 'fit')
% %                     subplot(2,2,3); plot(rf2.rotationtime, rotation*180/pi, 'b.-',rf2.rotationtime, rf2.fitrotation*180/pi,'r.-'); title ('rotation'); legend('init', 'fit')
% %                     subplot(2,2,4); plot(rf2.intensitytime, intensity, 'b.-',rf2.intensitytime, rf2.fitintensity,'r.-'); title ('intensity'); legend('init', 'fit')
% %                     pause(0.01);
% %                 end
%                 params = rf.adaptParams(rf2);
%                 
%             end
% 
%             if (rf.verbose > 0)
%                 fprintf(1, '%sfinished parameter increase and blur decrease: elapsed time = %.0fs\n', rf.prefix, toc(ts1));
%             end
% 
%             rf = rf.fitRaster(params(rf.inds_r), params([rf.inds_xo;rf.inds_yo]), params(rf.inds_i));
%             if (rf.verbose > 0)
%                 fprintf(1, '%sfinished fit: elapsed time = %.0fs\n', rf.prefix, toc(ts1));
%                 
%             end
% 
%         end
%         
        function rf = fitRasterCascade (rf, varargin)
            % rf = fitRasterCascade(rf, initrotation, initoffset, initintensity)
            
            %just make sure we're not fitting outside the range by accident
    %        trange = [min([rf.photons(3,:) rf.dwell(3,:)]) max([rf.photons(3,:) rf.dwell(3,:)])];
            trange = rf.photons(3,[1 end]);
            rf.offsettime = rf.offsettime(rf.offsettime >= trange(1) & rf.offsettime <= trange(2));
            rf.intensitytime = rf.intensitytime(rf.intensitytime >= trange(1) & rf.intensitytime <= trange(2));
            rf.rotationtime = rf.rotationtime(rf.rotationtime >= trange(1) & rf.rotationtime <= trange(2));
            
            if (length(rf.offsettime) < 2 || length(rf.rotationtime) < 2)
                error ('offset and/or rotation times don''t match data -- run setParamSteps()');
            end
            
                
            
            rf = rf.configureProblem();
            params = getInitialParameters(rf, varargin{:});
            maxedge = min(diff(rf.template.x([1 end])), diff(rf.template.y([1 end])))/6;
            
            ts1 = tic;
            
            
            nsteps = max(3,ceil(0.5*max(log2(max([rf.Mo, rf.Mi, rf.Mr])))));
            no = round(2.^(linspace(1, log2(rf.Mo), nsteps)));
            nr = round(2.^(linspace(1, log2(rf.Mr), nsteps)));
            blur = 2.^(linspace(log2(5), -1, nsteps));
            fmin = min(1,5e4/length(rf.photons));
            frac = 2.^(linspace(log2(fmin), 0, nsteps));
            if (rf.Mi > 0)
                ni = round(2.^(linspace(1, log2(rf.Mi), nsteps)));
            else
                ni = 2*ones(size(no));
            end
            if (rf.verbose > 0)
                fprintf(1, '%sincreasing parameter resolution while decreasing blur, total steps = %d\n',rf.prefix, nsteps);
            end
            for j = 1:nsteps
                if (rf.verbose > 1)
                    fprintf(1, '%sincreasing parameter resolution while decreasing blur, step %d of %d\n',rf.prefix, j, nsteps);
                end
                
                rf2 = rf.blurAndCrop(blur(j),  min(3*blur(j), maxedge)).downsample(frac(j));
             %   rf2 = rfbc.downsample(numdwell(j));
                rf2 = rf2.setParamSteps(no(j), nr(j), ni(j));
                rf2.Doffset = rf.Doffset * rf.Mo/rf2.Mo;
                rf2.Drotation = rf.Drotation * rf.Mr/rf2.Mr;
                if (rf2.Mi > 0)
                    rf2.Dintensity = rf.Dintensity * rf.Mi/rf2.Mi;
                end
                rf2 = rf2.configureProblem();
                [~,offset, rotation, intensity] = rf2.adaptParams(rf, params);
                if (rf.verbose > 1)
                    fprintf(1,'%sinitial mean offset = %g,%g ; initial mean rotation = %.1f\n',rf.prefix, mean(offset,2), polarmean(rotation));
                end
                rf2 = rf2.fitRaster(rotation, offset, intensity);
                if (rf.verbose > 1)
                    fprintf(1,'%sfit mean offset = %g,%g ; fit mean rotation = %.1f\n',rf.prefix, mean(rf2.fitoffset,2), polarmean(rf2.fitrotation));
                end
%                 if (rf.verbose > 2)
%                     figure(1); clf(1);
%                     subplot(2,2,1); plot(rf2.offsettime, offset(1,:),'b.-', rf2.offsettime, rf2.fitoffset(1,:),'r.-'); title ('x-offset'); legend('init', 'fit')
%                     subplot(2,2,2); plot(rf2.offsettime, offset(2,:),'b.-', rf2.offsettime, rf2.fitoffset(2,:),'r.-'); title ('y-offset'); legend('init', 'fit')
%                     subplot(2,2,3); plot(rf2.rotationtime, rotation*180/pi, 'b.-',rf2.rotationtime, rf2.fitrotation*180/pi,'r.-'); title ('rotation'); legend('init', 'fit')
%                     subplot(2,2,4); plot(rf2.intensitytime, intensity, 'b.-',rf2.intensitytime, rf2.fitintensity,'r.-'); title ('intensity'); legend('init', 'fit')
%                     pause(0.01);
%                 end
                params = rf.adaptParams(rf2);
                
            end

            if (rf.verbose > 0)
                fprintf(1, '%sfinished parameter increase and blur decrease: elapsed time = %.0fs\n', rf.prefix, toc(ts1));
            end

            rf = rf.fitRaster(params(rf.inds_r), params([rf.inds_xo;rf.inds_yo]), params(rf.inds_i));
            if (rf.verbose > 0)
                fprintf(1, '%sfinished fit: elapsed time = %.0fs\n', rf.prefix, toc(ts1));
                
            end

        end
        
       
        
        function rf = fitRaster(rf, initrotation, initoffset, initintensity)
            blurmicron = 5;
            if (nargin > 1 && ischar(initrotation) && strcmpi(initrotation, 'usefit'))
                initrotation = rf.fitrotation;
                initoffset = rf.fitoffset;
                initintensity = rf.fitintensity;
            else

                edgebuffer = min(diff(rf.template.x([1 end])), diff(rf.template.y([1 end])))/8; %central 75%
                if (nargin == 1 || isempty(initrotation))                
                    [initrotation, initoffset] = rf.initialRotationAndOffset(blurmicron,edgebuffer);
                end
                if (nargin == 2)
                    rf2 = rf.initialAlignmentRasterFitter(blurmicron, edgebuffer).fitRaster(initrotation,[0;0]);
                    initrotation = polarmean(rf2.fitrotation);
                    initoffset = mean(rf2.fitoffset,2);
                end
            end
            if (numel(initoffset) ~= 2*rf.Mo)        
                ix = repmat(initoffset(1,1),[rf.Mo 1]);
                iy = repmat(initoffset(2,1),[rf.Mo 1]);
            else
                ix = initoffset(1,:); 
                iy = initoffset(2,:);
            end
            if (numel(initrotation) ~= rf.Mr)
                ir = repmat(initrotation, [rf.Mr 1]);
            else
                ir = initrotation;
            end
            if (nargin < 4)
                initintensity = 1;
            end
            if (numel(initintensity) ~= rf.Mi)
                ii = repmat(initintensity, [rf.Mi 1]);
            else
                ii = initintensity;
            end
            if (isrow(ix)), ix = ix'; end
            if (isrow(iy)), iy = iy'; end
            if (isrow(ir)), ir = ir'; end
            if (isrow(ii)), ii = ii'; end
           % ii = ones([rf.Mi 1]);
            params_init = [ix;iy;ir;ii];
            objfun = @(x) rf.negLogP(x);
            [vinit, ginit, hinit, checks] = rf.negLogP(params_init);
            if (~checks.allok)
                checks
                error ('initial conditions produce nonfinite values');
            end
            if (2*rf.Mo + rf.Mr + rf.Mi) > 20
                options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective', 'Display', rf.display);
             else
                 options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true, 'Display', rf.display);
             end
            
            rf.fitparams = fminunc(objfun, params_init, options);
            rf.fitLogP = -objfun(rf.fitparams);
            rf.fitoffset = rf.fitparams([rf.inds_xo;rf.inds_yo]);
            rf.fitrotation = rf.fitparams(rf.inds_r)';
            rf.fitintensity = rf.fitparams(rf.inds_i)';
            if (rf.verbose > 2)             
                figure(2); clf(2);
                subplot(2,2,1); plot(rf.offsettime, ix,'b.-', rf.offsettime, rf.fitoffset(1,:),'r.-'); title ('x-offset'); ylabel('um'); legend('init', 'fit')
                subplot(2,2,2); plot( rf.offsettime, iy,'b.-', rf.offsettime, rf.fitoffset(2,:),'r.-'); title ('y-offset'); ylabel('um');legend('init', 'fit')
                subplot(2,2,3); plot(rf.rotationtime,  ir*180/pi, 'b.-',rf.rotationtime, rf.fitrotation*180/pi,'r.-'); ylabel ('degrees'); title ('rotation'); legend('init', 'fit')
                if (~isempty(ii))
                    subplot(2,2,4); plot(rf.intensitytime, ii, 'b.-',rf.intensitytime, rf.fitintensity,'r.-'); title ('intensity'); legend('init', 'fit')             
                end
            end
        end
        function [initrotation, initoffset] = initialRotationAndOffset(rf, blurmicrons, edgemicrons)
            existsAndDefault('numdwell', length(rf.dwell));
            rf2 = rf.initialAlignmentRasterFitter(blurmicrons, edgemicrons);
           
            theta = linspace(-3*pi/4, 3*pi/4, 4);
            thetafinal = theta;
            offset = zeros(2,length(theta));
            logP = theta;
            for j = 1:length(theta)
                rf2 = rf2.fitRaster(theta(j),[0;0]);
                thetafinal(j) = polarmean(rf2.fitrotation);
                logP(j) = rf2.fitLogP;
                offset(:,j) = mean(rf2.fitoffset, 2);
                if (rf.verbose > 1)
                    fprintf(1, '%sinitial theta -> %.0f ; final theta -> %.0f ; logP -> %g\n',rf.prefix, theta(j)*180/pi, thetafinal(j)*180/pi, logP(j)');
                end
            end
            
            logP = logP - min(logP);
            if (rf.verbose > 2)
                figure(1); 
                subplot(2,2,1); plotyy (theta, thetafinal, theta, logP); xlabel('theta initial');
                subplot(2,2,2); plot(thetafinal, logP, 'bo-'); xlabel('theta final'); ylabel('logP');
                subplot(2,2,3); plot(theta, offset); xlabel('theta initial'); ylabel('offset');
                subplot(2,2,4); plot(thetafinal, offset, 'ko-');xlabel('theta final'); ylabel('offset');
                pause(0.01);
            end
            [~,I] = max(logP);
            initrotation = thetafinal(I);
            initoffset = offset(:,I);
            if (rf.verbose > 0)
                fprintf(1, '%sinitial theta -> %.0f deg; initial offset -> %.1f,%.1f\n',rf.prefix, initrotation*180/pi, initoffset);
            end
        end

        function [im, counts, dwell] = correctedIm(rf, timerange, xe, ye, correctionlevel, blursigma)
            existsAndDefault('timerange', [-Inf Inf]);
            existsAndDefault('xe', rf.template.xe);
            existsAndDefault('ye', rf.template.ye);
            existsAndDefault('correctionlevel', 2);
            existsAndDefault('blursigma', 1);
            tvalid = @(x) (x >= min(timerange)) & (x < max(timerange));
            pv = tvalid(rf.photons(3,:));
            dv = tvalid(rf.dwell(3,:));
            if (correctionlevel)
                [up,vp] = rf.adjustedCoords(rf.photons(1:2,:), rf.fitparams, rf.weights_photon);
                [ud,vd] = rf.adjustedCoords(rf.dwell(1:2,:), rf.fitparams, rf.weights_dwell);
                
            else
                up = rf.photons(1,:); vp = rf.photons(2,:);
                ud = rf.dwell(1,:); vd = rf.dwell(2,:);
            end
            if (correctionlevel > 1 && ~isempty(rf.fitintensity))
                ic = interp1(rf.intensitytime, rf.fitintensity, rf.dwell(3,:),'linear',1);
            else
                ic = ones(size(rf.dwelltime));
            end
            up = up(pv); vp = vp(pv); ud = ud(dv); vd = vd(dv);
            dwt = ic(dv)    .*rf.dwelltime(dv);

            [~,~,xp] = histcounts(up, xe);
            [~,~,yp] = histcounts(vp, ye);
            sz = [length(xe)-1 length(ye)-1];
            counts = accumarray([xp(xp&yp)' yp(xp&yp)'], 1, sz);
            [~,~,xd] = histcounts(ud, xe);
            [~,~,yd] = histcounts(vd, ye);
            dwell = accumarray([xd(xd&yd)' yd(xd&yd)'], dwt(xd&yd), sz);
            if (blursigma > 0)
                im = conv2(gaussKernel(blursigma), gaussKernel(blursigma), counts, 'same')./conv2(gaussKernel(blursigma), gaussKernel(blursigma), dwell, 'same');
            else
                im = counts./dwell;
            end
            im(~isfinite(im)) = 0;
        end

        function rf = blurAndCrop(rf, blurmicrons, edgemicrons)
            sigmax = blurmicrons/median(diff(rf.template.x), 'omitnan');
            sigmay = blurmicrons/median(diff(rf.template.y), 'omitnan');
            gkx = gaussKernel(sigmax);
            gky = gaussKernel(sigmay);
            rf.template.rate = conv2(gkx, gky, padarray(rf.template.rate, [floor(length(gkx)/2) floor(length(gky)/2)], 'replicate'), 'valid');
            if (~isempty(rf.photons))
                pv = rf.photons(1,:) > min(rf.template.x) + edgemicrons & rf.photons(1,:) < max(rf.template.x) - edgemicrons & ...
                    rf.photons(2,:) > min(rf.template.y) + edgemicrons & rf.photons(2,:) < max(rf.template.y) - edgemicrons;
                rf.photons = rf.photons(:,pv);
            end
            if (~isempty(rf.dwell))
                dv = rf.dwell(1,:) > min(rf.template.x) + edgemicrons & rf.dwell(1,:) < max(rf.template.x) - edgemicrons & ...
                    rf.dwell(2,:) > min(rf.template.y) + edgemicrons & rf.dwell(2,:) < max(rf.template.y) - edgemicrons;
                rf.dwell = rf.dwell(:,dv);
                rf.dwelltime = rf.dwelltime(dv);
            end
        end
        function rf = initialAlignmentRasterFitter (rf, blurmicrons, edgemicrons)
            frac = min(1, 4e4/length(rf.photons));
            rf = rf.blurAndCrop(blurmicrons,edgemicrons);
            rf = rf.downsample(frac);
            if (rf.incldwell)
                 trange = [min(rf.photons(3,1), rf.dwell(3,1)),max(rf.photons(3,end), rf.dwell(3,end))];
            else
                trange = rf.photons(3,[1 end]);
            end
           
            rf.offsettime = trange;
            rf.rotationtime = trange;
            rf.intensitytime = trange;
            rf.display = 'off';
%             rf.Dintensity = 1e-8;
%             rf.Drotation = 1e-8;
%             rf.Doffset = 1e-8; 
%             rf.useintensity = false;
            rf = rf.configureProblem();
        end
            

    end
    
    
    methods (Static)
        function weight = getWeightMatrix (pttime, offsettime, thetatime, scalingtime, usescaling)
            %weights is 4N x M
            Mo = numel(offsettime);
            Mt = numel(thetatime);
            if (usescaling)
                Ms = numel(scalingtime);
            else
                Ms = 0;
            end
            N = numel(pttime);
            M = 2*Mo+Mt+Ms;
            [mo, no, wo] = RasterFitter.getWeights(pttime, offsettime);
            [mt, nt, wt] = RasterFitter.getWeights(pttime, thetatime);
            if (usescaling)
                [ms, ns, ws] = RasterFitter.getWeights(pttime, scalingtime);
                m = [mo mo+Mo mt+2*Mo ms+2*Mo+Mt];
                n = [no no+N nt+2*N ns+3*N];
                weight = sparse(n,m,[wo wo wt ws], 4*N, M);
            else
                m = [mo mo+Mo mt+2*Mo];
                n = [no no+N nt+2*N];
                weight = sparse(n,m,[wo wo wt], 3*N, M);
            end
            
        end
        function [m,n,w] = getWeights (pttime, paramtime)
            %pttime, paramtime should be 1XN, 1xM respectively
            if (~isrow(pttime))
                pttime = pttime';
            end
            if (~isrow(paramtime))
                paramtime = paramtime';
            end
            prev = interp1(paramtime, 1:length(paramtime), pttime, 'previous', 'extrap');
            next =  interp1(paramtime, 1:length(paramtime), pttime, 'next', 'extrap');
            pv = isfinite(prev);
            nv = isfinite(next) & (next ~= prev); 
            
            mp = prev(pv);
            np = find(pv);
            mn = next(nv);
            nn = find(nv);
            
            wp = zeros(size(prev));
            wn = zeros(size(next));
            wp(~nv) = 1;
            wn(~pv) = 1;
            valid = pv & nv;
            prev = prev(valid); next = next(valid);
            interval = paramtime(next)-paramtime(prev);
            wp(valid) = (paramtime(next)-pttime(valid))./interval;
            wn(valid) = 1 - wp(valid);
            
            m = [mp mn];
            n = [np nn];
            w = [wp(pv) wn(nv)];
            
        end
        
        function [val,grad] = logDiffProbability (x, v)
            dx = diff(x);
            if (isrow(v))
                v = v';
            end
            val = sum(-0.5*dx.^2./v);
            if (nargout < 2)
                return;
            end
            dx = dx./v;
            grad = [dx;0] - [0;dx];
        end
        function [i,j,val] = logDiffProgbabilityHessian(v)
            m = length(v+1);
            inds = 1:m;
            i = [inds inds+1 inds inds+1];
            j = [inds inds+1 inds+1 inds];
            val = [-v -v v v];
        end
        
        function [u,v] = adjustedCoords(pts, paramvec, weights)
             N = size(pts, 2);
            offthetascale = (weights*paramvec)';
            ind_x = 1:N;
            ind_y = (N+1):(2*N);
            ind_t = (2*N+1):3*N;
            
            offset = offthetascale([ind_x;ind_y]);
            theta = offthetascale(ind_t);

            u = cos(theta).*(pts(1,:) + offset(1,:)) -  sin(theta).*(pts(2,:) + offset(2,:));
            v = sin(theta).*(pts(1,:) + offset(1,:)) +  cos(theta).*(pts(2,:) + offset(2,:));
        end
        
        function [funval,grad, hess] = rasterLineThroughImageParameterizedValuesAndDerivatives(imdata, pts, paramvec, weights, islog, dwelltime, useintensity)
            % function [val,grad, hess] = rasterLineThroughImageParameterizedValuesAndDerivatives(pts, paramvec, weights, imdata)
            % calculates values and derivatives with respect to offset and theta given
            % scan along pts through image defined by imdata
            % calls rasterLineThroughImageValuesAndDerivatives and functions similarly,
            % but instead of pts, offset, theta, all being same length
            %
            % pts is 2xN
            % paramvec is M x 1
            % weights is 4N x M
            % weights*paramvec = [pts_x';pts_y';theta';scaling']
            %
            % val is sum(F(u,v).*scaling) or sum(F(u,v) + log(scaling)); u,v given below
            % grad is M x 1; grad(j) = dval / d paramvec(j)
            % hess is M x M; hess(j,k) = d^2val / d paramvec(j) d paramvec (k)
            % hess is a sparse matrix
            %
             existsAndDefault('dwelltime', []);
            N = size(pts, 2);
            offthetascale = (weights*paramvec)';
            ind_x = 1:N;
            ind_y = (N+1):(2*N);
            ind_t = (2*N+1):3*N;
            ind_s = (3*N+1):4*N;
            
            offset = offthetascale([ind_x;ind_y]);
            theta = offthetascale(ind_t);
            if (useintensity)    
                scaling = offthetascale(ind_s);
            else
                scaling = ones(size(theta));
            end
            if (nargout == 1)
                funval =  RasterFitter.rasterLineThroughImageValuesAndDerivatives(imdata, pts, offset, theta, scaling, islog, dwelltime);
                return
            end
            if (nargout == 2)
                [funval, grad] =  RasterFitter.rasterLineThroughImageValuesAndDerivatives(imdata, pts, offset, theta, scaling, islog, dwelltime);
            else
                [funval, grad, h] =  RasterFitter.rasterLineThroughImageValuesAndDerivatives(imdata, pts, offset, theta, scaling, islog, dwelltime);
            end
            if (~useintensity)
                grad = grad(1:3,:);
            end
            %from 4xN [dFdx;dFdy;dFdtheta;dFscaling] to 3Nx1 [dFdx';dFdy';dFdtheta';dFdscaling']
            grad = reshape(grad', [], 1);
            grad = transpose(weights)*grad;
            
            if (nargout < 3)
                return;
            end
            
            dxx = h.d2Fdx2;
            dyy = h.d2Fdy2;
            dtt = h.d2Fdtheta2;
            dxy = h.d2Fdxdy;
            dxt = h.d2Fdxdtheta;
            dyt = h.d2Fdydtheta;
            dss = h.d2Fdscaling2;
            dxs = h.d2Fdxdscaling;
            dys = h.d2Fdydscaling;
            dts = h.d2Fdthetadscaling;
            
            %sparse matrix operations are sloooowwww.....
            %hx = sparse(i,j,hessval,4*N,4*N) takes 47 seconds,  
            %hess = transpose(w)*hx*w takes 20 seconds 
            %out of total 100 second run time
            if (useintensity)
                i =   [ind_x, ind_y, ind_x, ind_y, ind_t, ind_x, ind_t, ind_y, ind_t, ind_s, ind_x, ind_s, ind_y, ind_s, ind_t, ind_s];
                j =   [ind_x, ind_y, ind_y, ind_x, ind_t, ind_t, ind_x, ind_t, ind_y, ind_s, ind_s, ind_x, ind_s, ind_y, ind_s, ind_t];
                hessval = [dxx,   dyy,   dxy,   dxy,   dtt,   dxt,   dxt,   dyt,   dyt,   dss,   dxs,   dxs,   dys,   dys,   dts,   dts];
  
                hx = sparse(i,j,hessval,4*N,4*N);
            else
                i =   [ind_x, ind_y, ind_x, ind_y, ind_t, ind_x, ind_t, ind_y, ind_t];
                j =   [ind_x, ind_y, ind_y, ind_x, ind_t, ind_t, ind_x, ind_t, ind_y];
                hessval = [dxx,   dyy,   dxy,   dxy,   dtt,   dxt,   dxt,   dyt,   dyt];            
                hx = sparse(i,j,hessval,3*N,3*N);
            end
            w = sparse(weights); %weights may be sparse already
            
            hess = transpose(w)*hx*w;
            
            
        end
        
        function [val, grads, hess] = rasterLineThroughImageValuesAndDerivatives(imdata, pts, offset, theta, scaling, islog,dwelltime)
            %function [val, grads, hess] = rasterLineThroughImageValuesAndDerivatives(pts, offset, theta, scaling, islog)
            %calculates values and derivatives with respect to offset and theta given
            %scan along pts through image defined by imdata
            %
            %pts 2xN vector of scan locations
            %offset 2xN vector of offsets
            %theta 1xN vector of rotations
            %scaling 1xN vector of intensity scaling
            %
            %islog - if true, instead of multiplying F by scaling at each point,
            %        we add log(scaling) at each point
            %
            %dwelltime - (1xN) optional: multiplier for rate function and
            %           derivatives at each point, accounts for conversion
            %           from rate to expected counts in dwell calculation
            %           NOT a fit parameter, no derivatives are taken w.r.t
            %           dwelltime -- dwelltime should only be used if islog
            %           is false
            % transformation from scan coordingates to stabilized coordinates
            % u = cos(theta).*(pts(1,:) + offset(1,:)) -  sin(theta).*(pts(2,:) + offset(2,:))
            % v = sin(theta).*(pts(1,:) + offset(1,:)) +  cos(theta).*(pts(2,:) + offset(2,:))
            %
            %assumes imdata has the following gridded interpolants
            %F - value
            %dFdu, dFdv - spatial derivatives
            %d2Fdu2, d2Fdv2, d2Fduduv - spatial second derivatives
            %
            %F = F(u,v): functional form, not matrix form.
            
            u = cos(theta).*(pts(1,:) + offset(1,:)) -  sin(theta).*(pts(2,:) + offset(2,:));
            v = sin(theta).*(pts(1,:) + offset(1,:)) +  cos(theta).*(pts(2,:) + offset(2,:));
            existsAndDefault('dwelltime', 1);
            if(islog)
                dwelltime = 1;
            end
            %F = dwelltime.*interp2(imdata.u, imdata.v, imdata.F', u,v, 'linear', imdata.Fextrap);
            F = dwelltime.*imdata.F(u,v);
            if (islog)
                val = sum(F + log(scaling));
            else
                val = sum(scaling.*F);
            end
            if (nargout < 2)
                return;
            end
            dFdu = dwelltime.*imdata.dFdu(u,v);
            dFdv = dwelltime.*imdata.dFdv(u,v);
            
            
            dudx = cos(theta);
            dvdx = sin(theta);
            dudy = -sin(theta);
            dvdy = cos(theta);
            dudtheta = -sin(theta).*(pts(1,:) + offset(1,:)) -  cos(theta).*(pts(2,:) + offset(2,:));
            dvdtheta = cos(theta).*(pts(1,:) + offset(1,:)) -  sin(theta).*(pts(2,:) + offset(2,:));
            
            
            dFdx = dFdu.*dudx + dFdv.*dvdx;
            dFdy = dFdu.*dudy + dFdv.*dvdy;
            dFdtheta = dFdu.*dudtheta + dFdv.*dvdtheta;
            if (islog)
                grads = [dFdx;dFdy;dFdtheta;1./scaling];
            else
                grads = [dFdx.*scaling;dFdy.*scaling;dFdtheta.*scaling;F];
            end
            if (nargout < 3)
                return;
            end
            
            
            d2Fdu2 = dwelltime.*imdata.d2Fdu2(u,v);
            d2Fdv2 =  dwelltime.*imdata.d2Fdv2(u,v);
            d2Fdudv =  dwelltime.*imdata.d2Fdudv(u,v);
            d2udtheta2 = -u;
            d2vdtheta2 = -v;
            d2udxdtheta = -sin(theta);
            d2vdxdtheta = cos(theta);
            d2udydtheta = -cos(theta);
            d2vdydtheta = -sin(theta);
            
            
            hess.d2Fdx2 = d2Fdu2.*(dudx.^2) + 2*d2Fdudv.*dudx.*dvdx + d2Fdv2.*(dvdx.^2);
            hess.d2Fdy2 = d2Fdu2.*(dudy.^2) + 2*d2Fdudv.*dudy.*dvdy + d2Fdv2.*(dvdy.^2);
            hess.d2Fdxdy = d2Fdu2.*(dudy.*dudx) + d2Fdudv.*(dudy.*dvdx + dudx.*dvdy) + d2Fdv2.*(dvdy.*dvdx);
            hess.d2Fdtheta2 = d2Fdu2.*(dudtheta.^2) + 2*d2Fdudv.*dudtheta.*dvdtheta + d2Fdv2.*(dvdtheta.^2) + ...
                dFdu.*d2udtheta2 + dFdv.*d2vdtheta2;
            hess.d2Fdxdtheta = d2Fdu2.*(dudtheta.*dudx) + d2Fdudv.*(dudtheta.*dvdx + dudx.*dvdtheta) + d2Fdv2.*(dvdtheta.*dvdx) + ...
                dFdu.*d2udxdtheta + dFdv.*d2vdxdtheta;
            hess.d2Fdydtheta = d2Fdu2.*(dudtheta.*dudy) + d2Fdudv.*(dudtheta.*dvdy + dudy.*dvdtheta) + d2Fdv2.*(dvdtheta.*dvdy) + ...
                dFdu.*d2udydtheta + dFdv.*d2vdydtheta;
            if (islog)
                hess.d2Fdxdscaling = zeros(size(hess.d2Fdx2));
                hess.d2Fdydscaling = zeros(size(hess.d2Fdx2));
                hess.d2Fdthetadscaling = zeros(size(hess.d2Fdx2));
                hess.d2Fdscaling2 = -1./scaling.^2;
            else
                hess.d2Fdx2 =  hess.d2Fdx2.*scaling;
                hess.d2Fdy2 = hess.d2Fdy2.*scaling;
                hess.d2Fdxdy =  hess.d2Fdxdy.*scaling;
                hess.d2Fdtheta2 = hess.d2Fdtheta2.*scaling;
                hess.d2Fdxdtheta  = hess.d2Fdxdtheta.*scaling;
                hess.d2Fdydtheta =  hess.d2Fdydtheta.*scaling;
                hess.d2Fdscaling2 = zeros(size(dFdx));
                hess.d2Fdxdscaling = dFdx;
                hess.d2Fdydscaling = dFdy;
                hess.d2Fdthetadscaling = dFdtheta;
                
            end
            
        end
        
        function immodel = createImModel(immodel)
            %immodel should have fields
            %u
            %v
            %F(u,v) -- functional, not matrix form
            %Fextrap
            %
            %all u should be equally spaced, all v should be equally spaced
            u = immodel.u;
            v = immodel.v;
            if ~isfield(immodel, 'Fextrap')
                immodel.Fextrap = min(immodel.F,[1 2]);
            end

            newu = linspace(min(u), max(u), ceil(max(u)-min(u))*10); %100 nm resampled pixels
            newv = linspace(min(v), max(v), ceil(max(v)-min(v))*10);
            expandFirst = false;
            if (expandFirst)
                immodel.u = newu;
                immodel.v = newv;
                [immodel.uu, immodel.vv] = ndgrid(newu, newv);
                
                immodel.F_im = interpn(u, v, immodel.F, immodel.uu, immodel.vv, 'cubic');
                
                deltau = median(diff(immodel.u),'omitnan');
                deltav = median(diff(immodel.v),'omitnan');
                
                %Fadid and Simoncelli coefficients:
                %from https://en.wikipedia.org/wiki/Image_derivatives#Farid_and_Simoncelli_Derivatives
                %see also paper at doi 10.1109/TIP.2004.823819
                k  = [0.030320  0.249724  0.439911  0.249724  0.030320];
                d  = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
                d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
                immodel.dFdu_im = conv2(d, k, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./deltau;
                immodel.dFdv_im = conv2(k, d, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./deltav;
                immodel.d2Fdu2_im = conv2(d2,k, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./(deltau.^2);
                immodel.d2Fdv2_im = conv2(k, d2, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./(deltav.^2);
                immodel.d2Fdudv_im = conv2(k, d, padarray(immodel.dFdu_im, [2 2], 'replicate'), 'valid')./deltav;

                [uu,vv] = ndgrid([2*newu(1)-newu(2) newu 2*newu(end)-newu(end-1)],[2*newv(1)-newv(2) newv 2*newv(end)-newv(end-1)]);
        
                fields = {'F', 'dFdu', 'dFdv', 'd2Fdu2', 'd2Fdv2', 'd2Fdudv'};
                extrapval = [immodel.Fextrap, 0, 0, 0, 0, 0];
                for j = 1:length(fields)
                    immodel.(fields{j}) = griddedInterpolant(uu, vv, padarray(immodel.([fields{j} '_im']), [1 1], extrapval(j)), 'linear', 'nearest');
                end
            else
               
                
                deltau = median(diff(immodel.u),'omitnan');
                deltav = median(diff(immodel.v),'omitnan');
                
                %Fadid and Simoncelli coefficients:
                %from https://en.wikipedia.org/wiki/Image_derivatives#Farid_and_Simoncelli_Derivatives
                %see also paper at doi 10.1109/TIP.2004.823819
                k  = [0.030320  0.249724  0.439911  0.249724  0.030320];
                d  = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
                d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
                immodel.F_im = immodel.F;
                immodel.dFdu_im = conv2(d, k, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./deltau;
                immodel.dFdv_im = conv2(k, d, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./deltav;
                immodel.d2Fdu2_im = conv2(d2,k, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./(deltau.^2);
                immodel.d2Fdv2_im = conv2(k, d2, padarray(immodel.F_im, [2 2], 'replicate'), 'valid')./(deltav.^2);
                immodel.d2Fdudv_im = conv2(k, d, padarray(immodel.dFdu_im, [2 2], 'replicate'), 'valid')./deltav;

                [uu,vv] = ndgrid([2*newu(1)-newu(2) newu 2*newu(end)-newu(end-1)],[2*newv(1)-newv(2) newv 2*newv(end)-newv(end-1)]);
        
                fields = {'F', 'dFdu', 'dFdv', 'd2Fdu2', 'd2Fdv2', 'd2Fdudv'};
                extrapval = [immodel.Fextrap, 0, 0, 0, 0, 0];
                immodel.u = newu;
                immodel.v = newv;
                [immodel.uu, immodel.vv] = ndgrid(newu, newv);
                
                
                for j = 1:length(fields)
                    immodel.([fields{j} '_im']) = interpn(u, v, immodel.([fields{j} '_im']), immodel.uu, immodel.vv, 'cubic');
                    immodel.(fields{j}) = griddedInterpolant(uu, vv, padarray(immodel.([fields{j} '_im']), [1 1], extrapval(j)), 'linear', 'nearest');
                end

            end

            
        end
        
    end
    
    
end

