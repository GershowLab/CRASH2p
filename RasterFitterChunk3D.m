classdef RasterFitterChunk3D < RigidAligner3D
    % Calculates image values, derivatives, and Hessians 
    % at specified locations where photon was observed
    
    
    
    properties
        k_theta = 0.01; %angle continuity penalty  = -k(cos(theta_j - theta_j-1))
        k_delta = 0.002; %continuity weight for delta vector -- translation = .5 k (x_j - x_j-1)^2
        k_weight_method = 0; %0 = 1 / 1 = mult by # photons / 2 = mult by log (# photons)
        fixedBefore = [];
        fixedAfter = [];
        fit2D = false;
    end
    
    
    
    methods
       
        
        function [avec, pvec] = interpolateFit (rfc, ti)
            % function avec = interpolateFit (rfc, ti)
            % interpolates the fit represented by rfc to a new time window,
            % respecting nonlinear constraints
            % rfc < RasterFitterChunk vector
            % ti - vector of time points
            % avec - length(ti) x navec(6) interpolated affine transform
            
            p = cat(1,rfc.fitParameters);
            p(:,1:3) = unwrap(p(:,1:3),[],2);
          
            avec = zeros(length(ti), 12);
            
            tx = 0.5 * ([rfc.startTime] + [rfc.endTime]);
            pvec = interp1(tx, p, ti, 'linear');
            pvec(:,1:3) = mod(p(:,1:3), 2*pi);
            for j = 1:size(pvec,1)
                avec(j,:) = reshape(RigidAligner3D.rigidTransform(pvec(j,:)), 1, []);
            end
            
        end

        function [rfc, val] = initialDivisionAndAlignment(rfc, frameTime)
            rfc = rfc.splitIntoFrames(frameTime);
            for j = 1:length(rfc)
                [~,rfc(j)] = rfc(j).fitRigidFromScratch();
            end
            val = [rfc.funPerPhoton];
        end
        
        
        function [rfc, pvec] = fitChunks (rfc, pvec_init, binsize, chunksize, overlap)
            if (nargin < 2)
                pvec_init = [];
            end
            if (nargin < 3)
                binsize = [];
            end
            nrfc = length(rfc);
            nchunks = ceil((nrfc-chunksize)/(chunksize-overlap));
            if (nchunks <= 1)
                [rfc,pvec] = rfc.fit(pvec_init, binsize);
                return;
            end
            si = floor(linspace(1,nrfc,nchunks+1));
            si = si(1:end-1);
            for j = 1:length(si)
                
                inds = si(j):min(si(j)+chunksize-1, length(rfc));
                chunk = rfc(inds);
                if (inds(1) > 1)
                    chunk(1).fixedBefore = rfc(inds(1)-1).fitParameters;
                end
                if (~isempty(pvec_init) && size(pvec_init,1) >= length(rfc))
                    pvi = pvec_init(inds, :);
                else
                    pvi = pvec_init;
                end
                chunk = chunk.fit(pvi, binsize);
                for k = 1:length(chunk)
                    rfc(inds(k)).fitParameters = chunk(k).fitParameters;
                    rfc(inds(k)).fitAffine = chunk(k).fitAffine;
                    rfc(inds(k)).initAffine = chunk(k).initAffine;
                end
            end
            pvec = cat(1,rfc.fitParameters);
            
        end
        
        function rfc2 = combine(rfc,n)
            %function rfc2 = combine(rfc,n)
            %
            %combines a set of rfc into n rfcs, 
            %e.g. if rfc has length 10 and n = 2, then
            %rfc2(1) combines rfc(1:5) and rfc2(2) combines rfc(6:10)
            if (nargin < 2 || n == 1)
                rfc2 = rfc(1).setPhotonLoc([rfc.ploc], [rfc.pphase], [rfc.ptime]);% RasterFitterChunk3D(rfc(1).immodel, , rfc(1).z_scale);
                return
            end
            if (n >= length(rfc))
                rfc2 = rfc; %don't bother with new allocation/copying
                return;
            end
            ii = ceil(linspace(1,length(rfc),n+1));
            for j = 1:n
                if (j == n)
                    inds = ii(j):ii(end);
                else
                    inds = ii(j):(ii(j+1)-1);
                end
                rfc2(j) = combine(rfc(inds),1); %#ok<AGROW>
            end
        end
            
        
        function [theta_z, fval, extras] = fitInitialThetaZ(rfc, binsize, display)
          
           if (nargin < 2)
                binsize = [];
           end
           if (nargin < 3)
               display = false;
           end
           
            
           n = min(length(rfc), 2.^(0:ceil(log2(length(rfc)))));
           tzprev = [];
           ixprev = [];
           for j = 1:length(n)
               rfc2 = combine(rfc,n(j));
               if (j == 1)
                   %find initial best orientation for whole chunk 
                   theta_z = mod(rfc2.orientRotationBounded([-0.1 2*pi],60), 2*pi);
               else
                   if (n(j - 1) <= 1)
                       theta_z = repmat(theta_z, [1 n(j)]);
                   else
                       theta_z = mod(interp1(linspace(1,length(rfc),n(j-1)), unwrap(theta_z),linspace(1,length(rfc),n(j)), 'linear'),2*pi);
                   end
                   for k = 1:length(rfc2)
                       rfc2(k).lb(1) = -pi;
                       rfc2(k).ub(1) = 3*pi;
                   end
                   
                   [theta_z,fval,extras] = rfc2.fitThetaZOnly(theta_z, binsize);
                   theta_z = mod(theta_z, 2*pi);
                   if (display)
                       figure(1); clf(1);
                       plot(ixprev, tzprev, 'b.-', linspace(1,length(rfc),n(j)), unwrap(theta_z)*180/pi, 'r.-'); drawnow; 
                   end
                   ixprev = linspace(1,length(rfc),n(j));
                   tzprev = unwrap(theta_z)*180/pi;
               end
           end
           
           
           
        end
        
        function [theta, fval, extras] = fitThetaZOnly(rfc, theta_z, binsize)
            if (nargin < 3)
                binsize = [];
            end
            [rfc.fit2D] = deal(true);
            lb = cat(1,rfc.lb);
            ub = cat(1,rfc.ub);
            lb = lb(:,1);
            ub = ub(:,1);
             
            for j = 1:length(rfc)
                rfc(j) = rfc(j).binPhotonsByBinsize(binsize);
            end
            [v,g] = rfc.evalAndContinuityJustThetaZ(theta_z);
            op = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Display', 'off', 'MaxFunctionEvaluations', 3e4);
            [theta, fval,extras.exitflag,extras.output,extras.lambda,extras.grad,extras.hessian] = fmincon(@rfc.evalAndContinuityJustThetaZ, theta_z, [], [], [], [], lb, ub, [], op);        
        end
        
        function [rfc, pvec, fval, extras, pvec_init, fval_image_only] = fit (rfc, pvec_init, binsize)
            % function [rfc, pvec, fval, extras, avec_init] = fit (rfc, theta0)
            % function [rfc, pvec, fval, extras, avec_init] = fit (rfc, pvec_init)
            
            %inputs
            % rfc < RasterFitterChunk - vector
            % pvec_init -  length(rfc) x nparams(6 or 3): initial guess for
            %       affine transformation for each rfc (or 1 x nparams)
           
            %avec - length(rfc) x 12 - vectorized version of affine
            %   matrix for each rfc
            %rfc - vector of raster fitter chunks
            %fval - output of optimizer
            %extras - output of optimizer
            %avec_init - initial input to optimizer
            
           
            
            if (nargin < 3 || isempty(binsize))                     
                binsize = [median(diff(rfc(1).immodel.u)) median(diff(rfc(1).immodel.v)) median(diff(rfc(1).immodel.w))];            
            end
            
            if (nargin < 2 || isempty(pvec_init))
                pvec_init = [0 0 0 0 0 0];
            end
            
            if (size(pvec_init,1) == 1)
                pvec_init = repmat(pvec_init, [length(rfc), 1]);
            end
            
            [rfc.fit2D] = deal(size(pvec_init,2) == 3);
            
            
            for j = 1:length(rfc)
                rfc(j) = rfc(j).binPhotonsByBinsize(binsize);
            end
                        
            lb = cat(1,rfc.lb);
            ub = cat(1,rfc.ub);
            
            if (rfc(1).fit2D)
                lb = lb(:,[1 4 5]);
                ub = ub(:,[1 4 5]);
            end
            
%             
%             
%             lb = repmat([-pi/4 -pi/4 -pi/4 -100 -100 -100], [length(rfc), 1]); 
%             ub = repmat([2*pi+pi/4 pi/4 pi/4  100 100 100], [length(rfc), 1]);
%             
            
            %convert 2D avec to 1D parameter vector for fitting
            pvec_init_1 = RasterFitterChunk3D.pvec2to1(pvec_init);
            lb = RasterFitterChunk3D.pvec2to1(lb);
            ub = RasterFitterChunk3D.pvec2to1(ub);
            
          
            
            %run nonlinear minimizer 
            op = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Display', 'notify', 'MaxFunctionEvaluations', 3e4);
            [pvec_1, fval,extras.exitflag,extras.output,extras.lambda,extras.grad,extras.hessian] = fmincon(@rfc.evalAndContinuity, pvec_init_1, [], [], [], [], lb, ub, [], op);
            
            %convert parameter vectors back to 2D avec
            pvec = RasterFitterChunk3D.pvec1to2(pvec_1, rfc(1).fit2D);
            pvec_init = RasterFitterChunk3D.pvec1to2(pvec_init_1, rfc(1).fit2D);
            if (rfc(1).fit2D)
                z = pvec(:,1)*0;
                pvec = [pvec(:,1) z z pvec(:,2:3) z];
                pvec_init = [pvec_init(:,1) z z pvec(:,2:3), z];
            end
                
            
            pvec(:,1) = mod(pvec(:,1), 2*pi);
            %load results back into rfc
            for j = 1:length(rfc)
                %aa = avec(j, :);
                
                rfc(j).fitParameters = pvec(j,:);
                rfc(j).fitAffine = RigidAligner3D.rigidTransform(pvec(j,:));
                rfc(j).initAffine =RigidAligner3D.rigidTransform(pvec_init(j,:));
                fval_image_only(j) = eval(rfc(j),pvec(j,:)); %#ok<AGROW>
            end
            
        end 
        
        
        
        function [g,gd,ge, ged, gc, gcd] = gradientTest(rfc, pvec, delta, nparams)
            if (nargin < 2 || isempty(pvec))
                pvec = cat(1,rfc.fitParameters);
            end
            if (nargin < 3 || isempty(delta))
                delta = 1e-3;
            end
            
            pvec = RasterFitterChunk3D.pvec2to1(pvec); %no effect if already 1D
            
            if (nargin < 4)
                nparams = length(pvec);
            end
            [v,g] = rfc.evalAndContinuity(pvec);
            
            gd = 0*g;
            for j = 1:nparams
                d = 0*pvec;
                d(j) = delta;
                gd(j) = (rfc.evalAndContinuity(pvec+d)-v)/delta;
            end
            
            [v,ge] = rfc.eval(pvec);
            
            ged = 0*g;
            for j = 1:nparams
                d = 0*pvec;
                d(j) = delta;
                ged(j) = (rfc.eval(pvec+d)-v)/delta;
            end
            
            [v,gc] = rfc.continuity(pvec);
            
            gcd = 0*g;
            for j = 1:nparams
                d = 0*pvec;
                d(j) = delta;
                gcd(j) = (rfc.continuity(pvec+d)-v)/delta;
            end
            
        end
        
        
       
        
        function [val, grad] = continuity (rfc, pvec_1)
            % [val, grad, hess] = continuity (rfc, avec_1)
            % part of objective function for fitter, 
            
            if (length(rfc) == 1)
                val = 0;
                grad = 0*pvec_1;
                return;
            end
            pvec = RasterFitterChunk3D.pvec1to2(pvec_1, rfc(1).fit2D);
            
            for j = 1:length(rfc)
                if (rfc(j).k_weight_method == 0)
                    continue;
                end
                pw = sum(rfc(j).pweightvec);
                if (rfc(j).k_weight_method == 2)
                    pw = log(pw);
                end
                rfc(j).k_delta = rfc(j).k_delta*pw;
                rfc(j).k_theta = rfc(j).k_theta*pw;
            end
            pvarr = pvec;

            if (size(pvec, 2) == 3) %2D
                %weight is nparams x n
                weight = [cat(2,rfc(1:end-1).k_theta);repmat(cat(2,rfc(1:end-1).k_delta), [2 1])];
                gi = 1:size(pvec,1);
                if (~isempty(rfc(1).fixedBefore))
                    pvarr = [rfc(1).fixedBefore([1 4 5]); pvarr];
                    weight = weight(:,[1 1:end]);
                    gi = gi+1;
                end
                if (~isempty(rfc(end).fixedAfter))
                    pvarr = [pvarr;rfc(end).fixedAfter([1 4 5])];
                    weight = [weight rfc(end).k_theta*[1 0 0 ]'+ rfc(end).k_delta*[0 1 1]']; 
                end
                
                ai = 1; %angle inds
                pi = 2:3;%position inds
                n = 3;

            else
                
                
                %weight is nparams x n
                weight = [repmat(cat(2,rfc(1:end-1).k_theta), [3 1]);repmat(cat(2,rfc(1:end-1).k_delta), [3 1])];
                gi = 1:size(pvec,1);
                if (~isempty(rfc(1).fixedBefore))
                    pvarr = [rfc(1).fixedBefore; pvarr];
                    weight = weight(:,[1 1:end]);
                    gi = gi+1;
                end
                if (~isempty(rfc(end).fixedAfter))
                    pvarr = [pvarr;rfc(end).fixedAfter];
                    weight = [weight rfc(end).k_theta*[1 1 1 0 0 0 ]'+ rfc(end).k_delta*[0 0 0 1 1 1]']; 
                end
                ai = 1:3;
                pi = 4:6;
                n = 6;

% 
%                 delta = diff(pvarr,[],1);
% 
%                 wt = weight.';
% 
%                 val = delta; %n+/-1 x nparams
%                 val(:,1:3) = -wt(:,1:3).*cos(delta(:,1:3)); 
%                 val(:,4:end) = 0.5*wt(:,4:end).*delta(:,4:end).^2;
%                 val = sum(val, 'all');
% 
%                 wd = delta;
%                 wd(:,1:3) = wt(:,1:3).*sin(delta(:,1:3)); 
%                 wd(:,4:end) = wt(:,4:end).*delta(:,4:end);
%                 grad = [zeros([1 RasterFitterChunk.navec]); wd] - [wd; zeros([1 RasterFitterChunk.navec])];
%                 grad = RasterFitterChunk3D.pvec2to1(grad(gi,:));
            end
            
                delta = diff(pvarr,[],1);

                wt = weight.';

                
                val = delta; %n+/-1 x nparams
                val(:,ai) = -wt(:,ai).*cos(delta(:,ai)); 
                val(:,pi) = 0.5*wt(:,pi).*delta(:,pi).^2;
                val = sum(val, 'all');

                wd = delta;
                wd(:,ai) = wt(:,ai).*sin(delta(:,ai)); 
                wd(:,pi) = wt(:,pi).*delta(:,pi);
                grad = [zeros([1 n]); wd] - [wd; zeros([1 n])];
                grad = RasterFitterChunk3D.pvec2to1(grad(gi,:));
        end
        
        function [vlist, g, g2] = gradTest(rfc, pvec, delta)
           
            
            if (rfc.invertEval)
                m = -rfc.eval_weight;
            else
                m = rfc.eval_weight;
            end
            
            [tmat, gradmat] = RasterFitterChunk3D.rotAndStetchInZ(pvec);
            
            
            newpts = tmat*rfc.plocvec; %3x4 * 4xN --> 3xN
            x = newpts(1,:);
            y = newpts(2,:);
            z = newpts(3,:);
            
            vlist = m * rfc.pweightvec.*rfc.immodel.F(x,y,z);
            
            tmat2 = RasterFitterChunk3D.rotAndStetchInZ((pvec + delta)./ts);
            newpts2 = tmat2*rfc.plocvec; %3x4 * 4xN --> 3xN
            x2 = newpts2(1,:);
            y2 = newpts2(2,:);
            z2 = newpts2(3,:);
            
            vlist2 = m * rfc.pweightvec.*rfc.immodel.F(x2,y2,z2);
            
            g2 = vlist2 - vlist;
            
            gu = m*(rfc.pweightvec.*rfc.immodel.dFdu(x,y,z));
            gv = m*(rfc.pweightvec.*rfc.immodel.dFdv(x,y,z));
            gw = m*(rfc.pweightvec.*rfc.immodel.dFdw(x,y,z));
            
            imgrad = [gu(:)';gv(:)';gw(:)'];
            
            g = 0*g2;
            for j = 1:RasterFitterChunk3D.nparams
                g = g + delta(j)*sum(imgrad .* (gradmat(:,:,j)*rfc.plocvec),1);
            end
            
            
        end
        
        function [val, grad] = eval(rfc, pvec_1)
            pvec = RasterFitterChunk3D.pvec1to2(pvec_1, rfc(1).fit2D);
            val = 0;
            grad = 0*pvec;
            for j = 1:length(rfc)
                if (nargout > 1)
                    [v,g] = eval@RigidAligner3D(rfc(j),pvec(j,:));
                    val = val + v;
                    grad(j,:) = g;
                else
                    val = val + eval@RigidAligner3D(rfc(j),pvec(j,:));
                end
            end
            grad = RasterFitterChunk3D.pvec2to1(grad);
        end
        
        function [val,grad] = evalAndContinuityJustThetaZ(rfc, theta_z)
            pvec_1 = zeros(1, 3*length(theta_z));
            pvec_1(1:3:end) = theta_z;
            if (nargout <= 1)
                val = rfc.evalAndContinuity(pvec_1);
                return;
            end
            [val,g] = rfc.evalAndContinuity(pvec_1);
            grad = g(1:3:end);
        end
        
        function [val, grad] = evalAndContinuity(rfc, pvec_1)
            if (nargout <= 1)
                val = rfc.eval(pvec_1) + rfc.continuity(pvec_1);
                return;
            end
            if (nargout == 2)
                [v1,g1] = rfc.eval(pvec_1);
                [v2,g2] = rfc.continuity(pvec_1);
                val = v1+v2;
                grad = g1 + g2;
                return;
            end
        end
        
        function rfc = splitAtTime(rfc, splitTime)
            %rfc = splitAtTime(rfc, splitTime)
            %rfc out is length(splitTime)+1;  split at each splitTime
            if (size(rfc.ptime,2) ~= size(rfc.ploc,2))
                warning ('cannot split by time because time is missing or sized incorrectly');
                size(rfc.ptime)
                size(rfc.plocvec)          
                return;
            end
            if (any(splitTime < rfc.ptime(1) | splitTime > rfc.ptime(end)))
                warning ('split time is not in range of photon times');
                return;
            end
            for j = 1:length(splitTime)
                ind(j) = find(rfc.ptime < splitTime(j), 1, 'last'); %#ok<AGROW>
            end
            rfc = rfc.splitAtInd(ind);
        end
        
        function rfc = setBoundsNearFit(rfc, varargin)
            % rfc = setBoundsNearFit(rfc, fitparameters, delta)
           
            if (length(rfc) > 1)
                for j = 1:length(rfc)
                    rfc(j) = rfc(j).setBoundsNearFit(varargin{:});
                end
                return;
            end
 
            
            if (length(varargin) < 1 || isempty(varargin{1}))
                fitparams = rfc.fitParameters;
            else
                fitparams = varargin{1};
            end
             if (length(varargin) < 2 || isempty(varargin{2}))
                delta = [deg2rad(10) deg2rad(5) deg2rad(5) 2 2 2];
            else
                delta = varargin{2};
             end
             delta = abs(delta);
             rfc.ub = fitparams + delta;
             rfc.lb = fitparams - delta;
             
        end     
        function rfc = splitAtInd(rfc, ind)
            %rfc = splitAtInd(rfc, ind)
            %splits at number of photons defined by ind
            if (any(ind < 1 | ind >= size(rfc.ploc,2)))
                warning ('split is out of range')
                return;
            end
            ii = [0 ind length(rfc.ptime)];
            
            ploc = rfc.ploc;
            ptime = rfc.ptime;
            pphase = rfc.pphase;
            
            dloc = rfc.dloc;
            dtime = rfc.dtime;
            dphase = rfc.dphase;
            dweight = rfc.dweight;

            %diffusion argument
            %p(a1,a2) = exp(-(a1-a2)^2/4Dt^2)
            %log p = 1/2 k (a1-a2)^2; k = 1/2Dt^2
            %t' = t/N --> k' = N^2k            
            %divide spring into equal lengths
            %k d^2 = N k' (d/N)^2
            %k' = N k
            n = length(ii) - 1;
            rfc.k_theta = n*rfc.k_theta;
            rfc.k_delta = n*rfc.k_delta;
            
            for j = 1:n
                inds = (ii(j)+1):(ii(j+1));
                rfc(j) = rfc(1).setPhotonLoc(ploc(:,inds), pphase(inds), ptime(inds));
                ii_dwell = dtime >= ptime(inds(1)) & dtime <= ptime(inds(end));
                rfc(j) = rfc(j).setDwellLoc(dloc(:,ii_dwell), dphase(:,ii_dwell), dweight(:,ii_dwell),dtime(:,ii_dwell));
            end

        end
        
        function rfc = merge(rfcarr)
            photonloc = [rfcarr.ploc];
            phottime = [rfcarr.ptime];
            rfc = rfcarr(1);
            kfactor = (rfcarr(1).endTime - rfcarr(1).startTime) / (rfcarr(end).endTime - rfcarr(1).startTime);
            rfc.k_a = rfc.k_a * kfactor;
            rfc.k_delta = rfc.k_delta * kfactor;
            rfc = rfc.setPhotonLoc(photonloc, phottime);        
        end
        
        function rfc = splitIntoN(rfc, N)
            %rfc = splitIntoN(rfc, N)
            %splits into N parts with equal numbers of photons
           if (N <= 1)
               return;
           end
           ivec = linspace(1,size(rfc.plocvec,2), N+1);
           rfc = rfc.splitAtInd(ivec(2:(end-1)));
        end
        
        function rfc = splitIntoNEqualTimeArr(rfc,N)
            rfc2 = repmat(rfc(:)', [N 1]);
            for j = 1:length(rfc)
                rfc2(:,j) = reshape(rfc(j).splitIntoNEqualTime(N),[],1);
            end
            rfc = reshape(rfc2, 1, []);
            [~,I] = sort([rfc.startTime]);
            rfc = rfc(I);
        end
        
        function rfc = splitIntoNEqualTime(rfc, N)
            % rfc = splitIntoNEqualTime(rfc, N)
            % splits into N parts of equal duration
           if (N <= 1)
               return;
           end
           if (size(rfc.ptime,2) ~= size(rfc.ploc,2))
                warning ('cannot split by time because time is missing or sized incorrectly');
                size(rfc.ptime)
                size(rfc.plocvec)          
                return;
            end
           tvec = linspace(min(rfc.ptime),max(rfc.ptime), N+1);
           rfc = rfc.splitAtTime(tvec(2:(end-1)));
        end
        
        
        
        function rfc = splitIntoFrames(rfc, frameTime)
            rfc = rfc.splitIntoNEqualTime(round((rfc.endTime - rfc.startTime) / frameTime));
        end
        
        function [rfc,pvec_init] = initializeSegment(rfc, pvec_left, pvec_right, maxDeltaPerChunk, pvec_init)
            tbounds = unwrap(mod([pvec_left(1) pvec_right(1)], 2*pi));
            pvec_left(1,:) = tbounds(1);
            pvec_right(1,:) = tbounds(2);
            rfc(1).fixedBefore = pvec_left;
            rfc(end).fixedAfter = pvec_right;
            rfc(1).lb = pvec_left - maxDeltaPerChunk;
            rfc(1).ub = pvec_left + maxDeltaPerChunk;
            for j = 2:length(rfc)
                rfc(j).lb = rfc(j-1).lb - maxDeltaPerChunk;
                rfc(j).ub = rfc(j-1).ub + maxDeltaPerChunk;
            end
            rfc(end).lb = pvec_right - maxDeltaPerChunk;
            rfc(end).ub = pvec_right + maxDeltaPerChunk;
            for j = (length(rfc)-1):-1:1
                rfc(j).lb = max(rfc(j).lb, rfc(j+1).lb - maxDeltaPerChunk);
                rfc(j).ub = min(rfc(j).ub, rfc(j+1).ub + maxDeltaPerChunk);
            end
                    
            
            anglemin = [min(tbounds)-2*pi, -pi/4, -pi/4];
            anglemax = [max(tbounds)+2*pi, pi/4, pi/4];
            
            deltamin = min(min(pvec_left(4:6), pvec_right(4:6))-5, [-10 -10 -10]);
            deltamax = max(max(pvec_left(4:6), pvec_right(4:6))+5, [10 10 10]);
            
            lbglobal = [anglemin, deltamin];
            ubglobal = [anglemax, deltamax];
            
            for j = 1:length(rfc)
                rfc(j).lb = max(rfc(j).lb, lbglobal);
                rfc(j).ub = min(rfc(j).ub, ubglobal);
            end
            if (nargin < 5 || isempty(pvec_init))
               pvec_init = interp1([0 length(rfc)+1], [pvec_left; pvec_right], 1:length(rfc));
            end
            for k = 1:3
                %force everything within bounds
                for j = 1:length(rfc)
                    while(pvec_init(j,1) < rfc(j).lb(1))
                        pvec_init(j,1) = pvec_init(j,1) + 2*pi;
                    end
                    while(pvec_init(j,1) > rfc(j).ub(1))
                        pvec_init(j,1) = pvec_init(j,1) - 2*pi;
                    end
                    if (pvec_init(j,1) < rfc(j).lb(1))
                        pvec_init(j,1) = 0.5*(rfc(j).ub(1) + rfc(j).lb(1));
                    end
                end
                %unwrap first two iterations; last iteration make sure
                %everything remains in bound
                if (k < 3)
                    pvec_init(:,1) = unwrap(pvec_init(:,1));
                end
            end
            
        end
        
        
    end
    
    methods(Static)
        
        
        
        function pvec2 = pvec1to2 (pvec1, is2D)
            
            %reshapes 1D parameter list to be nrfcXnparams sized
            if (nargin >= 2 && is2D)
                Np = 3;
            else
                Np = 6;%RasterFitterChunk3D.nparams;
            end
            if (numel(pvec1) <= Np) %allow for possibility of simple 2D rotation
                pvec2 = pvec1;
                return;
            end
            pvec2 = reshape(pvec1, Np, [])';          
        end
        function pvec1 = pvec2to1 (pvec2)
            pvec1 = reshape(pvec2', 1, []);
        end
        function [theta_z, theta_y, theta_x, dx, dy, dz] = pvecToParams(pvec2)
            theta_z = pvec2(:,1);
            theta_y = pvec2(:,2);
            theta_x = pvec2(:,3);
            dx = pvec2(:,5);
            dy = pvec2(:,6);
            dz = pvec2(:,7);
        end
        function pvec2 = paramsToPvec(theta_z, theta_y, theta_x, dx, dy, dz)
            pvec2 = [theta_z, theta_y, theta_x, dx, dy, dz];
        end
        
        
    end
    
end

