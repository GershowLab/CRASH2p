classdef NonRigidAlignerBSpline < RigidAligner3D
    % Calculates image values, derivatives, and Hessians 
    % at specified locations where photon was observed
    
    properties (Constant)
       

        

    end
    
    properties
        
        fixedAffine;
        npts = [6,6,4];
        geometry_weight = 0; %set > 0 to regularize grid using elastic springs 
        anchorlength = 10000; %don't constrain points to remain near initial position beyond 0.4 rule 
        %fitparameters = 4D - 3 x npts(1) x npts(2) x npts(3)
        dweight_correction = 1; %updated in adjustDwellToMatchRate()
    end
    
    properties (Access = public)
        pB %{nx x Np ,ny x Np,nz x Np} 
        dB %{nx x Nd ,ny x Nd,nz x Nd} 
        x0 %3 components - origin of mesh
        delta %3 components - step size of mesh

    end
    
    
    methods
       
        

        function nra = NonRigidAlignerBSpline(ra) %immodel, photonloc, photonphase, photontime,  z_scale)
            %rfc = RasterFitChunk(photonloc, immodel)
            %photonloc = [xphoton;yphoton;time] (time is optional but used
            %to split rfc into smaller components)
            %immodel has fields u, v, F - function of (u,v) -- u is
            %horizontal axis (x-axis) -- ie F is functional, not matrix form
            %immodel can have Fextrap 

            nra.z_scale = ra(1).z_scale;
            nra = nra.setPhotonLoc(ra);
            nra.immodel = ra(1).immodel;
            nra.xrange = ra(1).xrange;
            nra.yrange = ra(1).yrange;
            nra.czrange = ra(1).czrange;           
        end
        

        %converts 1D parameter vector into set of multi-level displacement
        %matrices, according to cells npts{}, delta{}, and x0{} in nra
        %if v isn't provided, initializes disp mats with all 0s
        function s = vecToDisplacmentStructs(nra, v)
            if (~iscell(nra.npts))
                s.npts = nra.npts;
                s.numelems = 3*prod(s.npts);
                s.delta = nra.delta;
                s.x0 = nra.x0;
                if (nargin >= 2)
                    s.dispmat = reshape(v, [3 s.npts]);
                else
                    s.dispmat = zeros([3 s.npts]);
                end
                s.dweight_correction = nra.dweight_correction;
                return
            end
            offset = 0;
            for j = 1:length(nra.npts)
                s(j).npts = nra.npts{j};
                s(j).numelems = 3*prod(s(j).npts);
                s(j).delta = nra.delta{j};
                s(j).x0 = nra.x0{j};
                if (nargin >= 2)
                    s(j).dispmat = reshape(v(offset + (1:s(j).numelems)), [3 s(j).npts]);
                    offset = offset + s(j).numelems;
                else
                    s(j).dispmat = zeros([3 s(j).npts]);
                end
                s(j).dweight_correction = nra.dweight_correction;
            end
        end
        
        %converts set of multi-level displacement into single parameter
        %vector
        function v = displacementStructsToVec(nra, s) %#ok<INUSL> 
            v = zeros([1 sum([s.numelems])]);
            offset = 0;
            for j = 1:length(s)
                v(offset + (1:s(j).numelems)) = s(j).dispmat(:);
                offset = offset + s(j).numelems;
            end
        end


        %finds locations of original and displaced control points at each
        %level
        function [oldPts, newPts] = getControlPointLocations(nra, j) 
            if (nargin < 2 && iscell(nra.npts))
                for j = 1:length(nra.npts)
                    [oldPts{j}, newPts{j}] = nra.getControlPointLocations(j); %#ok<AGROW> 
                end
                return;
            end
            if (nargin >= 2)
                np = nra.npts{j};
                del = nra.delta{j};
                X0 = nra.x0{j};
                if (nargout > 1)
                    dm = nra.fitParameters(j).dispmat;
                end
            else
                np = nra.npts;
                del = nra.delta;
                X0 = nra.x0;
                if (nargout > 1)
                    dm = nra.fitParameters.dispmat;
                end
            end

            [xx,yy,zz] = ndgrid((0:(np(1)-1))*del(1), (0:(np(2)-1))*del(2), (0:(np(3)-1))*del(3));
            oldPts = zeros([3 size(xx)]);
            oldPts(1,:,:,:) = xx + X0(1);
            oldPts(2,:,:,:) = yy + X0(2);
            oldPts(3,:,:,:) = zz + X0(3);
            if (nargout > 1)
                newPts = oldPts + dm;
                newPts = reshape(newPts, 3, []);
            end

            oldPts = reshape(oldPts, 3, []);

        end

        %calculates parameters required for fitting
        function nra = calculateBSPlineParameters(nra)
            ppts = nra.plocvec(1:3,:);
            dpts = nra.dlocvec(1:3,:);
            mn = min([ppts dpts],[], 2);
            mx = max([ppts dpts],[], 2);
            if (~iscell(nra.npts))
                nra.npts = {nra.npts};
                nra.delta = {};
                nra.x0 = {};
            end
            for j = 1:length(nra.npts)
                nra.delta{j} = (mx-mn)./(nra.npts{j}' - 1);
                nra.x0{j} = mn;
            end
                
                
            nra.pB = nra.getBs(ppts);
            nra.dB = nra.getBs(dpts);
            

            s = nra.vecToDisplacmentStructs(); %empty creates all zeros

            for j = 1:length(s)
                for k = 1:3
                    s(j).dispmat(k,:,:,:) = 0.4*s(j).delta(k);
                end
            end
            nra.ub = nra.displacementStructsToVec(s);
            nra.lb = -nra.ub;
        end


        %finds weights of each control point corresponding to each element
        %of pts
        function B = getBs(nra, pts)
           
            for k = 1:length(nra.npts)
                normpts = (pts - repmat(nra.x0{k},[1 size(pts,2)]))./repmat(nra.delta{k},[1 size(pts,2)]); %pts runs from 0 to N-1
            
                B{k} = {zeros(nra.npts{k}(1),size(pts,2)), zeros(nra.npts{k}(2),size(pts,2)), zeros(nra.npts{k}(3),size(pts,2))};
                for j = 1:3
                  for qq = (-1):(nra.npts{k}(j)+2)
                      i = min(max(1,qq), nra.npts{k}(j));
                      s = normpts(j,:) - qq + 1;
                      B{k}{j}(i,:) = B{k}{j}(i,:) + NonRigidAlignerBSpline.chanCubicPolynomial(s); 

                       %out of range points' weights are added to boundary points weights - 
                      % i runs from 1 to N
                      % qq runs from -1 to N+2
                      % (pts - qq) : 0 runs from -N+2 to 1 (needs to run to 2)
                      %            : N-1 runs from -3 to N (needs to run from
                      %                                        -2)
                      % add 1 to make s = pts - qq + 1
                      
                      %check: for x = 2, should have B(1) = B3(2) = 0, B(2) = B3(1) = 1/6, B(3) = B3(0) = 2/3, B(4) = B3(-1) = 1/6, B(5...) = B3(s<=-2) = 0
                      %for x = 2, s runs from +2 to N-1; when qq = 3, i = 3, s=2 - 3 + 1 = 0
                  end
                
                end
            end

        end



        
        function [newpts, delta] = getTransform(nra, pts, y, z)
            %function [newpts, delta] = getTransform(nra, x, y, z)
            %function [newpts, delta] = getTransform(nra, pts)
            %applies transformation encoded by nra.fitParameters to pts

            if (nargin > 2)
                pts = [pts(:) y(:) z(:)]';
            end
            B = nra.getBs(pts);
            delta = pts*0;
            for j = 1:length(nra.fitParameters)
                delta = delta + NonRigidAlignerBSpline.calculateDisplacement(nra.fitParameters(j).dispmat, B{j});
            end
            newpts = pts + delta;
        end
        
        
        %transcribes photon and dwell locations from ra to nra
        %points are stored as initial locations, affine transform to be
        %applied (result of ra.fit) is stored with them
        function nra = setPhotonLoc(nra, ra)
            for j = 1:length(ra)
                nra.ploc{j} = ra(j).ploc;  
                nra.pphase{j} = ra(j).pphase;
                nra.ptime{j} = ra(j).ptime;
                nra.fixedAffine{j} = ra(j).fitAffine;
            end
            nra.startTime = min([nra.ptime{:}], [],'all');
            nra.endTime = max([nra.ptime{:}], [],'all');
            if all([ra(:).hasdwell]) % only set nra dwell if all ra has dwell
                for j = 1:length(ra)
                    nra.dtime{j} = ra(j).dtime;                    
                    nra.dloc{j} = ra(j).dloc;  
                    nra.dphase{j} = ra(j).dphase;
                    nra.dweight{j} = ra(j).dweight;
                end
                nra.hasdwell = true;
            end
                       
        end
        
        %TODO
        function [im, ax] = makeImage(nra, ax, usefit)
            % function [im, ax] = makeImage(ra, ax, usefit)
            % if ax not provided, uses immodel ax
            % if usefit is false, does affine transformation but not
            % nonrigid
            %may not work - binCountsAndDwell may be more appropriate
            if (nargin < 2 || isempty(ax))
                ax.x = nra.immodel.u;
                ax.y = nra.immodel.v;
                ax.z = nra.immodel.w;
            end
            if (nargin < 3)
                usefit = true;
            end

            fa = FastAssembler(ax.x, ax.y, ax.z,[],true);
            im = binCountsAndDwell(nra, fa, usefit);
%             fa.zOfP = @(p) nra.z_scale*cos(p);
% 
%             locvec = zeros([3 0]);
%             pp = [];
%             for j = 1:length(nra.ploc)
%                 valid = nra.ploc{j}(1,:) >= min(nra.xrange) & nra.ploc{j}(1,:) <= max(nra.xrange) & nra.ploc{j}(2,:) >= min(nra.yrange) & nra.ploc{j}(2,:) <= max(nra.yrange) & cos(nra.pphase{j}) >= min(nra.czrange) & cos(nra.pphase{j}) <= max(nra.czrange);
%                 pl = nra.ploc{j}(:,valid);
%                 pp = [pp nra.pphase{j}(valid)]; %#ok<AGROW> 
%                 pl(3,:) = nra.z_scale*cos(nra.pphase{j}(valid));
%                 pl(4,:) = 1;
%                 pl = nra.fixedAffine{j}*pl;
%                 locvec = [locvec pl]; %#ok<AGROW> 
%             end
% 
%             if (usefit)
%                 locvec = nra.getTransform(locvec);
%             end
% 
%             im = fa.accumVals(fa.depthCorrection(pp), locvec(1,:), locvec(2,:), locvec(3,:));            
                
        end
        
        function [counts,dwell,photonloc,photonphase,dwellloc,dwellweight] = binCountsAndDwell (nra, fa, usefit)
            % function  [counts,dwell,photonloc,dwellloc,dwellweight] = binCountsAndDwell (nra, fa, usefit)
            %
            % bins counts and dwell with
            
            if (nargin < 2 || isempty(fa))
                fa = FastAssembler(nra.immodel.u, nra.immodel.v, nra.immodel.w,[],true);
            end
            if (nargin < 3 || isempty(usefit))
                usefit = true;
            end
            
            locvec = zeros([3 0]);
            dlocvec =  zeros([3 0]);
            pp = [];
            dw = [];
            for j = 1:length(nra.ploc)
                valid = nra.ploc{j}(1,:) >= min(nra.xrange) & nra.ploc{j}(1,:) <= max(nra.xrange) & nra.ploc{j}(2,:) >= min(nra.yrange) & nra.ploc{j}(2,:) <= max(nra.yrange) & cos(nra.pphase{j}) >= min(nra.czrange) & cos(nra.pphase{j}) <= max(nra.czrange);
                pl = nra.ploc{j}(:,valid);
                pp = [pp nra.pphase{j}(valid)]; %#ok<AGROW>
                pl(3,:) = nra.z_scale*cos(nra.pphase{j}(valid));
                pl(4,:) = 1;
                pl = nra.fixedAffine{j}*pl;
                locvec = [locvec pl]; %#ok<AGROW>
                if (nra.hasdwell)
                    valid = nra.dloc{j}(1,:) >= min(nra.xrange) & nra.dloc{j}(1,:) <= max(nra.xrange) & nra.dloc{j}(2,:) >= min(nra.yrange) & nra.dloc{j}(2,:) <= max(nra.yrange) & cos(nra.dphase{j}) >= min(nra.czrange) & cos(nra.dphase{j}) <= max(nra.czrange);
                    valid = valid&(nra.dloc{j}(1,:).^2 + nra.dloc{j}(2,:).^2 < nra.rsqrange);
                    dl = nra.dloc{j}(:,valid);
                    dw = [dw nra.dweight{j}(valid)]; %#ok<AGROW>
                    dl(3,:) = nra.z_scale*cos(nra.dphase{j}(valid));
                    dl(4,:) = 1;
                    dl = nra.fixedAffine{j}*dl;
                    dlocvec = [dlocvec dl]; %#ok<AGROW>
                end
            end
            
            if (usefit)
                locvec = nra.getTransform(locvec);
                dlocvec = nra.getTransform(dlocvec);
            end
            counts = fa.accumVals(fa.depthCorrection(pp), locvec(1,:), locvec(2,:), locvec(3,:));
            if (nargout > 1)
                dwell = fa.accumVals(dw, dlocvec(1,:), dlocvec(2,:), dlocvec(3,:));
            end
            if nargout>2
                photonloc = locvec;
                photonphase = pp;
                dwellloc = dlocvec;
                dwellweight = dw;
            end
        end
                
        function nra = binPhotons (nra, nbins)
            % nra = binPhotons (nra, nbins)
            %
            % bins photons and dwell locations for use by fitter

            fa = nra.getFastAssembler(nbins, nra.fixedAffine{1});
            fadwell = nra.getFastAssembler(nbins/2, nra.fixedAffine{1}); %put dwell in larger bins to speed up computation

            nra.binnedRaw.im = fa.binPhotons([],[],[]); %creates all 0s
            nra.binnedRaw.dim = fadwell.binPhotons([],[],[]); %creates all 0s
            
            [nra.binnedRaw.x, nra.binnedRaw.y, nra.binnedRaw.z] = fa.getBinCenters();
            [nra.binnedRaw.dwx, nra.binnedRaw.dwy, nra.binnedRaw.dwz] = fadwell.getBinCenters();

            for j = 1:length(nra.ploc) 
                valid = nra.ploc{j}(1,:) >= min(nra.xrange) & nra.ploc{j}(1,:) <= max(nra.xrange) & nra.ploc{j}(2,:) >= min(nra.yrange) & nra.ploc{j}(2,:) <= max(nra.yrange) & cos(nra.pphase{j}) >= min(nra.czrange) & cos(nra.pphase{j}) <= max(nra.czrange);
                valid = valid&(nra.ploc{j}(1,:).^2 + nra.ploc{j}(2,:).^2 < nra.rsqrange);
                nra.binnedRaw.im = nra.binnedRaw.im + fa.binPhotons(nra.ploc{j}(1,valid), nra.ploc{j}(2,valid), nra.pphase{j}(valid),[],nra.fixedAffine{j});

                if (nra.hasdwell)
                    valid = nra.dloc{j}(1,:) >= min(nra.xrange) & nra.dloc{j}(1,:) <= max(nra.xrange) & nra.dloc{j}(2,:) >= min(nra.yrange) & nra.dloc{j}(2,:) <= max(nra.yrange) & cos(nra.dphase{j}) >= min(nra.czrange) & cos(nra.dphase{j}) <= max(nra.czrange);
                    valid = valid&(nra.dloc{j}(1,:).^2 + nra.dloc{j}(2,:).^2 < nra.rsqrange);
    
                    nra.binnedRaw.dim = nra.binnedRaw.dim + fadwell.binWithWeightAndAffine (nra.dloc{j}(1,valid), nra.dloc{j}(2,valid), nra.dphase{j}(valid), nra.dweight{j},nra.fixedAffine{j});
    
                
                   
                end
            end
            
            
            w = nra.binnedRaw.im;
            [x,y,z] = ndgrid(nra.binnedRaw.x, nra.binnedRaw.y, nra.binnedRaw.z);
            
            w = w(:);
            nra.plocvec = [x(w>0)';y(w>0)';z(w>0)']; nra.plocvec(4,:) = 1;
            nra.pweightvec = w(w>0)';

            [x,y,z] = ndgrid(nra.binnedRaw.dwx, nra.binnedRaw.dwy, nra.binnedRaw.dwz);
            w = nra.binnedRaw.dim(:);
            nra.dlocvec = [x(w>0)';y(w>0)';z(w>0)']; nra.dlocvec(4,:) = 1;
            nra.dweightvec = w(w>0)';

            
            
        end
        
        
        function [val, grad] = eval(nra, pvec)
            %[val, grad] = eval(nra, pvec)
            %
            % part of objective function for fitter
            % calculates image energy and gradients 
            % does not calculate continuity terms
            %
                       
            if (isempty(nra.pweightvec) || sum(nra.pweightvec) == 0 || ~isfinite(sum(nra.pweightvec)))
                val = 0;
                grad = 0*pvec;
                return;
            end
            
            eval2D = false;%numel(pvec) <= 12;
            
            if (nra.invertEval)
                m = -1;
            else
                m = 1;
            end
            if (nra.normByPhotons)
                m = m/sum(nra.pweightvec);
            end

            ppts = nra.plocvec(1:3,:);
            delta = 0*ppts;

            s = nra.vecToDisplacmentStructs(pvec);

            for j = 1:length(s)
                delta = delta + NonRigidAlignerBSpline.calculateDisplacement(s(j).dispmat,  nra.pB{j});
            end
            ppts = ppts+delta;

            if (eval2D)
                f = nra.immodel.F_2D(ppts(1,:),ppts(2,:));
            else
                f = nra.immodel.F(ppts(1,:),ppts(2,:),ppts(3,:));
            end
           
            if (nra.hasdwell)
                dpts = nra.dlocvec(1:3,:);
                delta = 0*dpts;
                for j = 1:length(s)
                    delta = delta + NonRigidAlignerBSpline.calculateDisplacement(s(j).dispmat,  nra.dB{j});
                end
                 if (eval2D)
                    dwellf = nra.immodel.F_2D(dpts(1,:),dpts(2,:));
                else
                    dwellf = nra.immodel.F(dpts(1,:),dpts(2,:),dpts(3,:));
                end
                val = m * (sum(nra.pweightvec.*log(f),'all','omitnan') - sum(nra.dweightvec.*dwellf,'all','omitnan'));

            else
                val = m * sum(nra.pweightvec.*f);
            end

            spring_gs = s;
            for j = 1:length(s)          
                [au,ag] = NonRigidAlignerBSpline.calculateAnchorEnergyAndGradient(s(j).dispmat, nra.anchorlength);
                [springu,springg] = NonRigidAlignerBSpline.calculateSpringEnergyAndGradient(s(j).dispmat, s(j).delta);
                val = val + (au + springu)*nra.geometry_weight;
                spring_gs(j).dispmat = (ag + springg)*nra.geometry_weight;
            end

            
            if (nargout < 2)
                return;
            end
            if (nra.hasdwell)
                if (eval2D)
                    gu = m*(nra.pweightvec.*nra.immodel.dFdu_2D(ppts(1,:),ppts(2,:))./f);
                    gv = m*(nra.pweightvec.*nra.immodel.dFdv_2D(ppts(1,:),ppts(2,:))./f);
                    gw = 0*gv;
                else
                    gu = m*(nra.pweightvec.*nra.immodel.dFdu(ppts(1,:),ppts(2,:),ppts(3,:))./f);
                    gv = m*(nra.pweightvec.*nra.immodel.dFdv(ppts(1,:),ppts(2,:),ppts(3,:))./f);
                    gw = m*(nra.pweightvec.*nra.immodel.dFdw(ppts(1,:),ppts(2,:),ppts(3,:))./f);
                end
                pimgrad = [gu(:)';gv(:)';gw(:)'];
                if (eval2D)
                    gu = m*(nra.dweightvec.*nra.immodel.dFdu_2D(dpts(1,:),dpts(2,:))); 
                    gv = m*(nra.dweightvec.*nra.immodel.dFdv_2D(dpts(1,:),dpts(2,:)));
                    gw = 0*gv;
                else
                    gu = m*(nra.dweightvec.*nra.immodel.dFdu(dpts(1,:),dpts(2,:),dpts(3,:)));
                    gv = m*(nra.dweightvec.*nra.immodel.dFdv(dpts(1,:),dpts(2,:),dpts(3,:)));
                    gw = m*(nra.dweightvec.*nra.immodel.dFdw(dpts(1,:),dpts(2,:),dpts(3,:)));
                end
                dimgrad = [gu(:)';gv(:)';gw(:)'];

            else
                if (eval2D)
                    gu = m*(nra.pweightvec.*nra.immodel.dFdu_2D(ppts(1,:),ppts(2,:))); 
                    gv = m*(nra.pweightvec.*nra.immodel.dFdv_2D(ppts(1,:),ppts(2,:)));
                    gw = 0*gv;
                else
                    gu = m*(nra.pweightvec.*nra.immodel.dFdu(ppts(1,:),ppts(2,:),ppts(3,:)));
                    gv = m*(nra.pweightvec.*nra.immodel.dFdv(ppts(1,:),ppts(2,:),ppts(3,:)));
                    gw = m*(nra.pweightvec.*nra.immodel.dFdw(ppts(1,:),ppts(2,:),ppts(3,:)));
                end
                pimgrad = [gu(:)';gv(:)';gw(:)'];
            end

            im_gs = s;

            for j = 1:length(s)
                grad = NonRigidAlignerBSpline.calculateGradient(pimgrad, nra.pB{j});
                if (nra.hasdwell)
                    grad = grad - NonRigidAlignerBSpline.calculateGradient(dimgrad, nra.dB{j});
                end
                im_gs(j).dispmat = grad;
            end

            grad = nra.displacementStructsToVec(im_gs) + nra.displacementStructsToVec(spring_gs);

        end
        

        %sequential fitting appears slower than all at once with 2 levels
        %(1 example) -- but could be worth playing with
        function [nra, pvec] = fitNonRigidSequential(nra, binsize)
            disp('warning: this may be slower than fitNonRigid - make sure you know what you''re doing')
             if (nargin < 2)
                binsize = []; 
             end
             pvec = [];
             for j = 1:length(nra.npts)
                 nra2 = nra;
                 nra2.npts = nra2.npts(1:j);
                 [nra2, pvec] = nra2.fitNonRigid(pvec, binsize);
             end
             nra = nra2;
        end

        function [nra, pvec] = fitNonRigid (nra, pvec_init, binsize)
            % [nra, pvec] = fitNonRigid (nra, pvec_init, binsize)
            % fits to multi-level b-spline free form deformation
            % adjust npts, geometry_weight, anchorlength to control fit

            if (nargin < 3)
                binsize = []; 
            end

            nra = nra.binPhotonsByBinsize(binsize);
            nra = nra.adjustDwellToMatchRate(false);
            nra = nra.calculateBSPlineParameters();

            if (nargin < 2 || isempty(pvec_init))
                 pvec_init = nra.displacementStructsToVec(nra.vecToDisplacmentStructs()); %zeros
%                 pvec_init = zeros([3,nra.npts]);
            end
            if (isstruct(pvec_init))
                for j = 1:length(pvec_init)
                    if any(pvec_init(j).npts ~= nra.npts{j})
                        error ('initial guess structure does not match nra structure')
                    end
                end
                %initial guess struct should have initial guesses for earlier layers,
                %but does not have to have all layers \
                % (i.e. if len(nra.npts) is 3, then initial guess can have 2
                % elements and be ok

                pvi = nra.displacementStructsToVec(pvec_init);
     
                pvec_init = nra.displacementStructsToVec(nra.vecToDisplacmentStructs()); %zeros
                pvec_init(1:min(length(pvi),length(pvec_init))) = pvi(1:min(length(pvi),length(pvec_init))); %allow shortening and lengthening init guess
            end
            
            op = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Display', 'notify', 'MaxFunctionEvaluations', 3e4, 'StepTolerance', 1e-8);

            [v,g] = nra.eval(pvec_init);


            [pvec, val] = fmincon(@(x) nra.eval(x), pvec_init, [], [], [], [], nra.lb, nra.ub, [], op);
            nra.fitParameters = nra.vecToDisplacmentStructs(pvec);
            pvec = nra.fitParameters;

            if (nra.normByPhotons)
                nra.funPerPhoton = val;
            else
                nra.funPerPhoton = val/(sum(nra.pweightvec));
            end
        end

        function nra = adjustDwellToMatchRate(nra, eval2D)
            if (~nra.hasdwell)
                return;
            end
            dpts = nra.dlocvec;
            
            if (eval2D)
                dwellf = nra.immodel.F_2D(dpts(1,:),dpts(2,:));
            else
                dwellf = nra.immodel.F(dpts(1,:),dpts(2,:),dpts(3,:));
            end
            predicted_counts = sum(nra.dweightvec.*dwellf,'all','omitnan');
            total_counts = sum(nra.pweightvec,'all','omitnan');
            correction = total_counts/predicted_counts;
           % disp(['correction = ' num2str(correction)]);
            nra.dweightvec = nra.dweightvec*correction;
            nra.dweight_correction = correction;
        end
        
        function nra = setFitParameters(nra,pvec)
            nra.fitParameters = pvec;
            nra.npts = cell(size(pvec));
            nra.x0 = cell(size(pvec));
            nra.delta = cell(size(pvec));
            for i=1:length(nra.x0)
                nra.npts{i} = pvec(i).npts;
                nra.x0{i} = pvec(i).x0;
                nra.delta{i} = pvec(i).delta;
            end
            % ignores dweight_correction
        end
    end
    
    methods(Static)
        
%         function b = truncatedBernsteinPolynomial(n,k,s)
%             if (size(s,1) > 1 && size(s,2) > 2)                    
%                 b = 0*s;
%                 for j = 1:size(s,1)
%                     b(j,:) = truncatedBernsteinPolynomial(n,k,s(j,:));
%                 end
%                 return;
%             end
%             b = nchoosek(n,k)*s.^(k).*(1-s).^(n-k);
%             b(s < 0) = 0;
%             b(s >= 1) = 0;
%         end

        function b = chanCubicPolynomial(s)
            if (size(s,1) > 1 && size(s,2) > 2)                    
                b = 0*s;
                for j = 1:size(s,1)
                    b(j,:) = chanCubicPolynomial(s(j,:));
                end
                return;
            end
            b = (2+s).^3/6 .*(s >= -2 & s < -1) + (2/3 - s.^2 - s.^3/2).*(s >= -1 & s < 0) + (2/3 - s.^2 + s.^3/2).*(s >= 0 & s < 1) + (2-s).^3/6 .*(s >= 1 & s < 2);
        end

        function d = calculateDisplacement(disp_mat, B)
            %d = calculateDisplacement(disp_mat, B)
            %cpt_mat is 3 x nx x ny x nz
            %B is {nx x N ,ny x N,nz x N}
            %d is 3 x N
            
            d = zeros(size(disp_mat, 1), size(B{1},2));
            npts = [size(B{1},1), size(B{2}, 1), size(B{3},1)];
            for i = 1:npts(1)
                for j = 1:npts(2)
                    for k = 1:npts(3)
                        d = d+disp_mat(:,i,j,k).*B{1}(i,:).*B{2}(j,:).*B{3}(k,:);
                    end
                end
            end
        end

        function g = calculateGradient(df, B)
            %g = calculateGradient(df, B)
            %df is 3 x N
            %B is {nx x N ,ny x N,nz x N}
            %g is 3 x nx x ny x nz
            
 
            df(~isfinite(df)) = 0;
            npts = [size(B{1},1), size(B{2}, 1), size(B{3},1)];
            g = zeros([size(df,1) npts]);
            for i = 1:npts(1)
                for j = 1:npts(2)
                    for k = 1:npts(3)
                        g(:,i,j,k) = df*(B{1}(i,:).*B{2}(j,:).*B{3}(k,:))';
                    end
                end
            end
            
        end

        function [U,G] = calculateSpringEnergyAndGradient(M, ptspacing)
            %M is 3 x nx x ny x nz M(:,i,j,k) is (dx,dy,dz) @ point ijk
            %ptspacing is 3 element vector giving deltaX, deltaY, deltaZ
            %spacing between points
            %             if (nargin < 2)
            %                 ptspacing = [1 1 1];
            %             end
        
            ptspacing = ptspacing(:);
        
            DMI = diff(M(1,:,:,:),[],2)/ptspacing(1);
            DMJ = diff(M(2,:,:,:),[],3)/ptspacing(2);
            DMK = diff(M(3,:,:,:),[],4)/ptspacing(3);
        
            DMIJ = (M(1:2,2:end,2:end,:) - M(1:2,1:end-1,1:end-1,:)).*ptspacing(1:2)./sum(ptspacing(1:2).^2);
            DMIK = (M([1 3],2:end,:,2:end) - M([1 3],1:end-1,:,1:end-1)).*ptspacing([1 3])./sum(ptspacing([1 3]).^2);
            DMJK = (M(2:3,:,2:end,2:end) - M(2:3,:,1:end-1,1:end-1)).*ptspacing(2:3)./sum(ptspacing(2:3).^2);
        
            DMIJK = (M(:,2:end,2:end,2:end) - M(:,1:end-1,1:end-1,1:end-1)).*ptspacing./sum(ptspacing.^2);
        
            U = 0.5*(sum(DMI.^2,'all') + sum(DMJ.^2,'all') + sum(DMK.^2,'all') + (sum(sum(DMIJ,1).^2,'all') + sum(sum(DMIK,1).^2,'all') + sum(sum(DMJK,1).^2,'all')) +  sum(sum(DMIJK,1).^2, 'all'));
        
            DDMI = -padarray(DMI,[0 1 0 0],0,'post') + padarray(DMI,[0 1 0 0],0,'pre');
            DDMJ = -padarray(DMJ,[0 0 1 0],0,'post') + padarray(DMJ,[0 0 1 0],0,'pre');
            DDMK = -padarray(DMK,[0 0 0 1],0,'post') + padarray(DMK,[0 0 0 1],0,'pre');
        
        
            DDMIJ = -padarray(DMIJ,[0 1 1 0],0,'post') + padarray(DMIJ,[0 1 1 0],0,'pre');
            DDMIK = -padarray(DMIK,[0 1 0 1],0,'post') + padarray(DMIK,[0 1 0 1],0,'pre');
            DDMJK = -padarray(DMJK,[0 0 1 1],0,'post') + padarray(DMJK,[0 0 1 1],0,'pre');
        
            DDMIJK =  -padarray(DMIJK,[0 1 1 1],0,'post') + padarray(DMIJK,[0 1 1 1],0,'pre');
        
            GX = DDMI/ptspacing(1) + sum(DDMIJ,1)*ptspacing(1)/sum(ptspacing(1:2).^2) + sum(DDMIK,1)*ptspacing(1)/sum(ptspacing([1 3]).^2) + sum(DDMIJK,1)*ptspacing(1)/sum(ptspacing.^2);
            GY = DDMJ/ptspacing(2) + sum(DDMIJ,1)*ptspacing(2)/sum(ptspacing(1:2).^2) + sum(DDMJK,1)*ptspacing(2)/sum(ptspacing(2:3).^2)  + sum(DDMIJK,1)*ptspacing(2)/sum(ptspacing.^2);
            GZ = DDMK/ptspacing(3) + sum(DDMIK,1)*ptspacing(3)/sum(ptspacing([1 3]).^2) + sum(DDMJK,1)*ptspacing(3)/sum(ptspacing(2:3).^2) + sum(DDMIJK,1)*ptspacing(3)/sum(ptspacing.^2);
        
            G = cat(1,GX,GY,GZ);
        end

        function [U,G] = calculateAnchorEnergyAndGradient(M, lengthscale)
            U = 0.5*sum(M.^2,'all')/lengthscale.^2;
            G = M/lengthscale.^2;
        end
    

        

        
        
        

    end
    
end



