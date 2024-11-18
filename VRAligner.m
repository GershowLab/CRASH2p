classdef VRAligner
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vr
        templatecoarse
        templatefine
        pvec
        metrics
        
        facoarse
        fafine
        
        fadisplay = []
        
        template_range
        aligner_range
        
        z_scale
        
        horizontal_percentile_cut = 0.95;
        rate_rescale = 1e7;
        beamnum = 2;
        beamcolor = 'r';
        
        theta_init = [];
        
        immodel_coarse
        immodel_fine
        
        rmin = 5e4; %minimum rate for template images
        
        static_scan_box = false;
        
        
        pvec_subframe;
        subframe_time_edges;
        hassubframe;
        
        nsubframe = 4;
        
        additional_rotation = [0 0 0];
        
        pvec_non_rigid
        metrics_non_rigid
        
        valid
    end
    
    properties (Transient = true)
        accumulated_images
    end
    
    
    methods
        function vra = VRAligner(vr)
            % vra = VRAligner(vr)
            % vr < Volume Reconstructor
            vra = vra.setVR(vr);
        end
        
        function vra = setVR(vra, vr)
            % vra = setVR(vra, vr)
            % vra < VRAligner
            % vr < VolumeReconstructor
            vra.vr = vr;
            if(isempty(vr.scanBox))
                vra.vr = vra.vr.calculateScanBox(vra.beamnum, vra.beamcolor);
            end
            
            vra.static_scan_box = diff(prctile(vra.vr.scanBox.center(1,:), [1 99])) < 2 && diff(prctile(vra.vr.scanBox.center(2,:), [1 99])) < 2;
            
            is = vr.settings(end).imagingSettings;
            sr = max(is.xScanRange_Microns, is.yScanRange_Microns);
            vra.z_scale = vr.settings(end).analogScalingFactors(2).zTagMicronsPerTick * 160E6 / (2*pi*vr.settings(end).analogScalingFactors(2).tagLensFrequency);
            
            
%             if (vra.static_scan_box) 
%                 vra.template_range.x = sr/sqrt(2) * [-1 1];
%                 vra.template_range.y = sr/sqrt(2) * [-1 1];
%             else

            vra.template_range.x = sqrt(max(vra.vr.scanBox.globalBox(1:2).^2) + max(vra.vr.scanBox.globalBox(3:4).^2))*[-1 1];
            vra.template_range.y = vra.template_range.x;
            vra.template_range.z = [-1.5 1.5]*vra.z_scale;
            
            vra.aligner_range.edgeClearance = sr/20;
            
            %TODO: need to handle rsq in case where scan box is moving ---
            if (vra.static_scan_box)
                vra.aligner_range.rsq = (0.45*min(is.xScanRange_Microns, is.yScanRange_Microns)).^2;
            else
                vra.aligner_range.rsq = Inf;
            end
            vra.aligner_range.cz = [-0.8 0.8];
            
            vra = vra.createAssemblers();
            vra = vra.resetResults();
            
        end
        
        function tr = trackingTimeRange(vra)
            tx = (vra.vr.tsr.tx(vra.vr.tsr.tracker.tracking));   
            tr = [max(vra.vr.frame.edges(1),min(tx)) min(vra.vr.frame.edges(end-1),max(tx))];
        end
        
        function vra = resetResults(vra)
            % vra = resetResults(vra)
            % clears/initializes to 0 metrics and fit vectors
            vra.pvec = zeros(length(vra.vr.frame.edges)-1, 6);
            vra.hassubframe = false(size(vra.pvec,1), 1);
            vra.subframe_time_edges = zeros(size(vra.pvec,1), vra.nsubframe+1);
            vra.pvec_subframe = zeros([size(vra.pvec,1) vra.nsubframe size(vra.pvec,2)]);
            
            metric.fval_per_photon = 0;
            metric.fval_per_photon_in_radius = 0;
            metric.mi = 0;
            metric.mi2D = 0;
            metric.fval_sequence = 0;
            vra.metrics = repmat(metric, [length(vra.vr.frame.edges)-1 1]);
            
            vra = vra.resetAccumulators();
            vra.theta_init = zeros(1, length(vra.vr.frame.edges)-1); % not used?
            
            vra.pvec_non_rigid = cell(length(vra.vr.frame.edges)-1,1);
            metric = [];
            metric.fval_per_photon = 0;
            metric.mi = 0;
            metric.mi2D = 0;
            vra.metrics_non_rigid = repmat(metric, [length(vra.vr.frame.edges)-1 1]);
            
            vra.valid = false(1,length(vra.vr.frame.edges)-1);
         end
        
        function vra = createAssemblers(vra, bincoarse, binfine)
            %  vra = createAssemblers(vra, bincoarse, binfine)
            %  bincoarse - 1x3 [x y z] bin size for first coarse pass [2.5 2.5 4]
            %  binfine - 1x3 [x y z] bin size for second finer pass [1 1 1.5]
            if (nargin < 2 || isempty(bincoarse))
                bincoarse = [2.5 2.5 4];
            end
            if (nargin < 3 || isempty(binfine))
                binfine = [1 1 1.5];
            end
            if (numel(bincoarse) < 3)
                bincoarse = bincoarse*[1 1 1];
            end
            if (numel(binfine) < 3)
                binfine = binfine*[1 1 1];
            end
            
            xe = linspace(vra.template_range.x(1), vra.template_range.x(2), ceil(diff(vra.template_range.x)/bincoarse(1) + 1));
            ye = linspace(vra.template_range.y(1), vra.template_range.y(2), ceil(diff(vra.template_range.y)/bincoarse(2) + 1));
            ze = linspace(vra.template_range.z(1), vra.template_range.z(2), ceil(diff(vra.template_range.z)/bincoarse(3) + 1));
            vra.facoarse = FastAssembler(xe,ye,ze,[],false, vra.z_scale);
            
            xe = linspace(vra.template_range.x(1), vra.template_range.x(2), ceil(diff(vra.template_range.x)/binfine(1) + 1));
            ye = linspace(vra.template_range.y(1), vra.template_range.y(2), ceil(diff(vra.template_range.y)/binfine(2) + 1));
            ze = linspace(vra.template_range.z(1), vra.template_range.z(2), ceil(diff(vra.template_range.z)/binfine(3) + 1));
            vra.fafine = FastAssembler(xe,ye,ze,[],false, vra.z_scale);           
            
        end
        
        function vra = createDisplayAssembler(vra, buffer, voxelsize)
            existsAndDefault('buffer', [10 0.1]);
            existsAndDefault('voxelsize', vra.vr.settings(end).imagingSettings.pixelSize_Microns);
            im = vra.templatefine.F_2D;
            xx = im(im > vra.rmin/vra.rate_rescale);
            thresh = percentile(xx, vra.horizontal_percentile_cut);
            inds = im > thresh;
            [vv,uu] = meshgrid(vra.templatefine.v, vra.templatefine.u);
            x = uu(inds);
            y = vv(inds);
            
            xrange = [min(x) max(x)];
            yrange = [min(y) max(y)];
            
            lgdim = max(diff(xrange),diff(yrange));
            
            if (numel(buffer) == 1)
                if (buffer < 1) %assume fraction
                    buf = buffer*lgdim;
                else
                    buf = buffer;
                end
            else
                buf = max(buffer(1), buffer(2)*lgdim);
            end
            
            if length(voxelsize)==1
                voxelsize = repmat(voxelsize,1,3);
            end
            
            xe = (floor((xrange(1)-buf)/voxelsize(1)-0.5):ceil((xrange(2)+buf)/voxelsize(1)+0.5))*voxelsize(1);
            ye = (floor((yrange(1)-buf)/voxelsize(2)-0.5):ceil((yrange(2)+buf)/voxelsize(2)+0.5))*voxelsize(2);
            ze = vra.template_range.z(1):voxelsize(3):vra.template_range.z(2);
            vra.fadisplay = FastAssembler(xe, ye, ze, [], false, vra.z_scale);
        end
        
        function [ppos, pphase, ptime, poff] = getPhotonLocations(vra, trange, altcolor, uncorrected)
            % [ppos, pphase, ptime, poff] = getPhotonLocations(vra, trange, altcolor)
            % wrapper for VolumeReconstructor.getPhotonLocations,
            %   uses vra.beamnum, vra.beamcolor (if altcolor, uses other
            %   color
            %  ppos(3) is the value of the piezo
            %  poff(3) is the displacement of the tracked neuron from 0 /
            %  tagAmplitude and is not used for offset correction
            %
            %  if corrected (default), vr offsetCorrect and phaseCorrect
            %  are used, but scanCorrect is not; as scanCorrect moves the
            %  tracked neuron away from (0,0) [code updated 6/6/24]
            % 
            if (nargin >= 3)
                if (ischar(altcolor))
                    c = altcolor;
                else
                    if (altcolor)
                        if (vra.beamcolor == 'r')
                            c = 'g';
                        else
                            c = 'r';
                        end
                    else
                        c = vra.beamcolor;
                    end
                end
            else
                c = vra.beamcolor;
            end
            existsAndDefault('uncorrected', false);
            if (uncorrected)
                 [ppos, pphase, ptime, poff] = vra.vr.getPhotonLocations(vra.beamnum, c, trange, [], false, true, false);
            else
                %[pos, phase, time, offset] = getPhotonLocations(vr, beamnum, color, timerange, galvoDeltaT, offsetCorrect, phaseCorrect, scanCorrect)
               [ppos, pphase, ptime, poff] = vra.vr.getPhotonLocations(vra.beamnum, c, trange, [], true, true, false);
            end
        end
        
        %[pos, time, period, offset] = getDwellLocations(vr, beamnum, color, timerange, galvoDeltaT, offsetCorrect, scanCorrect)

        function [pos, time, period, offset] = getDwellLocations(vra, trange, altcolor, uncorrected)
            % [ppos, pphase, ptime, poff] = getPhotonLocations(vra, trange, altcolor)
            % wrapper for VolumeReconstructor.getDwellLocations,
            %   uses vra.beamnum, vra.beamcolor (if altcolor, uses other
            %   color)
            %
            % gets location of beam each time tag lens crosses z = 0
            % applies beam color misalignment but not offset
            % (phase = pi/4, 3 pi / 4; phase correction not used)
            % if aod, only gets when phase = pi/4
            %
            %
            % trange = time interval in which dwell location is returned returned
            % pos = x,y location of galvo (+ aod) at each dwell point
            % time = time when dwell location is recorded
            % period = period of tag lens at given dwell location (total
            % dwell time at each point is 1/2 period)
            %  pos(3) is the value of the piezo
            %  off(3) is the displacement of the tracked neuron from 0 /
            %  tagAmplitude and is not used for offset correction
            %
            %  if corrected (default), vr offsetCorrect and phaseCorrect
            %  are used, but scanCorrect is not; as scanCorrect moves the
            %  tracked neuron away from (0,0) 
            
            if (nargin >= 3)
                if (ischar(altcolor))
                    c = altcolor;
                else
                    if (altcolor)
                        if (vra.beamcolor == 'r')
                            c = 'g';
                        else
                            c = 'r';
                        end
                    else
                        c = vra.beamcolor;
                    end
                end
            else
                c = vra.beamcolor;
            end
            existsAndDefault('uncorrected', false);
            %        function [pos, time, period, offset] = getDwellLocations(vr, beamnum, color, timerange, galvoDeltaT, offsetCorrect, scanCorrect)

            if (uncorrected)
               [pos, time, period, offset] = vra.vr.getDwellLocations(vra.beamnum, c, trange, [], false, false);
            else
               [pos, time, period, offset] = vra.vr.getDwellLocations(vra.beamnum, c, trange, [], true, false);
            end
        end

        function [dpos, dphase, dtime, dwelltime, doffset] = getDwellLocationsSubsampled(vra, trange, deltaT, altcolor, uncorrected)
            % [dpos, dtime, dperiod, doffset] = getDwellLocations(vra, trange)
            % wrapper for VolumeReconstructor.getDwellLocations,
            %   uses vra.beamnum, vra.beamcolor
            if (nargin < 3 || isempty(deltaT))
                deltaT = 2.5e-8; %25 ns
            end
           if (nargin >= 4)
                if (ischar(altcolor))
                    c = altcolor;
                else
                    if (altcolor)
                        if (vra.beamcolor == 'r')
                            c = 'g';
                        else
                            c = 'r';
                        end
                    else
                        c = vra.beamcolor;
                    end
                end
            else
                c = vra.beamcolor;
           end
           existsAndDefault('uncorrected',false);
           if (uncorrected)
                [dpos, dphase, dtime, dwelltime, doffset] = vra.vr.getDwellLocationSubsampled(vra.beamnum,c, trange, deltaT,[], false, false);
           else
               [dpos, dphase, dtime, dwelltime, doffset] = vra.vr.getDwellLocationSubsampled(vra.beamnum,c, trange, deltaT,[], true, true);
           end

        end
        
        function vra = rotateTemplate(vra, theta_deg)      
            vra.templatecoarse.F_2D = imrotate(vra.templatecoarse.F_2D, theta_deg, 'bilinear', 'crop');
            mask = imrotate(true(size(vra.templatecoarse.F_2D)), theta_deg, 'nearest', 'crop');
            vra.templatecoarse.F_2D(~mask) = vra.rmin/vra.rate_rescale;
            
            vra.templatefine.F_2D = imrotate(vra.templatefine.F_2D, theta_deg, 'bilinear', 'crop');
            mask = imrotate(true(size(vra.templatefine.F_2D)), theta_deg, 'nearest', 'crop');
            vra.templatefine.F_2D(~mask) = vra.rmin/vra.rate_rescale;
            
            vra.templatecoarse.F = imrotate3(vra.templatecoarse.F, theta_deg(1), [0 0 1], 'linear', 'crop','FillValues',vra.rmin/vra.rate_rescale);
            vra.templatefine.F = imrotate3(vra.templatefine.F, theta_deg(1), [0 0 1], 'linear', 'crop','FillValues',vra.rmin/vra.rate_rescale);
        end
            
        
        function [vra, bbox] = makeTemplateHorizontal(vra, updatePvecs)
            %note that template must be centered around 0,0 (axis of
            %rotation is z-axis)
            %bbox = [xmin xmax; ymin ymax]
            
            if (abs(median(vra.templatefine.u)) > abs(median(diff(vra.templatefine.u)))) || (abs(median(vra.templatefine.v)) > abs(median(diff(vra.templatefine.v))))
                warning('template should be centered at 0,0');
            end
            
             im = vra.templatefine.F_2D;
             xx = im(im > vra.rmin/vra.rate_rescale);
             thresh = percentile(xx, vra.horizontal_percentile_cut);
             inds = im > thresh;
             [vv,uu] = meshgrid(vra.templatefine.v, vra.templatefine.u);
             x = uu(inds);
             y = vv(inds);
             [rx,ry] = minboundrect(x,y);
             
             %find long edges
             [~,I] = max(diff(rx).^2 + diff(ry).^2);
             theta = mod(-atan2(ry(I+1)-ry(I), rx(I+1)-rx(I)), 2*pi);
             
             im2 = imrotate(im, theta*180/pi, 'bilinear', 'crop');
             
             inds = im2 > thresh;
             x = uu(inds);
             y = vv(inds);
             y1 = y(x > mean(xx(:)));
             y2 = y(x < mean(xx(:)));
             if (std(y1) < std(y2))
                 theta = theta + pi;
             end
             
             vra = vra.rotateTemplate(theta*180/pi);
             
             bb = [cos(theta) -sin(theta); sin(theta) cos(theta)]*[rx';ry'];
             bbox = [min(bb,[],2) max(bb,[],2)];
             
             if (nargin >=2 && ~isempty(updatePvecs) && updatePvecs)
                at = RigidAligner3D.rigidTransform(theta, 0, 0, 0, 0, 0);
                for j = 1:size(vra.pvec,1)
                    if (all(vra.pvec(j,:) == 0))
                        continue;
                    end
                    atnew = [at;0 0 0 1]*[RigidAligner3D.rigidTransform(vra.pvec(j,:)); 0 0 0 1];
                    vra.pvec(j,:) = RigidAligner3D.rigidToPvec(atnew);
                end
             end
            
             
             
        end
        
        function at = manuallyOrientTemplate(vra, im) 
            u = vra.templatefine.u;
            v = vra.templatefine.v;
            w = vra.templatefine.w;
            
            figure(1); clf(1);


            pcolor(u,v,max(im,[],3)');shading flat; axis equal
            
            
            title('Please draw TARGET horizontal line from TARGET left to TARGET right')
            r = drawline;
            leftpt = r.Position(1,:);
            rightpt = r.Position(2,:);
            theta_z = -atan2(rightpt(2)-leftpt(2), rightpt(1)-leftpt(1));
            
            at = orientPoints_xyrot(u,v,w,im, theta_z);
            
            tform = affine3d([at' [0;0;0;1]]);
            im2 = permute(imwarp(permute(im,[2 1 3]), tform,'linear', 'OutputView', affineOutputView(size(im), tform, 'BoundsStyle','CenterOutput')),[2 1 3]);

            figure(1); clf(1);
            subplot(2,2,1)
            pcolor(u,v,max(im2,[],3)');shading flat; axis equal; axis tight;
            subplot(2,2,2)
            pcolor(w,v,squeeze(max(im2,[],1)));shading flat; axis equal; axis tight;
            subplot(2,2,3)
            pcolor(u,w,squeeze(max(im2,[],2))');shading flat; axis equal; axis tight;
            
            
        end
        
        function vra = manuallyAlignTemplate(vra, template_time_range)
            inds = find(vra.vr.frame.edges >= template_time_range(1) & vra.vr.frame.edges < template_time_range(2));
            [vra.templatefine.u, vra.templatefine.v, vra.templatefine.w] = vra.fafine.getBinCenters();
            [counts_fine, dwell_fine] = vra.binFrame(inds, vra.fafine, RigidAligner3D.rigidTransform([0 0 0 0 0 0]));
            mmat = imclose(dwell_fine > 0, ones([9 9 9]));
            rv = l1spline(counts_fine./dwell_fine, 0.01, 1, 100, 1, 1e-5);
            rv = mmat.*rv ;
             
            at = manuallyOrientTemplate(vra, rv);
            pause(0.01);
            
            vra = initializeTemplates(vra,  template_time_range, at);
        end
        
        function [vra,inds] = initializeTemplates(vra, template_time_range, orient)
            % vra = initializeTemplates(vra, template_time_range)
            % creates coarse and fine templates - bin sizes set by
            %    createAssembelrs
            % template_time_range - relatively still time range to
            %    initialize template
            % template will be rotated to make major axis horizontal
            % uses vra.horizontal_percentile_cut to decide which parts of
            %   the template to use in this calculation (only voxels above this percentile count) 

            existsAndDefault('orient',true);

            mindwell = 1/VolumeReconstructor.FPGA_CLOCK;
             
            inds = find(vra.vr.frame.edges >= template_time_range(1) & vra.vr.frame.edges < template_time_range(2));
                   

%              [vra.templatecoarse.u, vra.templatecoarse.v, vra.templatecoarse.w] = vra.facoarse.getBinCenters();
%              [counts_coarse, dwell_coarse] = vra.binFrame(inds, vra.fa
%              mmat = imclose(dwell_coarse > 0, ones([3 3 3]));
%              rv = l1spline(counts_coarse./dwell_coarse, 0.01, 1, 100, 1, 1e-5);
%              rv = mmat.*rv ;
             
             
%              rv(rv < vra.rmin) = vra.rmin;
%              vra.templatecoarse.F =rv/vra.rate_rescale;
%              vra.templatecoarse.F_2D = max(vra.templatecoarse.F,[],3);
             
             %make fine template
             [vra.templatefine.u, vra.templatefine.v, vra.templatefine.w] = vra.fafine.getBinCenters();
             [counts_fine, dwell_fine] = vra.binFrame(inds, vra.fafine);
%              counts_fine = vra.fafine.binPhotons(ppos(1,:), ppos(2,:), pphase,[],at);
%              dwell_fine = vra.fafine.binDwell(dpos(1,:), dpos(2,:), dperiod/2,[],at);
%              dwell_fine(dwell_fine < mindwell) = 0;
             mmat = imclose(dwell_fine > 0, ones([9 9 9]));
             rv = l1spline(counts_fine./dwell_fine, 0.01, 1, 100, 1, 1e-5);
             rv = mmat.*rv ;
             if (isempty(vra.rmin)) % rmin auto-calc moved down from above to fix tiny bug
                 vra.rmin = prctile (rv(imerode(mmat, ones([3 3 3]))), 10);
             end
             rv(rv < vra.rmin) = vra.rmin;
             vra.templatefine.F =rv/vra.rate_rescale;            
             vra.templatefine.F_2D = max(vra.templatefine.F,[],3);
             vra.templatefine.counts = counts_fine;
             vra.templatefine.dwell = dwell_fine;
             vra = makeCoarseTemplateFromFine(vra);

             if (orient)
                vra = vra.makeTemplateHorizontal();                
             end
             
             
        end
           
        function vra = setImmodels(vra)
            % vra = setImmodels(vra)
            % creates immodels from templates
            vra.immodel_coarse = Immodel(vra.templatecoarse);
            vra.immodel_fine = Immodel(vra.templatefine);
        end
        
        function vra = resetAccumulators(vra)
            vra.accumulated_images = [];

        end

        function at = getAffineTransform(vra, t)
            tx = reshape(vra.subframe_time_edges(:,1:4)',[],1);
            pv = reshape(permute(vra.pvec_subframe,[2 1 3]),[],6);
            pvi = interp1(tx(tx > 0), pv(tx>0,:), t, 'previous');
            for j = 1:length(t)
                at(:,:,j) = RigidAligner3D.rigidTransform(pvi(j,:));
            end

        end
        
        function [counts, dwell, im] = binSubFrame(vra, frameInd, fa, altcolor)
            % [counts, dwell, im] = binFrame(vra, ind, fa, at)
            % gets counts, dwell from frame(s) indicated by ind
            % use vra.fafine if fa not specified
            % if at not specified (usual use): aligns everything according to fit pvecs
            % if im is requested does spline interpolation/smoothing
            if (nargin < 3 || isempty(fa))
                fa = vra.fafine;
            end
            if (nargin < 4 || isempty(altcolor))
                altcolor = false;
            end
            if (~vra.hassubframe(frameInd))
                warning(['no subframe info for frame ' num2str(frameInd)])
                if (nargout < 3)
                    [counts, dwell] = vra.binFrame(frameInd, fa, [], altcolor);
                else
                    [counts, dwell, im] = vra.binFrame(frameInd, fa, [], altcolor);
                end
                return;
            end
            
            mindwell = 1/VolumeReconstructor.FPGA_CLOCK;
            [ppos, pphase, ptime] = vra.getPhotonLocations(vra.vr.frame.edges(frameInd + [0 1]), altcolor);
            [dpos, dphase, dtime, dtau] = vra.getDwellLocationsSubsampled(vra.vr.frame.edges(frameInd + [0 1]),[],altcolor);
            counts = zeros([[fa.getImSize] vra.nsubframe]);
            dwell = counts;
                
            for j = 1:vra.nsubframe
                pvec = vra.pvec_subframe(frameInd,j,:);
                if any(vra.additional_rotation)
                    rotation_angle_rad = deg2rad(vra.additional_rotation([3 2 1])); % match pvec order
                    rotation_angle_rad(3) = -rotation_angle_rad(3); % match image display axes direction
                    pvec_add = [rotation_angle_rad 0 0 0];
                    at = RigidAligner3D.rigidTransform(pvec);
                    at(4,:) = [0 0 0 1];
                    at_add = RigidAligner3D.rigidTransform(pvec_add);
                    at_add(4,:) = [0 0 0 1];
                    att = at_add*at; % order of operation: right to left
                    att = att(1:3,:);
                else
                    att = RigidAligner3D.rigidTransform(pvec);
                end
                valid = ptime >= vra.subframe_time_edges(frameInd,j) &  ptime <= vra.subframe_time_edges(frameInd,j+1);
%                 counts(:,:,:,j) =  fa.binWithWeightAndAffine(ppos(1,valid), ppos(2,valid), pphase(valid),[],att);
                counts(:,:,:,j) =  fa.binWithWeightAndAffine(ppos(1,valid), ppos(2,valid), pphase(valid),fa.depthCorrection(pphase(valid)),att);
                valid = dtime >= vra.subframe_time_edges(frameInd,j) &  dtime <= vra.subframe_time_edges(frameInd,j+1);
                dwell(:,:,:,j) = fa.binWithWeightAndAffine(dpos(1,valid), dpos(2,valid), dphase(valid),dtau(valid),att);
            end
            counts(dwell < mindwell) = 0;
            dwell(dwell < mindwell) = 0;
            if (nargout >= 3)
                im = l1spline(counts./dwell, 0.001, 1, 100, 1, 1e-5).*imclose(dwell > 0, ones([9 9 9 9]));
            end
            
        end
        
        function activity = getActivityInConvexVoi(vra, cv, timerange)
            existsAndDefault('timerange', [-Inf, Inf]);
            valid = vra.valid(:) & (vra.subframe_time_edges(:,1) > timerange(1)) & (vra.subframe_time_edges(:,end) < timerange(2));

            activity.tx = zeros([1 nnz(valid)*vra.nsubframe]);
            activity.rc = zeros([numel(cv) numel(activity.tx)]);
            activity.gc = activity.rc;
            activity.rd = activity.rc;
            activity.gd = activity.rc;

            j0 = find(valid(:),1,'first');
            for j = find(valid(:)')
                for k = 1:vra.nsubframe
                    n = (j-j0)*vra.nsubframe + k;
                    activity.tx(n) = vra.subframe_time_edges(j,k);
                    if (vra.beamcolor == 'r')
                        [rc,gc,rd,gd] = vra.integrateInConvexVoi(j,k,cv);
                    else
                        [gc,rc,gd,rd] = vra.integrateInConvexVoi(j,k,cv);
                    end
                    activity.rc(:,n) = rc(:);
                    activity.gc(:,n) = gc(:);
                    activity.rd(:,n) = rd(:);
                    activity.gd(:,n) = gd(:);                    
                end
            end

        end

        function [counts,altcounts,dwell,altdwell] = integrateInConvexVoi (vra, frameInd, subframeind, cv)
            % [counts,altcounts,dwell,altdwell] = integrateInConvexVoi (vra, frameInd, subframeind, cv)
            % gets number of photons and optionally dwell time within
            % convex hull voi(s) at specific frameInd,subframeind
            %
            % vois are transformed into the coordinates of acquisition
            %
            % note that photon and dwell locations are offset from scan
            % center by default
            % 
            

            %move convex hull vois from template to subframe coordinates
            for j = 1:numel(cv)
                cv(j) = cv(j).transform(RigidAligner3D.rigidTransform(vra.pvec_subframe(frameInd,subframeind,:)), true); %invert
            end

            timerange = vra.subframe_time_edges(frameInd, subframeind + [0 1]);

            [ppos, pphase] = vra.getPhotonLocations(timerange, false);
            [ppos_alt, pphase_alt] = vra.getPhotonLocations(timerange, true);
            
            %counts(:,:,:,j) =  fa.binWithWeightAndAffine(ppos(1,valid), ppos(2,valid), pphase(valid),fa.depthCorrection(pphase(valid)),att);

            ploc = [ppos(1:2,:);vra.z_scale*cos(pphase)]'; %photon locations as Nx3 vectors, offset from tracked neuron
            ploc_alt = [ppos_alt(1:2,:);vra.z_scale*cos(pphase_alt)]'; %photon locations as Nx3 vectors, offset from tracked neuron
            
            counts = zeros(size(cv));
            altcounts = counts;
            if (nargout > 2)
                dwell = counts;
                altdwell = counts;
            end
            for j = 1:numel(cv)
                valid = cv(j).inShape3D(ploc);
                counts(j) = sum(abs(sin(pphase(valid))))/0.7846; %correction used in fast assembler
                valid = cv(j).inShape3D(ploc_alt);
                altcounts(j) = sum(abs(sin(pphase_alt(valid))))/0.7846; %correction used in fast assembler
            end

            if (nargout < 3)
                return;
            end

           [pos, ~, period] = vra.getDwellLocations(timerange, false);
           [pos_alt, ~, period_alt] = vra.getDwellLocations(timerange, true);

           tau_per_um = period/(4*vra.z_scale);
           tau_per_um_alt = period_alt/(4*vra.z_scale);
           
           for j = 1:numel(cv)
               %h = cv(j).getZHeight(pos(1:2,:)');
               h = cv(j).height(pos(1,:), pos(2,:));
               dwell(j) = sum(tau_per_um(:).*h(:));
               %h = cv(j).getZHeight(pos_alt(1:2,:)');
               h = cv(j).height(pos_alt(1,:), pos_alt(2,:));
               altdwell(j) = sum(tau_per_um_alt(:).*h(:));
           end
        end
       

        function [counts, dwell, im] = binFrame(vra, ind, fa, at, altcolor, uncorrected)
            % [counts, dwell, im] = binFrame(vra, ind, fa, at)
            % gets counts, dwell from frame(s) indicated by ind
            % use vra.fafine if fa not specified
            % if at not specified (usual use): aligns everything according to fit pvecs
            % if im is requested does spline interpolation/smoothing
            if (nargin < 3 || isempty(fa))
                fa = vra.fafine;
            end
            if (nargin < 4 || isempty(at))
                at = [];
            end
            if (nargin < 5 || isempty(altcolor))
                altcolor = false;
            end
            existsAndDefault('uncorrected',false);
            if (uncorrected)
                at = eye(3,4);
            end

            mindwell = 1/VolumeReconstructor.FPGA_CLOCK;
            
            counts = 0;
            dwell = 0;
            if uncorrected || isempty(vra.pvec_non_rigid) || any(cellfun(@isempty,vra.pvec_non_rigid(ind)))
                for j = 1:length(ind)
                    if (mod(j,10) == 0), disp([num2str(j) ' / ' num2str(length(ind))]); end
                    [ppos, pphase] = vra.getPhotonLocations(vra.vr.frame.edges(ind(j) + [0 1]), altcolor, uncorrected);
                    [dpos, dphase, ~, dtau] = vra.getDwellLocationsSubsampled(vra.vr.frame.edges(ind(j) + [0 1]),[],altcolor,uncorrected);
                    if (isempty(at))
                        pvec = vra.pvec(ind(j),:); %#ok<*PROPLC>
                        if any(vra.additional_rotation)
                            rotation_angle_rad = deg2rad(vra.additional_rotation([3 2 1])); % match pvec order
                            rotation_angle_rad(3) = -rotation_angle_rad(3); % match image display axes direction
                            pvec_add = [rotation_angle_rad 0 0 0];
                            at = RigidAligner3D.rigidTransform(pvec);
                            at(4,:) = [0 0 0 1];
                            at_add = RigidAligner3D.rigidTransform(pvec_add);
                            at_add(4,:) = [0 0 0 1];
                            att = at_add*at; % order of operation: right to left
                            att = att(1:3,:);
                        else
                            att = RigidAligner3D.rigidTransform(pvec);
                        end
                    else
                        if (size(at,3) >= j)
                            att = at(:,:,j);
                        else
                            att = at(:,:,1);
                        end
                    end
%                     counts = counts + fa.binWithWeightAndAffine(ppos(1,:), ppos(2,:), pphase,[],att);
                    counts = counts + fa.binWithWeightAndAffine(ppos(1,:), ppos(2,:), pphase,fa.depthCorrection(pphase),att);
                    dwell = dwell + fa.binWithWeightAndAffine(dpos(1,:), dpos(2,:), dphase,dtau,att);
                end
            else % bin with non-rigid using nra object
                for j = 1:length(ind)
                    rfc = vra.createAligner(ind(j),vra.immodel_fine,true,true,false,altcolor);
                    nra = NonRigidAlignerBSpline(rfc);
                    rfc = []; % release memory
                    nra = nra.setFitParameters(vra.pvec_non_rigid{ind(j)});
                    [ct,dw] = nra.binCountsAndDwell(fa,true);
                    counts = counts + ct;
                    dwell = dwell + dw;
                    nra = []; % release memory
                end
            end
            dwell(dwell < mindwell) = 0;
            if (nargout >= 3)
                im = l1spline(counts./dwell, 0.01, 1, 100, 1, 1e-5).*imclose(dwell > 0, ones([9 9 9]));
            end
            
        end
        
        function [counts, dwell, im] = binTime(vra, trange, fa, at, altcolor)
            % same as vra.binFrame except using arbitrary time range
            % time range can be vector or N-by-2 matrix; if matrix then
            % treat each row as a separate pair of start and end times
            if (nargin < 3 || isempty(fa))
                fa = vra.fafine;
            end
            if (nargin < 4 || isempty(at))
                at = [];
            end
            if (nargin < 5 || isempty(altcolor))
                altcolor = false;
            end
            mindwell = 1/VolumeReconstructor.FPGA_CLOCK;
            
            counts = 0;
            dwell = 0;
            for j = 1:size(trange,1)
                [ppos, pphase] = vra.getPhotonLocations(trange(j,:),altcolor);
                [dpos, dphase, ~, dtau] = vra.getDwellLocationsSubsampled(trange(j,:),[],altcolor);
                if (isempty(at))
                    if range(trange)>=(1/vra.vr.frame.framerate_Hz)
                        jj = find(vra.vr.frame.edges>=trange(j,1) & vra.vr.frame.edges<trange(j,2));
                    else
                        jj = find(vra.vr.frame.edges>=trange(j,1),1,'first');
                    end
                    pvec = vra.pvec(ind(jj),:);
                    if any(vra.additional_rotation)
                        rotation_angle_rad = deg2rad(vra.additional_rotation([3 2 1])); % match pvec order
                        rotation_angle_rad(3) = -rotation_angle_rad(3); % match image display axes direction
                        pvec(:,1:3) = pvec(:,1:3)+rotation_angle_rad;
                    end
                    att = RigidAligner3D.rigidTransform(pvec);
                else
                    if (size(at,3) >= j)
                        att = at(:,:,j);
                    else
                        att = at(:,:,1);
                    end
                end
               counts = counts + fa.binWithWeightAndAffine(ppos(1,:), ppos(2,:), pphase,[],att);
               dwell = dwell + fa.binWithWeightAndAffine(dpos(1,:), dpos(2,:), dphase,dtau,att);
                %counts = counts + fa.binPhotons(ppos(1,:), ppos(2,:), pphase,[],RigidAligner3D.rigidTransform(vra.pvec(ind(j),:)));
                %dwell = dwell + fa.binDwell(dpos(1,:), dpos(2,:), dperiod/2, [], RigidAligner3D.rigidTransform(vra.pvec(ind(j),:)));
            end
            dwell(dwell < mindwell) = 0;
            if (nargout >= 3)
                im = l1spline(counts./dwell, 0.01, 1, 100, 1, 1e-5).*imclose(dwell > 0, ones([9 9 9]));
            end
        end
        
        function [vra,inds] = updateTemplatesFromAccumulators(vra,  accum_cutoff) 
             % vra = updateTemplatesFromAccumulators(vra)
             % recalculates coarse and fine templates using summed counts
             % and dwell in accumulators
             mindwell = 1/VolumeReconstructor.FPGA_CLOCK;             
             
            
             if (nargin > 1)
                 m = [vra.accumulated_images.metrics];
                 val = [m.(accum_cutoff.field)];
%                  if (~isfield(accum_cutoff, 'value'))
%                      try
%                          accum_cutoff.value = percentile(val, accum_cutoff.percentile);
%                      catch
%                          accum_cutoff.value = mean(val) - std(val);
%                      end
%                  end
%                  valid = val >= accum_cutoff.value;
                if ~isfield(accum_cutoff,'method')
                    accum_cutoff.method = 'percentile';
                end
                switch accum_cutoff.method
                     case 'value'
                         valid = val >= accum_cutoff.value;
                     case 'percentile'
                         accum_cutoff.value = percentile(val, accum_cutoff.percentile);
                         valid = val >= accum_cutoff.value;
                     case 'percentile by chunk'
                         frameinds = [vra.accumulated_images.ind];
                         nframe = range(frameinds);
                         nchunk = ceil(nframe/accum_cutoff.chunksize);
                         chunkstart = vra.accumulated_images(1).ind;
                         chunkstop = chunkstart+accum_cutoff.chunksize-1;
                         valid = [];
                         for i=1:nchunk
                             acindsinchunk = find(frameinds>=chunkstart & frameinds<=chunkstop);
                             cutoffval = percentile(val(acindsinchunk),accum_cutoff.percentile);
                             valid = [valid val(acindsinchunk) >= cutoffval]; %#ok<AGROW>
                             chunkstart = chunkstop+1;
                             chunkstop = chunkstart+accum_cutoff.chunksize-1;
                         end
                 end
             else
                 valid = true(size(vra.accumulated_images));
             end
             
             coarse_counts = 0*vra.accumulated_images(1).coarse_counts;
             coarse_dwell = 0*vra.accumulated_images(1).coarse_dwell;
             fine_counts = 0*vra.accumulated_images(1).fine_counts;
             fine_dwell = 0*vra.accumulated_images(1).fine_dwell;
             
             inds = find(valid);
             for i = (inds(:).')
                 coarse_counts = coarse_counts + vra.accumulated_images(i).coarse_counts;
                 coarse_dwell =  coarse_dwell + vra.accumulated_images(i).coarse_dwell;
                 fine_counts = fine_counts + vra.accumulated_images(i).fine_counts;
                 fine_dwell = fine_dwell + vra.accumulated_images(i).fine_dwell;
             end
            
%              [vra.templatecoarse.u, vra.templatecoarse.v, vra.templatecoarse.w] = vra.facoarse.getBinCenters();
%              coarse_dwell(coarse_dwell < mindwell) = 0;
%              mmat = imclose(coarse_dwell > 0, ones([3 3 3]));
%              rr = coarse_counts./coarse_dwell;
%              
%              rv = l1spline(coarse_counts./coarse_dwell, 0.01, 1, 100, 1, 1e-5);
%              rv = mmat.*rv ;
%              rv(rv < vra.rmin) = vra.rmin;
%              vra.templatecoarse.F =rv/vra.rate_rescale;
%              vra.templatecoarse.F_2D = max(vra.templatecoarse.F,[],3);
             
             %make fine template
             [vra.templatefine.u, vra.templatefine.v, vra.templatefine.w] = vra.fafine.getBinCenters();
            
             fine_dwell(fine_dwell < mindwell) = 0;
             mmat = imclose(fine_dwell > 0, ones([9 9 9]));
             rv = l1spline(fine_counts./fine_dwell, 0.01, 1, 100, 1, 1e-5);
             rv = mmat.*rv ;
             rv(rv < vra.rmin) = vra.rmin;
             vra.templatefine.F =rv/vra.rate_rescale;
             vra.templatefine.F_2D = mean(vra.templatefine.F,3); %mean here and max in coarse is intentional 
             vra.templatefine.counts = fine_counts;
             vra.templatefine.dwell = fine_dwell;
             
             %resample fine template to make coarse template
             vra = vra.makeCoarseTemplateFromFine();
%              sigma = [vra.facoarse.dx./vra.fafine.dx, vra.facoarse.dy./vra.fafine.dy, vra.facoarse.dz./vra.fafine.dz];
%              [vra.templatecoarse.u, vra.templatecoarse.v, vra.templatecoarse.w] = vra.facoarse.getBinCenters();
%              Fblur =  imblur(vra.templatefine.F,sigma);
%              [uu,vv,ww] = ndgrid(vra.templatecoarse.u, vra.templatecoarse.v, vra.templatecoarse.w);
%              vra.templatecoarse.F = interpn(vra.templatefine.u, vra.templatefine.v, vra.templatefine.w, Fblur,uu,vv,ww);
%              vra.templatecoarse.F_2D = max(vra.templatecoarse.F,[],3);   
             
        end
        
        function vra = makeCoarseTemplateFromFine(vra)
             sigma = [vra.facoarse.dx./vra.fafine.dx, vra.facoarse.dy./vra.fafine.dy, vra.facoarse.dz./vra.fafine.dz];
             [vra.templatecoarse.u, vra.templatecoarse.v, vra.templatecoarse.w] = vra.facoarse.getBinCenters();
             Fblur =  imblur(vra.templatefine.F,sigma);
             [uu,vv,ww] = ndgrid(vra.templatecoarse.u, vra.templatecoarse.v, vra.templatecoarse.w);
             vra.templatecoarse.F = interpn(vra.templatefine.u, vra.templatefine.v, vra.templatefine.w, Fblur,uu,vv,ww);
             vra.templatecoarse.F_2D = max(vra.templatecoarse.F,[],3);   
        end

        %function vra = updateTemplateFromFitResults (vra, metric_cutoff)

        function vra = updateTemplateFromFitResults (vra, metric_cutoff)
            mindwell = 1/VolumeReconstructor.FPGA_CLOCK;             
            m = vra.metrics;
            val = [m.(metric_cutoff.field)];
            if (~isfield(metric_cutoff, 'value'))
                try
                    metric_cutoff.value = percentile(val, metric_cutoff.percentile);
                catch
                    metric_cutoff.value = mean(val) - std(val);
                end
            end
            valid = val >= metric_cutoff.value;
            inds = find(valid);
            
            disp(['length ind = ' num2str(length(inds))]);
            
            ts = tic;
            disp('starting coarse binning');
%             [vra.templatecoarse.u, vra.templatecoarse.v, vra.templatecoarse.w] = vra.facoarse.getBinCenters();
%             [coarse_counts, coarse_dwell] = vra.binFrame(inds, vra.facoarse);
%             coarse_dwell(coarse_dwell < mindwell) = 0;
%             mmat = imclose(coarse_dwell > 0, ones([3 3 3]));
%             rv = l1spline(coarse_counts./coarse_dwell, 0.01, 1, 100, 1, 1e-5);
%             rv = mmat.*rv ;
%             rv(rv < vra.rmin) = vra.rmin;
%             vra.templatecoarse.F =rv/vra.rate_rescale;
%             vra.templatecoarse.F_2D = max(vra.templatecoarse.F,[],3);
%             toc(ts)
%             disp('starting fine binning');
            
            %make fine template
            [vra.templatefine.u, vra.templatefine.v, vra.templatefine.w] = vra.fafine.getBinCenters();
            [fine_counts, fine_dwell] = vra.binFrame(inds, vra.fafine);
            fine_dwell(fine_dwell < mindwell) = 0;
            mmat = imclose(fine_dwell > 0, ones([9 9 9]));
            rv = l1spline(fine_counts./fine_dwell, 0.01, 1, 100, 1, 1e-5);
            rv = mmat.*rv ;
            rv(rv < vra.rmin) = vra.rmin;
            vra.templatefine.F =rv/vra.rate_rescale;
            vra.templatefine.F_2D = max(vra.templatefine.F,[],3);
            vra = makeCoarseTemplateFromFine(vra);
            toc (ts);
        end
            
        function vra = initial2DAlignment(vra, inds, display)
            existsAndDefault('display', true);
            segmentWidth = 300;
            stitchHalfWidth = 30;

            theta_z = zeros(size(inds));
            fval = theta_z;
            
            vra = vra.setImmodels();
            ci = round(linspace(1,length(inds), ceil(length(inds)/segmentWidth)+1));
             
            for j = 1:(length(ci)-1)
                 vra = coarseAlignLinkedSegment(vra, inds(ci(j):ci(j+1)), display);
                 if (display)
                     figure(2); clf(2);
                     plotyy(inds(1:ci(j+1)), vra.pvec(inds(1:ci(j+1)),1)*180/pi, inds(1:ci(j+1)), vra.pvec(inds(1:ci(j+1)),[4 5])); drawnow
                 end
            end
            
            pvi = vra.pvec(inds,:);
            
            for j = 2:(length(ci)-1)
                %don't align translation
                ii = inds(max(1,ci(j)-stitchHalfWidth):min(ci(j)+stitchHalfWidth,length(inds)));
                rfc = vra.createRfcSegment(ii, vra.immodel_coarse, [-pi, -pi/4, -pi/4, -15, -15, -15], [3*pi, pi/4, pi/4, 15, 15, 15]); 
                
                rfc(1).fixedBefore = vra.pvec(ii(1)-1,:);
                rfc(end).fixedAfter = vra.pvec(ii(end)+1,:);
                
            
                rfc = rfc.fit(vra.pvec(ii,[1 4 5]));
                pvfinal = cat(1,rfc.fitParameters);            
                pvfinal(:,1) = mod(pvfinal(:,1),2*pi);
                vra.pvec(ii, :) = pvfinal;
                if (display)
                     figure(2); clf(2); subplot(2,1,1);
                     plot(inds, vra.pvec(inds,1)*180/pi,'r-', inds, pvi(:,1)*180/pi, 'b-', ii, pvfinal(:,1)*180/pi,'m-');
                     
                     subplot(2,1,2);
                      plot(inds, pvi(:,[4 5]), ii, pvfinal(:,[4 5]),'m-');
                end
                     
            end           
            
        end
        
        
        function rfc = createRfcSegment(vra, inds, immodel, lb, ub)
            existsAndDefault('immodel', vra.immodel_coarse);
            existsAndDefault('lb', []);
            existsAndDefault('ub', []);
            
            for j = 1:length(inds)
                frame_ind = inds(j);
                [ppos, pphase, ptime] = vra.getPhotonLocations(vra.vr.frame.edges(frame_ind + [0 1]));
                rfc(j) = RasterFitterChunk3D(immodel, ppos, pphase, ptime, vra.z_scale); %#ok<AGROW>
                
                rfc(j).xrange = vra.vr.scanBox.corners(1:2,frame_ind) + vra.aligner_range.edgeClearance*[1; -1]; %#ok<AGROW>
                rfc(j).yrange = vra.vr.scanBox.corners(3:4,frame_ind) + vra.aligner_range.edgeClearance*[1; -1];         %#ok<AGROW>   
                rfc(j).czrange = vra.aligner_range.cz; %#ok<AGROW>
                rfc(j).rsqrange = vra.aligner_range.rsq; %#ok<AGROW>
                rfc(j).k_weight_method = 0; %#ok<AGROW>
                
                if(~isempty(ub))
                    rfc(j).ub = ub; %#ok<AGROW>
                end
                if(~isempty(lb))
                    rfc(j).lb = lb; %#ok<AGROW>
                end
            end
        end
        
        function [vra, rfc] = fineAlignSubFrame(vra, inds)
            for j = 1:length(inds)
                try
                    rfc = vra.createRfcSegment(inds(j), vra.immodel_fine);
                    rfc = rfc.splitIntoNEqualTime(vra.nsubframe);
                    rfc = rfc.setBoundsNearFit(vra.pvec(inds(j),:));
                    if (inds(j) > 1 && any(vra.pvec(inds(j)-1,:)))
                        rfc(1).fixedBefore = vra.pvec(inds(j)-1,:);
                    end
                    if (inds(j) < size(vra.pvec,1) && any(vra.pvec(inds(j)+1,:)))
                        rfc(end).fixedAfter = vra.pvec(inds(j)+1,:);
                    end
                    pvi = repmat(vra.pvec(inds(j),:), [length(rfc) 1]);
                    rfc = rfc.fit(pvi);
                    vra.pvec_subframe(inds(j),:,:) = cat(1,rfc.fitParameters);
                    vra.subframe_time_edges(inds(j),:) = ([[rfc.startTime] rfc(end).endTime]);
                    vra.hassubframe(inds(j)) = true;
                    vra.valid(inds(j)) = true;
                catch me
                    warning('failed to register subframes for frame %d',inds(j));
                    disp(me.getReport);
                    pvec_subframe = permute(repmat(vra.pvec(inds(j),:),[1 1 4]),[1 3 2]);
                    vra.pvec_subframe(inds(j),:,:) = pvec_subframe;
                    time_edges = linspace(vra.vr.frame.edges(inds(j)),vra.vr.frame.edges(inds(j)+1),5);
                    vra.subframe_time_edges(inds(j),:) = time_edges;
                    vra.hassubframe(inds(j)) = true; % not sure if should do this
                    vra.valid(inds(j)) = false;
                end
            end
        end
        
        
        function [vra,rfc, pvfinal, pvi] = coarseAlignLinkedSegment(vra, inds, display)
            % [vra,rfc, pvfinal, pvi] = coarseAlignLinkedSegment(vra, inds, fittranslation, display)
            
            if (nargin < 3 || isempty(display))
                display = false;
            end
            
            rfc = vra.createRfcSegment(inds, vra.immodel_coarse, [-pi, -pi/4, -pi/4, -15, -15, -15], [3*pi, pi/4, pi/4, 15, 15, 15]); 
            
            %initial theta_z, no translation
            theta_z = rfc.fitInitialThetaZ([],false);
            
%             
            pvi = 0*vra.pvec(inds,:);
            pvi(:,1) = mod(theta_z,2*pi);      

            rfc = rfc.fit(pvi(:,[1 4 5]));
            pvfinal = cat(1,rfc.fitParameters);
            pvfinal(:,1) = mod(pvfinal(:,1),2*pi);
            vra.pvec(inds, :) = pvfinal;
        end
        
        function vra = alignFrame(vra, frame_ind, pvec_initial, updateaccumulators, output_stub, verbose)
            % vra = alignFrame(vra, frame_ind, pvec_initial, updateaccumulators, output_stub, verbose)
            % frame_ind - frame index
            % pvec_initial - initial guess; if missing or empty, fits from
            %   scratch with initial orientation guesses and trials
            %   if 'fit' - uses vra.pvec(frame_ind,:) as initial guess
            % updateaccumulators - if true, add fit frame to count and
            %   dwell accumulators
            % output_stub - filename prefix for frames saved to disk
            % verbose - 0: silent;
            %           1: toc every 200 frames
            %           2: toc every 20 frames
            %           3: (default) toc every frame, display image
            % accum_cutoff - if provided, only add to accumulator if
            %   metric(accum_cutoff.field) > accum_cutoff.value
            %   field choices: fval_per_photon, mi2d, mi
            
            if (nargin < 3)
                pvec_initial = [];
            end
             if (nargin < 4 || isempty(updateaccumulators))
                updateaccumulators = false;
             end
            if (nargin < 5)
                output_stub = [];
            end
            if (nargin < 6)
                verbose = 3;
            end
         
            if (length(frame_ind) > 1)
            	vra = vra.setImmodels();
                tic;
                for j = 1:length(frame_ind)
                    try
                        vra = vra.alignFrame(frame_ind(j), pvec_initial, updateaccumulators, output_stub, verbose);
                    catch me
                        disp(me.getReport())
                    end
                    switch verbose
                        case 1
                            if mod(j,200)==0
                                disp([num2str(j) '/' num2str(length(frame_ind))]);
                                toc;
                            end
                        case 2
                            if mod(j,20)==0
                                disp([num2str(j) '/' num2str(length(frame_ind))]);
                                toc;
                            end
                        case 3
                            disp([num2str(j) '/' num2str(length(frame_ind))]);
                            toc;
                    end
                end
                return;
            end

            [ppos, pphase, ptime] = vra.getPhotonLocations(vra.vr.frame.edges(frame_ind + [0 1]));
            
            if isempty(ppos)
                warning('no photons in frame %d!',frame_ind)
                return;
            end
            
            ra = RigidAligner3D(vra.immodel_coarse, ppos, pphase, ptime, vra.z_scale);
            
            ra.xrange = vra.vr.scanBox.corners(1:2,frame_ind) + vra.aligner_range.edgeClearance*[1; -1];
            ra.yrange = vra.vr.scanBox.corners(3:4,frame_ind) + vra.aligner_range.edgeClearance*[1; -1];
            
            ra.czrange = vra.aligner_range.cz;
            ra.rsqrange = vra.aligner_range.rsq;
            
            if (nargin < 3 || isempty(pvec_initial))
                ra = ra.oldfitRigidFromScratch();
            else
                if (strcmpi(pvec_initial, 'fit'))
                    pvec_initial = vra.pvec(frame_ind, :);
                end
                ra = ra.fitRigid(pvec_initial);
            end
            
            
            ra2 = RigidAligner3D(vra.immodel_fine, ppos, pphase, ptime, vra.z_scale);
            ra2.xrange = ra.xrange;
            ra2.yrange = ra.yrange;
            ra2.czrange = ra.czrange;
            ra2 = ra2.fitRigid(ra.fitParameters);
            ra2 = ra2.estimateMutualInformation();
            vra.pvec(frame_ind, :) = ra2.fitParameters;
            vra.metrics(frame_ind).fval_per_photon =  abs(ra2.funPerPhoton);
            vra.metrics(frame_ind).mi = ra2.mifit;
            vra.metrics(frame_ind).mi2D = ra2.mifit2D;
            
            ra3 = ra2;
            ra3.rsqrange = vra.aligner_range.rsq;
            ra3.normByPhotons = true;
            ra3 = ra3.binPhotonsByBinsize();
            vra.metrics(frame_ind).fval_per_photon_in_radius =  abs(ra3.eval(ra3.fitParameters));
            
            
            
            if (updateaccumulators)
                [dpos, dphase, ~, dtau] = vra.getDwellLocationsSubsampled(vra.vr.frame.edges(frame_ind + [0 1]));
                ai.ind = frame_ind;
                ai.metrics = vra.metrics(frame_ind);
                ai.coarse_counts = vra.facoarse.binWithWeightAndAffine(ppos(1,:), ppos(2,:), pphase,[],ra2.fitAffine);
                ai.coarse_dwell = vra.facoarse.binWithWeightAndAffine(dpos(1,:), dpos(2,:), dphase, dtau, ra2.fitAffine);
                ai.fine_counts = vra.fafine.binWithWeightAndAffine(ppos(1,:), ppos(2,:), pphase,[],ra2.fitAffine);
                ai.fine_dwell = vra.fafine.binWithWeightAndAffine(dpos(1,:), dpos(2,:), dphase, dtau, ra2.fitAffine);
                if isempty(vra.accumulated_images)
                    vra.accumulated_images = ai;
                else
                    vra.accumulated_images = [vra.accumulated_images ai];
                end
                if verbose==3
                    try
                        [x,y] = vra.fafine.getBinCenters();
                        subplot(1,2,1);
                        bob = mean(vra.templatefine.F,3);
                        pcolor(vra.templatefine.u, vra.templatefine.v, bob'); axis equal; shading flat; xlim([min(x) max(x)]); ylim([min(y) max(y)]);
                        subplot(1,2,2);
                        pcolor(x,y,(sum(ai.fine_counts,3)./sum(ai.fine_dwell,3))'.*(sum(ai.fine_dwell,3)>6e-9)'); axis equal; shading flat;
                        xlim([min(x) max(x)]); ylim([min(y) max(y)]);
                        caxis([vra.rmin vra.rate_rescale*max(bob(:))]);
                        drawnow
                    catch
                    end
                end
            end
            
        end
        
        function [vra,badframes] = alignFrameNonRigid(vra, frame_ind, options, updateaccumulators)
            existsAndDefault('updateaccumulators',false);
            
            % set up params
            cpg = options.control_point_grid;
            gw = options.geometry_weight;
            al = options.anchor_length;
            
            % set up result variables
            params = cell(size(frame_ind));
            metrics = repmat(vra.metrics_non_rigid(1),length(frame_ind),1);
            if updateaccumulators
                ai = struct('ind',[],'metrics',[],'coarse_counts',[],'coarse_dwell',[],'fine_counts',[],'fine_dwell',[]);
                ai = repmat(ai,[1 length(frame_ind)]);
            else
                ai = [];
            end
            
            % expand aligner range to include full FOV
            arange = vra.aligner_range;
            vra.aligner_range.edgeClearance = 0;
            vra.aligner_range.cz = [-1 1];
            vra.aligner_range.rsq = Inf;
            
            % do registration
            delete(gcp('nocreate')); % release memory
            if ispc
                parpool(2);
            end
            badframes = false(size(frame_ind));
            parfor i=1:length(frame_ind)
                ii = frame_ind(i);
                try
                    % make rigid aligner object
                    rfc = vra.createAligner(ii,vra.immodel_fine,true,true); %#ok<PFBNS>
                    % do non-rigid
                    nra = NonRigidAlignerBSpline(rfc);
                    rfc = []; %#ok<NASGU> % release memory
                    nra.npts = cpg;
                    nra.geometry_weight = gw;
                    nra.anchorlength = al;
                    nra = nra.fitNonRigid;
                    nra = nra.estimateMutualInformation;
                    % update non-rigid results to vra (and temp variables)
                    params{i} = nra.fitParameters;
                    metrics(i).fval_per_photon = nra.funPerPhoton;
                    metrics(i).mi = nra.mifit;
                    metrics(i).mi2D = nra.mifit2D;
                    % update accumulator
                    if updateaccumulators
                        ai(i).ind = ii;
                        ai(i).metrics.fval_per_photon = nra.funPerPhoton;
                        [ct,dw] = nra.binCountsAndDwell(vra.facoarse,true);
                        ai(i).coarse_counts = ct;
                        ai(i).coarse_dwell = dw/vra.rate_rescale;
                        [ct,dw] = nra.binCountsAndDwell(vra.fafine,true);
                        ai(i).fine_counts = ct;
                        ai(i).fine_dwell = dw/vra.rate_rescale;
                        ct = []; dw = []; % release memory
                    end
                    nra = []; % release memory
                catch me
                    fprintf('Error registering frame %d!\n',ii);
                    disp(me.getReport);
                    badframes(i) = true;
                end
            end
            badframes = frame_ind(badframes);
            
            vra.aligner_range = arange;
            vra.pvec_non_rigid(frame_ind) = params;
            vra.metrics_non_rigid(frame_ind) = metrics;
            if updateaccumulators
                vra.resetAccumulators;
                vra.accumulated_images = ai;
            end
        end
        
        function writeAndDisplayVolumes(vra, frame_ind, output_stub, display, dosubframe)
            if (isempty(output_stub) && ~display)
                return;
            end
            
            existsAndDefault('dosubframe', true);
            
            if (isempty(vra.fadisplay))
                fa = vra.fafine;
            else
                fa = vra.fadisplay;
            end
            if (dosubframe)
                if(~vra.hassubframe(frame_ind))
                    vra = vra.fineAlignSubFrame(frame_ind);
                end
                [rigid_red_counts, rigid_red_dwell] = vra.binSubFrame(frame_ind, fa, 'r');
                [rigid_green_counts, rigid_green_dwell] = vra.binSubFrame(frame_ind, fa, 'g');
                time_edges = vra.subframe_time_edges(frame_ind,:);
            else
                [rigid_red_counts, rigid_red_dwell] = vra.binFrame(frame_ind, fa, [], 'r');
                [rigid_green_counts, rigid_green_dwell] = vra.binFrame(frame_ind, fa, [], 'g');
                time_edges = vra.vr.frame.edges(frame_ind + [0 1]);
            end
            if (~isempty(output_stub)) && vra.valid(frame_ind)
                fname = sprintf('%s_frame%d.mat',output_stub, frame_ind);
                [xaxis,yaxis,zaxis] = fa.getBinCenters();
                save(fname, 'rigid_red_counts', 'rigid_red_dwell', 'rigid_green_counts', 'rigid_green_dwell', 'time_edges', 'xaxis', 'yaxis', 'zaxis','-v7.3');
            end
            if (display)
                try
                    figure(3); clf(3);
                    [x,y] = fa.getBinCenters();
                    subplot(2,1,1);
                    bob = mean(vra.templatefine.F,3);
                    pcolor(vra.templatefine.u, vra.templatefine.v, bob'); axis equal; shading flat; xlim([min(x) max(x)]); ylim([min(y) max(y)]);
                    subplot(2,1,2);
                    pcolor(x,y,(sum(rigid_red_counts,[3 4])./sum(rigid_red_dwell,[3 4]))'.*(sum(rigid_red_dwell,[3 4])>6e-9)'); axis equal; shading flat;
                    xlim([min(x) max(x)]); ylim([min(y) max(y)]);
                    caxis([vra.rmin vra.rate_rescale*max(bob(:))]);
                    drawnow
                catch
                end
            end            
        end
        
        function vra = alignSegmentInChunks(vra, inds, chunksize, output_stub, display)
            existsAndDefault('output_stub', []);
            existsAndDefault('display', false);
            n = ceil(length(inds)/chunksize);
            ci = round(linspace(1,length(inds), n+1));
            vra.pvec(:,1) = unwrap(mod(vra.pvec(:,1),2*pi));
            pvi = vra.pvec;
            for j = 1:(length(ci)-1)
                if (j == 1)
                    leftFixedInd = [];
                else
                    leftFixedInd = inds(ci(j))-1;
                end
                vra = vra.alignSegment(inds(ci(j)):inds(ci(j+1)), leftFixedInd, [], output_stub, display);
                if(display)
                    figure(2); clf(2);
                    plot(inds, pvi(inds,1)*180/pi, inds(1:ci(j+1)), vra.pvec(inds(1:ci(j+1)),1)*180/pi); drawnow;
                end
            end
        end

        function ra = createAligner(vra, frame_ind, immodel, isrfc, incldwell, rescaledwell, altcolor)
            %function ra = createAligner(vra, frame_ind, immodel, isrfc, incldwell, rescaledwell, altcolor)
            
            existsAndDefault('rescaledwell',true);
            existsAndDefault('altcolor',vra.beamcolor);
            
            [ppos, pphase, ptime] = vra.getPhotonLocations(vra.vr.frame.edges(frame_ind + [0 1]), altcolor);
            
            if (isrfc)
                ra = RasterFitterChunk3D(immodel, ppos, pphase, ptime, vra.z_scale);
            else
                ra = RigidAligner3D(immodel, ppos, pphase, ptime, vra.z_scale);
            end
            if (incldwell)
                [dpos, dphase, dtime, dwelltime, ~] = vra.getDwellLocationsSubsampled(vra.vr.frame.edges(frame_ind + [0 1]), [], altcolor);
                if rescaledwell
                    dwelltime = dwelltime*vra.rate_rescale;
                end
                ra = ra.setDwellLoc(dpos, dphase, dwelltime, dtime);
            end
            ra.xrange = vra.vr.scanBox.corners(1:2,frame_ind) + vra.aligner_range.edgeClearance*[1; -1];
            ra.yrange = vra.vr.scanBox.corners(3:4,frame_ind) + vra.aligner_range.edgeClearance*[1; -1];            
            ra.czrange = vra.aligner_range.cz;
            ra.rsqrange = vra.aligner_range.rsq;
            
            try
                ra.fitParameters = vra.pvec(frame_ind,:);
                ra.fitAffine = RigidAligner3D.rigidTransform(ra.fitParameters);
                if (isrfc && vra.hassubframe(frame_ind))
                    ra = ra.splitAtTime(vra.subframe_time_edges(frame_ind,2:end-1));
                    for j = 1:length(ra)
                        ra(j).fitParameters = squeeze(vra.pvec_subframe(frame_ind,j,:));
                        ra(j).fitAffine = RigidAligner3D.rigidTransform(ra(j).fitParameters);
                    end
                end
            catch me
                disp(me.getReport())
            end
            
            
        end
            

           

    
         
        
        function vra = alignSegment(vra, inds, leftFixedInd, rightFixedInd,output_stub, display)
           
            existsAndDefault('leftFixedInd', []);
            existsAndDefault('rightFixedInd', []);
            existsAndDefault('output_stub', []);
            existsAndDefault('display', false);
            %assumes initial alignment has already been done
            
           
                        
            %create rfcs for alignment
            for j = 1:length(inds)
                frame_ind = inds(j);
                rfc(j) = vra.createAligner(frame_ind, vra.immodel_coarse, true, false);
%                 [ppos, pphase, ptime] = vra.getPhotonLocations(vra.vr.frame.edges(frame_ind + [0 1]));
%                 rfc(j) = RasterFitterChunk3D(vra.immodel_coarse, ppos, pphase, ptime, vra.z_scale);
%                 
%                 rfc(j).xrange = vra.vr.scanBox.corners(1:2,frame_ind) + vra.aligner_range.edgeClearance*[1; -1];
%                 rfc(j).yrange = vra.vr.scanBox.corners(3:4,frame_ind) + vra.aligner_range.edgeClearance*[1; -1];            
%                 rfc(j).czrange = vra.aligner_range.cz;
%                 rfc(j).rsqrange = vra.aligner_range.rsq;
                rfc(j).k_weight_method = 0;
                
                rfc(j).ub([1 4 5]) = vra.pvec(frame_ind, [1 4 5]) + [deg2rad(15) 3 3];
                rfc(j).lb([1 4 5]) = vra.pvec(frame_ind, [1 4 5]) - [deg2rad(15) 3 3];
                %2,3,6 bounds remain at default
            end
            if(~isempty(leftFixedInd))
                rfc(1).fixedBefore = vra.pvec(leftFixedInd,:);
            end
            if(~isempty(rightFixedInd))
                rfc(end).fixedAfter = vra.pvec(rightFixedInd,:);
            end
            pvec_init = vra.pvec(inds, :);
            
            tic
            rfc = rfc.fit(pvec_init);
            toc
            pvec_intermediate = cat(1,rfc.fitParameters); 
            
            rfc = rfc.createImModel(vra.immodel_fine);
            [rfc, ~, fval, ~, ~, fval_image_only] = rfc.fit(pvec_intermediate);
            toc
            vra.pvec(inds,:) = cat(1,rfc.fitParameters); 
            vra.pvec(:,1) = unwrap(mod(vra.pvec(:,1),2*pi));
            
            if (display)
                figure(1); clf(1)
                subplot(2,2,1); plot(inds, pvec_init(:,1)*180/pi, inds, pvec_intermediate(:,1)*180/pi, inds,  vra.pvec(inds,1)*180/pi); title('theta_z'); legend('init', 'coarse', 'fine')
                subplot(2,2,2); plot(inds, pvec_intermediate(:,2:3)*180/pi, inds,  vra.pvec(inds,2:3)*180/pi); title('other angles'); legend('coarsex', 'coarsey', 'finex','finey')
                subplot(2,2,3); plot(inds, pvec_init(:,4:5), inds, pvec_intermediate(:,4:5), inds,  vra.pvec(inds,4:5)); title('dx dy'); legend('initx','inity', 'coarsex', 'coarsey', 'finex','finey')
                subplot(2,2,4); plot(inds, pvec_intermediate(:,6), inds,  vra.pvec(inds,6)); title('dz');legend('coarsez', 'finez')
            end
            try
                for j = 1:length(inds)
                    frame_ind = inds(j);
                    vra.metrics(frame_ind).fval_sequence = abs(fval);
                    vra.metrics(frame_ind).fval_per_photon =  abs(fval_image_only(j));
                    rfc(j) = rfc(j).estimateMutualInformation();
                    vra.metrics(frame_ind).mi = rfc(j).mifit;
                    vra.metrics(frame_ind).mi2D = rfc(j).mifit2D;
                    vra.writeAndDisplayVolumes(inds(j), output_stub, display && j == 1);
                end
            catch me
                warning('failed to calculate metrics or display volumes');
                disp(me.getReport());
            end
                
            
        end
        
        function [rigid_red_counts, rigid_red_dwell, rigid_green_counts, rigid_green_dwell, xaxis, yaxis, zaxis, time_edges] = loadSavedResults(vra, options, frame_ind, blur_sigma, xrange, yrange, zrange)
            
            existsAndDefault('blur_sigma', []);
            
            load(sprintf('%s_frame%d.mat',options.output_stub, frame_ind), 'xaxis', 'yaxis', 'zaxis', 'rigid_red_counts','rigid_red_dwell', 'rigid_green_counts', 'rigid_green_dwell', 'time_edges');
            function v = opval(fieldname, default)
                try 
                    if (~isstruct(options))
                        options = struct([]);
                    end
                catch
                    options = struct([]);
                end
                try
                    v = options.(fieldname);
                catch
                    v = default;
                end
                if (isempty(v))
                    v = default;
                end
                options(1).(fieldname) = v;
            end
            existsAndDefault('xrange',opval('final_template_xrange', []));
            existsAndDefault('yrange', opval('final_template_yrange', []));
            existsAndDefault('zrange', opval('final_template_zrange', []));
            
            
            if (~isempty(xrange))
                xi = xaxis >= min(xrange) & xaxis <= max(xrange);
            else
                xi = true(size(xaxis));
            end
            
            if (~isempty(yrange))
                yi = yaxis >= min(yrange) & yaxis <= max(yrange);
            else
                yi = true(size(yaxis));
            end
            
            if (~isempty(zrange))
                zi = zaxis >= min(zrange) & zaxis <= max(zrange);
            else
                zi = true(size(zaxis));
            end
            
            xaxis = xaxis(xi);
            yaxis = yaxis(yi);
            zaxis = zaxis(zi);
            
            if (~isempty(blur_sigma))
                sigma = blur_sigma./([median(diff(xaxis)) median(diff(yaxis)) median(diff(zaxis))]);
                sigma(4) = 0;
                rigid_red_counts = imblur(rigid_red_counts, sigma);
                rigid_red_dwell = imblur(rigid_red_dwell, sigma);
                rigid_green_counts = imblur(rigid_green_counts, sigma);
                rigid_green_dwell = imblur(rigid_green_dwell, sigma);
                
            end
            
            rigid_red_counts = rigid_red_counts(xi,yi,zi,:);
            rigid_red_dwell = rigid_red_dwell(xi,yi,zi,:);
            rigid_green_counts = rigid_green_counts(xi,yi,zi,:);
            rigid_green_dwell = rigid_green_dwell(xi,yi,zi,:);
            
            
        end
        
        function template = makeTemplateFromSavedResults(vra, options)
            function v = opval(fieldname, default)
                try 
                    if (~isstruct(options))
                        options = struct([]);
                    end
                catch
                    options = struct([]);
                end
                try
                    v = options.(fieldname);
                catch
                    v = default;
                end
                if (isempty(v))
                    v = default;
                end
                options(1).(fieldname) = v;
            end
            
            template_selection_field = opval('template_selection_field', 'fval_per_photon');
            
            final_template_num_frames = opval('final_template_num_frames', 20);
            
            [tsf,I] = sort([vra.metrics.(template_selection_field)], 'descend');
            I = I(tsf > 0);
            [red_counts, ~, ~, ~, xaxis, yaxis, zaxis] = loadSavedResults(vra, options, I(1));
%             
%             fname = @(frame_ind) sprintf('%s_frame%d.mat',options.output_stub, frame_ind);
%             load(fname(I(1)), 'xaxis', 'yaxis', 'zaxis', 'red_counts');
            template.counts = 0*red_counts;
            template.dwell = template.counts;
            
            template.rawcounts = template.counts;
            template.rawdwell = template.dwell;
            
            
            blur_sigma = opval('final_template_blur_sigma', []);
            for j = 1:min(final_template_num_frames, length(I))
                
               [red_counts, red_dwell, green_counts, green_dwell] = loadSavedResults(vra, options, I(j), blur_sigma);
                if (lower(vra.beamcolor) == 'r' || lower(vra.beamcolor) == 'a')
%                      load(fname(I(j)), 'red_counts', 'red_dwell');
                      template.counts = template.counts + red_counts;
                     template.dwell = template.dwell + red_dwell;
                end 
                if (lower(vra.beamcolor) == 'g' || lower(vra.beamcolor) == 'a')
%                      load(fname(I(j)), 'green_counts', 'green_dwell');
                      template.counts = template.counts + green_counts;
                     template.dwell = template.dwell + green_dwell;
                end
            end
            
            template.counts = sum(template.counts,4);
            template.dwell = sum(template.dwell, 4);
            template.u = xaxis;
            template.v = yaxis;
            template.w = zaxis;                      
        end
       
        function outputTemplateAndImagesForAlignment(vra, inds, options)
            template = vra.makeTemplateFromSavedResults(options);
            save([options.final_output_stub '_template.mat'], 'template');
            try 
                blur_sigma = options.final_template_blur_sigma;
            catch
                blur_sigma = [];
            end
            for j = 1:length(inds)
                if (~isempty(blur_sigma))
                    [red_counts, red_dwell, green_counts, green_dwell, xaxis, yaxis, zaxis, time_edges] = loadSavedResults(vra, options, inds(j), blur_sigma);
                    save(sprintf('%s_blur%d.mat',options.final_output_stub,inds(j)), 'red_counts', 'red_dwell', 'green_counts', 'green_dwell', 'xaxis', 'yaxis', 'zaxis', 'time_edges')
                end
                [red_counts, red_dwell, green_counts, green_dwell, xaxis, yaxis, zaxis, time_edges] = loadSavedResults(vra, options, inds(j));
                save(sprintf('%s_raw%d.mat',options.final_output_stub,inds(j)), 'red_counts', 'red_dwell', 'green_counts', 'green_dwell', 'xaxis', 'yaxis', 'zaxis', 'time_edges')
            end
            
            
            
            
        end
        
        function metric = saveFrameToDisk(vra,frame_ind,fstub,fa,rigid_pvec,non_rigid_pvec)
            existsAndDefault('rigid_pvec',[]); % don't do rigid if empty (ATTENTION: this behavior differs from elsewhere in VRAligner)
            existsAndDefault('non_rigid_pvec',[]); % don't do non-rigid if empty
            
            % use NonRigidAlignerBSpline.binCountsAndDwell() for everything
                        
            fpath = [fstub sprintf('_frame%d.mat',frame_ind)];
            a = matfile(fpath,'Writable',true);
            a.time_edges = vra.vr.frame.edges(frame_ind+[0 1]);
            [xc,yc,zc] = fa.getBinCenters;
            a.xaxis = xc;
            a.yaxis = yc;
            a.zaxis = zc;
            a.with_rigid = ~isempty(rigid_pvec);
            a.with_non_rigid = ~isempty(non_rigid_pvec);
            a.rigid_pvec = rigid_pvec;
            a.non_rigid_pvec = non_rigid_pvec;
            
            % get red volumes
            r_rfc = vra.createAligner(frame_ind,vra.immodel_fine,true,true,false,'r');
            if isempty(rigid_pvec) % if input is empty, reset rigid to all zero
                for i=1:length(r_rfc)
                    r_rfc(i).fitParameters = [0 0 0 0 0 0];
                    r_rfc(i).fitAffine = RigidAligner3D.rigidTransform([0 0 0 0 0 0]);
                end
            else
                for i=1:length(r_rfc)
                    r_rfc(i).fitParameters = rigid_pvec;
                    r_rfc(i).fitAffine = RigidAligner3D.rigidTransform(rigid_pvec);
                end
            end
            r_nra = NonRigidAlignerBSpline(r_rfc);
            if ~isempty(non_rigid_pvec)
                r_nra = r_nra.setFitParameters(non_rigid_pvec);
            end
            [rc,rd] = r_nra.binCountsAndDwell(fa,~isempty(non_rigid_pvec));
            
            % get green volumes
            g_rfc = vra.createAligner(frame_ind,vra.immodel_fine,true,true,false,'g');
            if isempty(rigid_pvec) % if input is empty, reset rigid to all zero
                for i=1:length(g_rfc)
                    g_rfc(i).fitParameters = [0 0 0 0 0 0];
                    g_rfc(i).fitAffine = RigidAligner3D.rigidTransform([0 0 0 0 0 0]);
                end
            else
                for i=1:length(g_rfc)
                    g_rfc(i).fitParameters = rigid_pvec;
                    g_rfc(i).fitAffine = RigidAligner3D.rigidTransform(rigid_pvec);
                end
            end
            g_nra = NonRigidAlignerBSpline(g_rfc);
            if ~isempty(non_rigid_pvec)
                g_nra = g_nra.setFitParameters(non_rigid_pvec);
            end
            [gc,gd] = g_nra.binCountsAndDwell(fa,~isempty(non_rigid_pvec));
            
            % save to disk
            if isempty(rigid_pvec) && isempty(non_rigid_pvec) % raw
                a.raw_red_counts = rc;
                a.raw_red_dwell = rd;
                a.raw_green_counts = gc;
                a.raw_green_dwell = gd;
            elseif ~isempty(rigid_pvec) && isempty(non_rigid_pvec) % rigid
                a.rigid_red_counts = rc;
                a.rigid_red_dwell = rd;
                a.rigid_green_counts = gc;
                a.rigid_green_dwell = gd;
            elseif ~isempty(rigid_pvec) && ~isempty(non_rigid_pvec) % non-rigid
                a.non_rigid_red_counts = rc;
                a.non_rigid_red_dwell = rd;
                a.non_rigid_green_counts = gc;
                a.non_rigid_green_dwell = gd;
            else % unexpected combination, save without prefix and together with pvecs for future reference
                a.red_counts = rc;
                a.red_dwell = rd;
                a.green_counts = gc;
                a.green_dwell = gd;
            end
            
            % calculate metric using red rate volume against template
            if nargout>0
                % get registered volume
                rv = rc./rd;
                rv(~isfinite(rv)) = 0;
                rv = imgaussfilt3(rv,1);
                % get template volume
                [XC,YC,ZC] = ndgrid(xc,yc,zc);
                v0 = vra.immodel_fine.F(XC,YC,ZC);
                v0 = v0*vra.rate_rescale;
                % get dwell mask for this frame
                bw = imclose(dw>0,ones([9 9 9]));
                % calculate multissim3
                metric = multissim3(rv.*bw,v0.*bw);
            end
        end
        
    end
    
    methods (Static)
        
        function op = processOptions() 
            function setval(fieldname, default)
               op.(fieldname) = default;
            end
            setval('scanbeam', 2); %which beam forms the images (usually beam 2 for hyperscope)
            setval('trange_galvoDT', []); %which time range to use to find galvo delta T - [] means automatic using parameters below
            setval('tw_lowpass', 0.1); %lowpass time window to determine when larva is still for galvoDT correction
            setval('tw_galvoDT', 5); %length of time window over which to determine galvoDT
            setval('template_timerange',[]); %time range to use to create initial template, [] means automatic determination
            setval('orient_template', true); %orient template automatically to align with axes
            setval('trange_alignment', []); %range of time to align, [] = default time tracker is runing
            setval('template_num_align', 100); %how many frames to align (spaced over entire time range) to determine final template
            setval('template_percentile', 0.75); %percetile cut for alignment agreement to make template (default of .75 means top 25 out of 100 frames used)
            setval('dstdir',[]); %where to save registration results
            setval('template_rmin', []); %template_rmin - [] means determine automagically
            setval('final_template_blur_sigma', []); %if template is to be blurred
            setval('save_frames',true); %whether to save registered frames to disk
        end
        function vra = processDirectory(srcdir, options)
            % vra = processDirectory(srcdir, options)
            % srcdir contains raw data (.bin file and associated text files
            % for tracker etc.)
            % options are optional -- see vra.processOptions for description
            % and default
            %
            % in addition to returning vra, saves following to
            % srcdir/registered
            % timestamp_vr.mat -- volume reconstruction with galvoDT and
            % phase corrections
            % timestamp_template.mat -- initial and final templates for
            % registration
            function v = opval(fieldname, default)
                try 
                    if (~isstruct(options))
                        options = struct([]);
                    end
                catch
                    options = struct([]);
                end
                try
                    v = options.(fieldname);
                catch
                    v = default;
                end
                if (isempty(v))
                    v = default;
                end
                options(1).(fieldname) = v;
            end
            scanbeam = opval('scanbeam', 2);
            savevr = false;
            timestamp = split(srcdir,filesep); %#ok<*UNRCH>
            timestamp = timestamp(~cellfun(@isempty, timestamp, 'uniformoutput', true));
            timestamp = timestamp{end};
            options.timestamp = timestamp;
            dstdir = opval('dstdir',fullfile(srcdir,'registered'));
            if (~exist(dstdir, 'dir'))
                mkdir(dstdir);
            end
            
            vrname =  [timestamp '_vr.mat'];
            try
                load(fullfile(srcdir, vrname), 'vr'); % see if vr.mat is saved during preview
                vr = vr.changeFilename(srcdir);
                vr.vrfa.vrfc.unload();
                disp('vr loaded from mat file');
            catch
                vr = VolumeReconstructor(srcdir);
%                 is = vr.settings(end).imagingSettings;
%                 tsr = vr.tsr;
                vr = vr.addBehaviorVideo;
                vr = vr.getStage;
                vr = vr.applyTrackerOffset;
                vr = vr.pongFrames();    
                savevr = true;
                disp('vr created from bin file');
                % TEMP FIX: since full vr saving doesn't work, load fields
                % calculated and saved during preview to save time - Rui
                %PLEASE DON'T DO THIS - If there is a bug, fix it
                %If you can't fix the bug, do this kind of kludge on your
                %own branch
            end
            if (~all(isfield(vr.phaseCorrection, {'t', 'phi'})))
                vr = vr.calculatePhaseCorrection();
                savevr = true;
                disp('phaseCorrection calculated');
            end
            if (isempty(vr.galvoDeltaT))
                trange_galvoDT = opval('trange_galvoDT', []);
                tw_lowpass = opval('tw_lowpass', 0.1);
                tw_galvoDT = opval('tw_galvoDT', 5);
                if isempty(trange_galvoDT)
                    trange_galvoDT = vr.getStillTimes2(tw_lowpass,tw_galvoDT);
                end
                vr = vr.calculateGalvoDeltaT(scanbeam,'a',trange_galvoDT);
                savevr = true;
                disp('galvoDeltaT calculated');
            end
            if (isempty(vr.scanBox) || (size(vr.scanBox.corners,2) < (length(vr.frame.edges) -1)))
                vr = vr.calculateScanBox;
                savevr = true;
                disp('scanBox calculated');
            end
           
            if (savevr)
                save(fullfile(srcdir, vrname), 'vr', '-v7.3');
                disp('vr saved');
            end
            
            vra = VRAligner(vr);
            templatename = [timestamp '_template.mat'];
            
            trange_alignment = opval('trange_alignment', vra.trackingTimeRange);
            
            try
                load(fullfile(dstdir, templatename), 'init_template');
                vra.templatefine = init_template.fine;
                vra.templatecoarse = init_template.coarse;                
                disp('initial template loaded from disk');
            catch
                template_timerange = opval('template_timerange',[]);
                if (isempty(template_timerange))
                    template_timerange = vr.getStillTimes2(0.1,1,false,trange_alignment);
                end
                orient_template = opval('orient_template', true);
                vra.rmin = opval('template_rmin', []);
                vra = vra.initializeTemplates(template_timerange, orient_template);
                init_template.fine = vra.templatefine;
                init_template.coarse = vra.templatecoarse;
                
                save(fullfile(dstdir, templatename), 'init_template');
                disp('initial template created');
            end
            
            figure(11); clf();
            subplot(2,2,1); pcolor(vra.templatefine.u, vra.templatefine.v, vra.templatefine.F_2D'); shading flat; axis equal; title ('initial fine template'); drawnow
            subplot(2,2,2); pcolor(vra.templatefine.u, vra.templatefine.w, squeeze(max(vra.templatefine.F,[],2))'); shading flat; axis equal; title ('initial fine template'); drawnow
            subplot(2,2,3); pcolor(vra.templatefine.w, vra.templatefine.v, squeeze(max(vra.templatefine.F,[],1))); shading flat; axis equal; title ('initial fine template'); drawnow
                        
            try
                load(fullfile(dstdir, templatename), 'template');
                vra.templatefine = template.fine;
                vra.templatecoarse = template.coarse;  
                disp('full template loaded from disk');
            catch
                template_num_align = opval('template_num_align', 100);
                template_percentile = opval('template_percentile', 0.75);
                inds = find(vra.vr.frame.edges > min(trange_alignment) & vra.vr.frame.edges < max(trange_alignment));
                di = max(1,length(inds)/template_num_align);
                inds = inds(round(1:di:length(inds)));
                vra = vra.resetAccumulators();
                vra = vra.alignFrame(inds, [], true, [], 2);
                ac.field = 'fval_per_photon_in_radius';
                ac.percentile = template_percentile;
                vra = vra.updateTemplatesFromAccumulators(ac);
                template.fine = vra.templatefine;
                template.coarse = vra.templatecoarse;
                
                save(fullfile(dstdir, templatename), 'template', '-append');
                disp('full template created');
            end
            inds = find(vra.vr.frame.edges > min(trange_alignment) & vra.vr.frame.edges < max(trange_alignment));
            
             figure(12); clf();
            subplot(2,2,1); pcolor(vra.templatefine.u, vra.templatefine.v, vra.templatefine.F_2D'); shading flat; axis equal; title ('full fine template'); drawnow
            subplot(2,2,2); pcolor(vra.templatefine.u, vra.templatefine.w, squeeze(max(vra.templatefine.F,[],2))'); shading flat; axis equal; title ('full fine template'); drawnow
            subplot(2,2,3); pcolor(vra.templatefine.w, vra.templatefine.v, squeeze(max(vra.templatefine.F,[],1))); shading flat; axis equal; title ('full fine template'); drawnow
           
            try 
                load (fullfile(dstdir, [timestamp '_initial2Dalignment.mat']), 'pvec_initial');
                vra.pvec = pvec_initial;
                vra = vra.setImmodels();

            catch
                vra = vra.initial2DAlignment(inds);                
                pvec_initial = vra.pvec;
                save(fullfile(dstdir, [timestamp '_initial2Dalignment.mat']), 'pvec_initial');
            end

            vra = vra.createDisplayAssembler(); 
            vra = vra.alignSegmentInChunks(inds,20,[], false);

            vra = vra.fineAlignSubFrame(inds);
            
            options.alignedInds = intersect(inds,find(vra.valid));
            
            ai = vra.accumulated_images;
            vr = vra.vr;
            vra.accumulated_images = [];
            vra.vr = [];
            save(fullfile(dstdir, [timestamp '_vra.mat']), 'vra','-v7.3');
            vraopts = options;
            save(fullfile(dstdir, [timestamp '_vra_options.mat']), 'vraopts');
            
            if ~options.save_frames
                disp('vra saved');
                return;
            end
            % stop here if not saving registered frames
            % (do all frame saving and zipping on hpc to get around google
            % drive quotas)
            
            disp('vra saved - now writing individual frames to disk');
            
            if (~exist(fullfile(dstdir, 'registered frames'), 'dir'))
                mkdir(fullfile(dstdir, 'registered frames'));
            end
            output_stub = fullfile(dstdir, 'registered frames', timestamp);
            options.output_stub = output_stub;
            
            vra.accumulated_images = ai;
            vra.vr = vr;
            
            for j = 1:length(inds)
                vra.writeAndDisplayVolumes(inds(j), output_stub, mod(j,20) == 0, true);
            end
            
            disp('all frames aligned');
            
            mainimage = imopen(vra.templatefine.F > percentile(vra.templatefine.F(vra.templatefine.F>vra.rmin/vra.rate_rescale),0.95),ones([7 7 7]));
            mainimage = imdilate(mainimage,ones([5 5 5]));
            options.final_template_xrange = opval('final_template_xrange', [min(vra.templatefine.u(any(mainimage,[2 3]))) max(vra.templatefine.u(any(mainimage,[2 3])))]);
            options.final_template_yrange = opval('final_template_yrange', [min(vra.templatefine.v(any(mainimage,[1 3]))) max(vra.templatefine.v(any(mainimage,[1 3])))]);
            options.final_template_zrange = opval('final_template_zrange', [min(vra.templatefine.w(any(mainimage,[1 2]))) max(vra.templatefine.w(any(mainimage,[1 2])))]);
            options.final_output_stub = fullfile(dstdir, 'ready for nonrigid', timestamp);
            if (~exist(fullfile(dstdir, 'ready for nonrigid'), 'dir'))
                mkdir(fullfile(dstdir, 'ready for nonrigid'));
            end
            
            vra.outputTemplateAndImagesForAlignment(inds, options);
            
%             ai = vra.accumulated_images;
%             vr = vra.vr;
%             vra.accumulated_images = [];
%             vra.vr = [];
%             save(fullfile(dstdir, [timestamp '_vra.mat']), 'vra','-v7.3');
%             disp('vra saved - done!');
%             vra.accumulated_images = ai;
%             vra.vr = vr;
            vraopts = options;
            save(fullfile(dstdir, [timestamp '_vra_options.mat']), 'vraopts');
           
           end
                
    end
    
   
    
    
end

