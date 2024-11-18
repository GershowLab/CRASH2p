classdef HMMFrameBuffer_VRA < HMMFrameBuffer
    % HMMFrameBuffer_VRA
    %   volumetric frame buffer for HyperMovieMaker_VR
    %   display options:
    %   - voxelsize: default = [1 1 2] um
    %   - blursigma: default = 1.3 for x/y/z/t
    %   - mindwell: default = FPGA clock tick time, shouldn't be smaller
    %               than this
    %   - fixFOV: default = false; for datasets with scan offset changing
    %             during recording, keep the neurons fixed and move the FOV
    
    properties
        vra
        do_ic = false;
        ic;
    end
    
    methods
        
        function fb = HMMFrameBuffer_VRA(vra,op)
            if (~isfield(op, 'fa'))
                op.fa = 'initlater';
                dofa = true;
            else
                dofa = false;
            end
            fb@HMMFrameBuffer(vra.vr, op);
            fb.vra = vra;
            if(dofa)
                fb.set_default_assembler;
            end
            fb.thisFrame.tx = [];
            fb.nextFrame.tx = [];
            fb.prevFrame.tx = [];


        end

        function im = get_projection(fb,t,clr,prj_ind)
            fb.set_time(t);
            [~,ind] = min(abs(t-fb.thisFrame.tx));
            ct = fb.thisFrameSmoothed.([clr(1) 'c']);
            dw = fb.thisFrameSmoothed.([clr(1) 'd']);
            ct = squeeze(sum(ct(:,:,:,ind),4-prj_ind));
            dw = squeeze(sum(dw(:,:,:,ind),4-prj_ind));
            alpha = fb.thisFrameSmoothed.alpha{prj_ind};
            im = ct./dw./alpha;
            if prj_ind==1 || prj_ind==2
                im = im';
            end

        end
        
        function set_time(fb,t)
            frameInd = find(fb.vra.valid' & (fb.vra.subframe_time_edges(:,1) < t), 1,'last');
            if (isempty(frameInd))
                frameInd = find(fb.vra.valid,1,'first');
            end
            fb.load_frame(frameInd);
           

        end
        
        function frame = load_frame_helper(fb,frameInd)
            frameInd = max(find(fb.vra.valid,1,'first'), min(frameInd, find(fb.vra.valid,1,'last')));
            [rc,rd] = fb.vra.binSubFrame(frameInd, fb.fa, false);
            [gc,gd] = fb.vra.binSubFrame(frameInd, fb.fa, true);

            frame.rc = rc;
            frame.rd = rd;
            frame.gc = gc;
            frame.gd = gd;
            frame.tedges = fb.vra.subframe_time_edges(frameInd,:);
            frame.tx = conv(frame.tedges, [.5 .5], 'valid');
        end

        function load_prev_frame(fb,frameInd)
            % loading frameInd-1 volumes into fb.prevFrame
            
            fb.prevFrame = fb.load_frame_helper(frameInd-1);
        end
        
        function load_this_frame(fb,frameInd)

            fb.thisFrame = fb.load_frame_helper(frameInd);
        end
        
        function load_next_frame(fb,frameInd)
           fb.nextFrame = fb.load_frame_helper(frameInd+1);
        end
                
        function get_smoothed_volume(fb)
            rc = imblur(cat(4,fb.prevFrame.rc,fb.thisFrame.rc,fb.nextFrame.rc),fb.blursigma);
            rd = imblur(cat(4,fb.prevFrame.rd,fb.thisFrame.rd,fb.nextFrame.rd),fb.blursigma);
            gc = imblur(cat(4,fb.prevFrame.gc,fb.thisFrame.gc,fb.nextFrame.gc),fb.blursigma);
            gd = imblur(cat(4,fb.prevFrame.gd,fb.thisFrame.gd,fb.nextFrame.gd),fb.blursigma);
            valid = isfinite(rd) & isfinite(gd) & rd > fb.mindwell & gd > fb.mindwell;
            rc(~valid) = 0;
            rd(~valid) = 0;
            gc(~valid) = 0;
            gd(~valid) = 0;
            inds = fb.vra.nsubframe+(1:fb.vra.nsubframe);
            fb.thisFrameSmoothed.rc = rc(:,:,:,inds);
            fb.thisFrameSmoothed.rd = rd(:,:,:,inds);
            fb.thisFrameSmoothed.gc = gc(:,:,:,inds);
            fb.thisFrameSmoothed.gd = gd(:,:,:,inds);
            for prj_ind = 1:3

                if (fb.do_ic)
                    ct = squeeze(sum(rc,[4 4-prj_ind]));
                    dw = squeeze(sum(rd, [4 4-prj_ind]));
                    fb.thisFrameSmoothed.alpha{prj_ind} = fb.ic(prj_ind).calculate_correction(ct,dw);
                else
                    fb.thisFrameSmoothed.alpha{prj_ind} = 1;
                end
            end

        end

        function set_intensity_correction(fb, do_ic)
            fb.do_ic = do_ic;
            if (do_ic)
                [xc,yc,zc] = fb.fa.getBinCenters();
                xi = interp1(fb.vra.templatefine.u, 1:length(fb.vra.templatefine.u), xc, 'nearest');
                yi = interp1(fb.vra.templatefine.v, 1:length(fb.vra.templatefine.v), yc, 'nearest');
                zi = interp1(fb.vra.templatefine.w, 1:length(fb.vra.templatefine.w), zc, 'nearest');
                
                nx = max(3,ceil(diff(xc([1 end]))/20));
                ny = max(3,ceil(diff(yc([1 end]))/20));
                nz = max(3,ceil(diff(zc([1 end]))/20));

                F = fb.vra.templatefine.F*fb.vra.rate_rescale;

            
                fb.ic = [IntensityCorrectorBSpline(squeeze(mean(F(xi,yi,:),3)), [nx ny]);
                    IntensityCorrectorBSpline(squeeze(mean(F(xi,:,zi),2)), [nx nz]);
                    IntensityCorrectorBSpline(squeeze(mean(F(:,yi,zi),1)), [ny nz])];      
                for j = 1:3
                    fb.ic(j).dwell_min = 2.5e-7; %longer minimum dwell because of projection
                end
            end
        end
        

        function fa = get_default_assembler(fb)

            temp = sort(fb.vra.templatefine.F_2D(:),'descend');
            thresh = temp(1000); %1000th brightest voxel
            
            bwim = fb.vra.templatefine.F_2D > thresh;
            bwim = imdilate(imopen(bwim, ones(3)),ones(9));



            irng = find(any(bwim,[2 3]),1,'first'):find(any(bwim,[2 3]),1,'last');
            jrng = find(any(bwim,[1 3]),1,'first'):find(any(bwim,[1 3]),1,'last');


            val = squeeze(max(fb.vra.templatefine.F(irng, jrng,:),[],[1 2]));
            [~,I] = max(val);
            
            krng = find(val(1:I)<median(val),1,'last'):(I-1+find(val(I:end)<median(val),1,'first'));


            xc = fb.vra.templatefine.u(irng);
            yc = fb.vra.templatefine.v(jrng);
            zc = fb.vra.templatefine.w(krng);

            fa = FastAssembler(xc,yc,zc,[],true,fb.vr.z_scale);
        end

    end
    
    
end