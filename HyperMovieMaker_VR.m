classdef HyperMovieMaker_VR < HyperMovieMaker
    %HyperMovieMaker_VR makes preview movie for hyperscope datasets
    %   Syntax: hmm = HyperMovieMaker_VR(vr,tsr,options);
    %
    %   VR-specific options:
    %   - fa: default range determined by imaging settings
    %   - vrframerate: default determined by vr.pongFrames()
    %   - voxelsize: default = [1 1 2] um, only used to create default
    %                assembler; should rerun hmm = hmm.set_default_assembler
    %                if modified after object creation
    %   - blursigma: default = 1.3 for x/y/z/t
    %   - mindwell: default = FPGA clock tick time, shouldn't be smaller
    %               than this
    %   - fixFOV: default = false; for datasets with scan offset changing
    %             during recording, keep the neurons fixed and move the FOV
    %
    %   To change assembler after object creation:
    %       hmm = hmm.set_assembler(new_fa);
    %   To change VR frame rate (not the same as movie frame rate) after
    %   object creation:
    %       hmm = hmm.set_vr_framerate(new_framerate);
    %   To change all other options listed above after object creation,
    %   modify the corresponding field in hmm.frameBuffer directly (changes
    %   won't be reflected in hmm.op fields though)

    properties
        frameBuffer
    end

    methods
        function obj = HyperMovieMaker_VR(vr,tsr,opp)
            if (nargin < 3)
                opp = [];
            end
            obj@HyperMovieMaker(tsr,opp);
            % for vr backward compatibility
            if ~isfield(vr.scanBox,'x')
                vr = vr.applyScanOffset;
            end
            vr.z_scale = vr.settings(end).analogScalingFactors(2).zTagMicronsPerTick*160E6/(2*pi*vr.settings(end).analogScalingFactors(2).tagLensFrequency);
            if isfield(obj.op,'vrframerate')
                vr = vr.pongFrames(1/obj.op.vrframerate);
            end
            obj.frameBuffer = HMMFrameBuffer(vr,opp);
        end
        
        function [im, uaxis, vaxis] = get_projection(obj, t, clr, prj_ind)
            % t < time (according to FPGA clock) at which to get projections
            % clr < color, either 'r'/'red', or 'g'/'green'
            % prj_ind < 1 = xy, 2 = xz, 3 = zy
                        
            % get u/vaxis first
            [xc,yc,zc] = obj.frameBuffer.fa.getBinCenters;
            switch prj_ind
                case 1
                    uaxis = xc; vaxis = yc;
                case 2
                    uaxis = xc; vaxis = zc;
                case 3
                    uaxis = zc; vaxis = yc;
            end
            im = NaN(length(vaxis),length(uaxis),length(t));
            
            % get binned and 4D smoothed volume
            frametimes = obj.frameBuffer.vr.frame.edges(1:end-1)+obj.frameBuffer.vr.frame.frametime/2;
            for j=1:length(t)
                [~,frameInd] = min(abs(frametimes-t(j)));
                obj.frameBuffer.load_frame(frameInd);
                ct = obj.frameBuffer.thisFrameSmoothed.([clr(1) 'c']);
                dw = obj.frameBuffer.thisFrameSmoothed.([clr(1) 'd']);
                thisim = squeeze(sum(ct,4-prj_ind)./sum(dw,4-prj_ind));
                % match display xy order
                if prj_ind==1 || prj_ind==2
                    thisim = thisim';
                end
                im(:,:,j) = thisim;
            end
        end
        
        function tx = get_default_time_axis(obj)
            tx = obj.frameBuffer.vr.frame.edges(1:end-1) + obj.frameBuffer.vr.frame.frametime/2;
        end
        
        function obj = get_default_proj_inds(obj)
            if isempty(obj.xi)
                obj.xi = 1:obj.frameBuffer.fa.nx;
            end
            if (isempty(obj.yi))
                obj.yi = 1:obj.frameBuffer.fa.ny;
            end
            if (isempty(obj.zi))
                obj.zi = 1:obj.frameBuffer.fa.nz;
            end
        end
        
        function fname = get_default_fname(obj)
            fname = strrep(obj.frameBuffer.vr.filename,'photontiming.bin','vr_movie.mp4');
        end
        
        %-----------------------
        % vr-specific functions
        %-----------------------
        
        function obj = set_assembler(obj,fa)
            obj.frameBuffer.fa = fa;
            obj.op.fa = fa;
        end
        
        function obj = set_default_assembler(obj)
            obj.frameBuffer.set_default_assembler;
            obj.op.fa = obj.frameBuffer.fa;
        end
        
        function fa = get_default_assembler(obj)
            fa = obj.frameBuffer.get_default_assembler;
        end
        
        function obj = set_vr_framerate(obj,framerate)
            obj.frameBuffer.vr = obj.frameBuffer.vr.pongFrames(1/framerate);
            obj.op.vrframerate = framerate;
        end
        
        function obj = set_projection_box(obj)
            obj.proj_overlay.do_overlay = false;
        end
        
                
    end
end