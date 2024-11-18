classdef HyperMovieMaker_VRA < HyperMovieMaker_VR
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
        
    end

    methods
        function obj = HyperMovieMaker_VRA(vra,tsr,opp)
            if (nargin < 3)
                opp = [];
            end
            obj@HyperMovieMaker_VR(vra.vr, tsr,opp);
            obj.frameBuffer = HMMFrameBuffer_VRA(vra,opp); %check with Rui about possible bug in _VR
        end

        function obj = set_intensity_correction(obj, do_ic)
            obj.frameBuffer.set_intensity_correction(do_ic);
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
            
            
            for j=1:length(t)
                im(:,:,j) = obj.frameBuffer.get_projection(t(j), clr, prj_ind);
            end
        end
        
        function tx = get_default_time_axis(obj)
            tx = reshape(conv2(obj.frameBuffer.vra.subframe_time_edges(obj.frameBuffer.vra.valid,:), [.5 .5], 'valid')',[],1);
        end
         
        function fname = get_default_fname(obj)
            fname = strrep(obj.frameBuffer.vr.filename,'photontiming.bin','vra_movie.mp4');
        end
        
        %-----------------------
        % vr-specific functions
        %-----------------------
        

        
        function obj = set_vr_framerate(obj,framerate)
            obj.frameBuffer.vr = []; %try to avoid making memory copies
            obj.frameBuffer.vra.vr = obj.frameBuffer.vra.vr.pongFrames(1/framerate);
            obj.frameBuffer.vr = obj.frameBuffer.vra.vr;
            obj.op.vrframerate = framerate;
        end
        
        
                
    end
end