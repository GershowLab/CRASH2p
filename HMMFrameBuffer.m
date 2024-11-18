classdef HMMFrameBuffer < handle
    % HMMFrameBuffer
    %   volumetric frame buffer for HyperMovieMaker_VR
    %   display options:
    %   - voxelsize: default = [1 1 2] um
    %   - blursigma: default = 1.3 for x/y/z/t
    %   - mindwell: default = FPGA clock tick time, shouldn't be smaller
    %               than this
    %   - fixFOV: default = false; for datasets with scan offset changing
    %             during recording, keep the neurons fixed and move the FOV
    
    properties
        vr
        fa
        voxelsize = [1 1 2]; % unit: um, only used to create default assembler
        blursigma = [1.3 1.3 1.3 1.3]; % x/y/z/t Gaussian smoothing kernel size
        mindwell = 1/VolumeReconstructor.FPGA_CLOCK; % min valid dwell time per voxel
        fixFOV = false; % for datasets with scan offset changing during recording, keep the neurons fixed and move the FOV
        thisFrameInd
        thisFrame
        prevFrame
        nextFrame
        thisFrameSmoothed
    end
    
    methods
        
        function fb = HMMFrameBuffer(vr,op)
            fb.vr = vr;
            % update default options with input
            fn = fieldnames(fb);
            for j=1:length(fn)
                if isfield(op,fn{j})
                    fb.(fn{j}) = op.(fn{j});
                end
            end
            if length(fb.voxelsize)==1
                fb.voxelsize = repmat(fb.voxelsize,[1 3]);
            end
            if length(fb.blursigma)==1
                fb.blursigma = repmat(fb.blursigma,[1 4]);
            end
            if isempty(fb.fa)
                fb.set_default_assembler;
            end
            % initialize volumes as empty
            fb.thisFrameInd = [];
            fb.thisFrame = struct('rc',[],'rd',[],'gc',[],'gd',[]);
            fb.prevFrame = struct('rc',[],'rd',[],'gc',[],'gd',[]);
            fb.nextFrame = struct('rc',[],'rd',[],'gc',[],'gd',[]);
            fb.thisFrameSmoothed = struct('rc',[],'rd',[],'gc',[],'gd',[]);
        end
        
        function set_default_assembler(fb)
            fb.fa = fb.get_default_assembler;
        end
        
        function fa = get_default_assembler(fb)
            ze = -fb.vr.z_scale:fb.voxelsize(3):fb.vr.z_scale;
            if fb.fixFOV
                % get assembler range from vr imaging settings
                is = fb.vr.settings(end).imagingSettings;
                xe = (-is.xScanRange_Microns/2):fb.voxelsize(1):is.xScanRange_Microns/2;
                ye = (-is.yScanRange_Microns/2):fb.voxelsize(2):is.yScanRange_Microns/2;
            else
                % get assembler range from vr.scanBox
                xe = fb.vr.scanBox.globalBox(1):fb.voxelsize(1):fb.vr.scanBox.globalBox(2);
                ye = fb.vr.scanBox.globalBox(3):fb.voxelsize(2):fb.vr.scanBox.globalBox(4);
            end
            fa = FastAssembler(xe,ye,ze,[],false,fb.vr.z_scale);
        end
        
        function load_frame(fb,frameInd)
            if frameInd == fb.thisFrameInd % frame already loaded, do nothing
                return;
            end
            
            if frameInd == fb.thisFrameInd+1 % load next frame only
                fb.prevFrame = fb.thisFrame;
                fb.thisFrame = fb.nextFrame;
                fb.load_next_frame(frameInd);
            elseif frameInd == fb.thisFrameInd-1 % load prev frame only
                fb.nextFrame = fb.thisFrame;
                fb.thisFrame = fb.prevFrame;
                fb.load_prev_frame(frameInd);
            else
                fb.load_prev_frame(frameInd);
                fb.load_this_frame(frameInd);
                fb.load_next_frame(frameInd);
            end
            fb.get_smoothed_volume;
            fb.thisFrameInd = frameInd;
        end
        
        function load_prev_frame(fb,frameInd)
            % loading frameInd-1 volumes into fb.prevFrame
            if frameInd==1 % no prev frame for first frame, replicate (following imblur() default)
                trange = fb.vr.frame.edges(frameInd+[0 1]);
            else
                trange = fb.vr.frame.edges(frameInd-1+[0 1]);
            end
            [rc,rd] = fb.vr.binXYZ_fa(2,'r',trange,fb.fa,fb.fixFOV);
            [gc,gd] = fb.vr.binXYZ_fa(2,'g',trange,fb.fa,fb.fixFOV);
            fb.prevFrame.rc = rc;
            fb.prevFrame.rd = rd;
            fb.prevFrame.gc = gc;
            fb.prevFrame.gd = gd;
        end
        
        function load_this_frame(fb,frameInd)
            % loading frameInd volumes into fb.prevFrame
            trange = fb.vr.frame.edges(frameInd+[0 1]);
            [rc,rd] = fb.vr.binXYZ_fa(2,'r',trange,fb.fa,fb.fixFOV);
            [gc,gd] = fb.vr.binXYZ_fa(2,'g',trange,fb.fa,fb.fixFOV);
            fb.thisFrame.rc = rc;
            fb.thisFrame.rd = rd;
            fb.thisFrame.gc = gc;
            fb.thisFrame.gd = gd;
        end
        
        function load_next_frame(fb,frameInd)
            % loading frameInd+1 volumes into fb.prevFrame
            if frameInd==length(fb.vr.frame.edges)-1 % no next frame for last frame, replicate (following imblur() default)
                trange = fb.vr.frame.edges(frameInd+[0 1]);
            else
                trange = fb.vr.frame.edges(frameInd+1+[0 1]);
            end
            [rc,rd] = fb.vr.binXYZ_fa(2,'r',trange,fb.fa,fb.fixFOV);
            [gc,gd] = fb.vr.binXYZ_fa(2,'g',trange,fb.fa,fb.fixFOV);
            fb.nextFrame.rc = rc;
            fb.nextFrame.rd = rd;
            fb.nextFrame.gc = gc;
            fb.nextFrame.gd = gd;
        end
                
        function get_smoothed_volume(fb)
            rc = imblur(cat(4,fb.prevFrame.rc,fb.thisFrame.rc,fb.nextFrame.rc),fb.blursigma);
            rd = imblur(cat(4,fb.prevFrame.rd,fb.thisFrame.rd,fb.nextFrame.rd),fb.blursigma);
            rc = rc(:,:,:,2);
            rd = rd(:,:,:,2);
            gc = imblur(cat(4,fb.prevFrame.gc,fb.thisFrame.gc,fb.nextFrame.gc),fb.blursigma);
            gd = imblur(cat(4,fb.prevFrame.gd,fb.thisFrame.gd,fb.nextFrame.gd),fb.blursigma);
            gc = gc(:,:,:,2);
            gd = gd(:,:,:,2);
            valid = isfinite(rd) & isfinite(gd) & rd > fb.mindwell & gd > fb.mindwell;
            rc(~valid) = 0;
            rd(~valid) = 0;
            gc(~valid) = 0;
            gd(~valid) = 0;
            fb.thisFrameSmoothed.rc = rc;
            fb.thisFrameSmoothed.rd = rd;
            fb.thisFrameSmoothed.gc = gc;
            fb.thisFrameSmoothed.gd = gd;
        end
        
    end
    
    
end