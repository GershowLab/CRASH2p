classdef HyperMovieMaker_Projections < HyperMovieMaker
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties       
        mov_projection
       
    end

    methods

         function obj = get_default_proj_inds(obj)
            if isempty(obj.xi)
                obj.xi = 1:length(obj.mov_projection.xaxis);
            end
            if (isempty(obj.yi))
                obj.yi = 1:length(obj.mov_projection.yaxis);
            end
            if (isempty(obj.zi))
                obj.zi = 1:length(obj.mov_projection.zaxis);
            end
        end

        function obj = HyperMovieMaker_Projections(mov_projection, tsr, opp)
            if (nargin < 3)
                opp = [];
            end
            if (nargin < 2)
                if (ischar(mov_projection))
                    fpath = mov_projection;
                    mov_projection = load(fpath);
                    mov_projection.fpath = fpath;
                else
                    try
                        fpath = mov_projection.fpath;
                    catch
                        error('please add fpath to mov_projection and recreate hmm');
                        
                    end
                end
                srcdir = fileparts(fileparts(fpath));
                [~,timestamp,~] = fileparts(srcdir);
                vr_name = [timestamp '_vr.mat'];
                load(fullfile(srcdir,vr_name), 'vr');
                vr = vr.changeFilename(srcdir);
                tsr = vr.tsr.conditionData();
                tsr = tsr.addBehaviorVideo();
                try
                    tsr = tsr.addSleapResult();
                catch me
                    disp(me.getReport());
                end
            end   
            obj@HyperMovieMaker(tsr,opp);
            obj.mov_projection = mov_projection;           
            obj.op.proj_options = mov_projection.op;
        end
                    
        function [im, uaxis, vaxis] = get_projection(obj, t, clr, prj_ind)
            %prj_ind: 1 = xy, 2 = xz, 3 = zy
            ind = 0*t;
            if (obj.mip(min(prj_ind, length(obj.mip))))
                flag = '_mip';
            else
                flag = '';
            end
            p = {'xy', 'xz', 'zy'};
            u = {'x','x','z'};
            v = {'y', 'z', 'y'};

            if (clr == 'r')
                clr = 'red';
            end
            if (clr == 'g')
                clr = 'green';
            end
            
            ui = {obj.xi, obj.xi, obj.zi};
            vi = {obj.yi, obj.zi, obj.yi};
            for j = 1:length(ind)
                [~,ind(j)] = min(abs(obj.mov_projection.frametimes-t(j)),[],'all','omitnan');
            end
            
            im = obj.mov_projection.([clr '_' p{prj_ind} flag])(vi{prj_ind}, ui{prj_ind}, ind(1));
            uaxis = obj.mov_projection.([u{prj_ind} 'axis'])(ui{prj_ind});
            vaxis = obj.mov_projection.([v{prj_ind} 'axis'])(vi{prj_ind});
            im = repmat(im, [1 1 length(t)]);
            for j = 2:length(t)
                im(:,:,j) = obj.mov_projection.([clr '_' p{prj_ind} flag])(vi{prj_ind}, ui{prj_ind}, ind(j));
            end
        end


        function tx = get_default_time_axis(obj)
            tx = sort(obj.mov_projection.frametimes(:));
        end

        function fpath = get_default_fname(obj)
            dstdir = fileparts(obj.tsr.tracker.filename);
            [~,fname] = fileparts(obj.mov_projection.op.savename);
            fpath = fullfile(dstdir, [fname '.mp4']);
            %obj.mov_projection.fname = [fname '.mp4'];

        end

       
      
        function obj = set_projection_box(obj)
            proj_overlay.x = NaN([1 2]);
            proj_overlay.y = NaN([1 2]);
            proj_overlay.z = NaN([1 2]);
            %if (obj.op.overlay_projection_lines)
            try
                [proj_overlay.x(1), proj_overlay.x(2)] = bounds(obj.mov_projection.xaxis(obj.mov_projection.op.xi));
            catch
                [proj_overlay.x(1), proj_overlay.x(2)] = bounds(obj.mov_projection.xaxis);
            end
            
            try
                [proj_overlay.y(1), proj_overlay.y(2)] = bounds(obj.mov_projection.yaxis(obj.mov_projection.op.yi));
            catch
                [proj_overlay.y(1), proj_overlay.y(2)] = bounds(obj.mov_projection.yaxis);
            end
            
            try
                [proj_overlay.z(1), proj_overlay.z(2)] = bounds(obj.mov_projection.zaxis(obj.mov_projection.op.zi));
            catch
                [proj_overlay.z(1), proj_overlay.z(2)] = bounds(obj.mov_projection.zaxis);
            end
            %end
            proj_overlay.box{1} = [proj_overlay.x([1 2 2 1 1]); proj_overlay.y([1 1 2 2 1])];
            proj_overlay.box{2} = [proj_overlay.x([1 2 2 1 1]); proj_overlay.z([1 1 2 2 1])];
            proj_overlay.box{3} = [proj_overlay.z([1 2 2 1 1]); proj_overlay.y([1 1 2 2 1])];
            obj.proj_overlay = proj_overlay;
            obj.proj_overlay.do_overlay = obj.op.overlay_projection_lines;
        end



    end
end