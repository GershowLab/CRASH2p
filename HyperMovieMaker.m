classdef HyperMovieMaker
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        axbehav
        axred
        axgreen
        axratio

        xi
        yi
        zi

        op

        tsr
        mip = false;
        
        time_range
        
        behav_xl
        behav_yl
        
        tx
        traj
        
        glim;
        rlim;
        alpharlim;
        ratiolim;
        
        tref = 0;
        
        gcp %green colorbar position
        gco %green colorbar orientation
        rcp
        rco

        outputname = ''

        vidr;

        firstrun;
        
        proj_overlay
    end

    methods(Abstract)

        % [im, uaxis, vaxis] = get_projection(obj, t, clr, prj_ind);
        % im < 2D matrix
        % uaxis, vaxis < 1D vectors, such that
        % pcolor(uaxis, vaxis, im) is a legal command
        %
        % obj < HyperMovieMaker
        % t < time (according to FPGA clock) at which to get projections
        % clr < color, either 'r'/'red', or 'g'/'green'
        % prj_ind < 1 = xy, 2 = xz, 3 = zy
        %
        % obj.xi, obj.yi, obj.zi are used to select ranges of full
        % projections to return
        [im, uaxis, vaxis] = get_projection(obj, t, clr, prj_ind);


        % tx = get_default_time_axis(obj);
        % tx < times at which to make a movie frame
        % ob < HyperMovieMaker
        tx = get_default_time_axis(obj);

        % obj = get_default_proj_inds(obj);
        % sets obj.xi, obj.yi, obj.zi to be the full range of the projected
        % image
        obj = get_default_proj_inds(obj);

        % fpath = default_fname(obj)
        % fpath < path to a file (.mp4 extension) generated from data
        % available to HyperMovieMaker object (e.g. from timestamp and
        % projection indices)
        fpath = get_default_fname(obj);
        
        obj = set_projection_box(obj);
    end

    methods
        function obj = HyperMovieMaker(tsr, opp)
            
            if nargin > 1 && iscell(opp)
                opp = struct(opp{:});
            end


            obj.op.ratio_white = true;
            obj.op.log_ratio = true;

            fn = fieldnames(obj.op);
            for j = 1:length(fn)
                try
                    obj.op.(fn{j}) = opp.(fn{j});
                catch
                end
            end

            obj.tsr = tsr;
            
            obj.op.cmap_r = inferno(256);
            obj.op.cmap_g = viridis(256);
            if (obj.op.ratio_white)
                obj.op.cmap_ratio = jet(256);
            else
                ht = hot(256);
                obj.op.cmap_ratio = ht(:,3:-1:1);
            end
            
            obj.op.fwdcolor = [0 .8 .6];
            obj.op.bckcolor = [.8 0 0.2];
            obj.op.othercolor = [0 0 0];
            
            obj.op.trajectory_scale_bar = 1; %mm
            obj.op.im_scale_bar = 50; %um
            obj.op.fignum = 1;            
                       
            obj.op.figwidth = 1280;
            obj.op.figheight = 720;
            
            obj.op.tfwd = [];
            obj.op.tbck = [];
            
            obj.op.gtext = {'Green','GCaMP/GFP'};
            obj.op.rtext = {'Red','mCherry'};

            obj.op.rattext = {'Green/Red'};
            obj.op.genotype = '';
            obj.op.show_trajectory = true;
            
            obj.op.mpind = 1;
            
            obj.op.reverse_z = false;
            
            obj.op.overlay_projection_lines = false;
            
            obj.op.max_behav_dim = 15e3; %micron
            
            obj.op.show_ratio = true;
            
            try
                fn = fieldnames(opp);
                for j = 1:length(fn)
                    obj.op.(fn{j}) = opp.(fn{j});
                end
            catch
            end
            
            obj.proj_overlay.do_overlay = false;

        end

        function obj = open_behavior_video(obj)
            obj.vidr = VideoReader(obj.tsr.behavior(1).videoFileName);
        end

        function [im,framenum, xaxis, yaxis] = get_behavior_image(obj,t)

                
            [im,framenum, xaxis, yaxis] = obj.tsr.getBehaviorImage(obj.vidr,t);
        end
            


       

        function obj = setup_figure(obj)
            
            obj = obj.open_behavior_video();
            obj = get_default_proj_inds(obj);

            obj.firstrun = true;

            fignum = obj.op.fignum;
            figure(fignum);
            close(fignum);
            figure(fignum);
            
            u = get(fignum, 'Units');
            set(fignum, 'Units','pixels');
            p = get(fignum, 'Position');
            
            ss = get(0,'ScreenSize');
            
            
            p([1 2]) = max(ss([3 4])/4 - [obj.op.figwidth obj.op.figheight]/4, 0);
            p([3 4]) = [obj.op.figwidth obj.op.figheight];
            
            set(fignum, 'Position', p);
            set(fignum, 'Units',u,'color','k');
            
            if (obj.op.show_ratio)
                vw = obj.op.figwidth/2;
                vh = obj.op.figheight/2;

                obj.axbehav = axes('Position',[0 .5 .5 .5]); %top left corner

                [~,x,y] = obj.get_projection(0,'r',1);
                [~,z] = obj.get_projection(0,'r',3);

                xr = diff(x([1 end]));
                yr = diff(y([1 end]));
                zr = diff(z([1 end]));

                ar = (xr + zr)/(yr + zr);

                h = min(vw/ar, vh)/obj.op.figheight;
                w = h*ar*obj.op.figheight/obj.op.figwidth;
                a = 0.005; %small spacing adjustment
                if (h > 0.45 && w > 0.45)
                    h = 0.45;
                    w = h*ar*obj.op.figheight/obj.op.figwidth;
                end


                if (vw - w*obj.op.figwidth) > (vh - h*obj.op.figheight)
                    %put scale bar to the side
                    sbw = 0.25*(0.5 - w);
                    sbh = .5 - 2*a;
                    obj.gcp = [.5-a-sbw a sbw sbh];
                    obj.gco = 'eastoutside';
                else 
                    %put scale bar below 
                    sbh = 0.25*(0.5 - h);
                    sbw = .4;% - 2*a;
                    obj.gcp = [.05 a sbw sbh];
                    obj.gco = 'southoutside';
                end
                obj.rcp = obj.gcp + [0.5 0 0 0];
                obj.rco = obj.gco;

                xw = w*xr/(xr + zr);
                zw = w*zr/(xr + zr);
                yh = h*yr/(yr + zr);
                zh = h*zr/(yr + zr);


                if (obj.op.ratio_white)
                    obj.axratio.background = axes('Position',[.5 .5 .5 .5], 'Color','w');
                    obj.axratio.background.XTick = [];
                    obj.axratio.background.YTick = [];
                    obj.axratio.background.XColor = 'w';
                    obj.axratio.background.YColor = 'w';
                end            

                obj.axratio.xy = axes('Position',[.5 1-yh xw yh]);
                obj.axratio.xz = axes('Position',[.5 1-h xw zh]);
                obj.axratio.zy = axes('Position',[.5+xw 1-yh zw yh]);

                obj.axred.xy = axes('Position',obj.axratio.xy.Position - [0 0.5 0 0]);
                obj.axred.xz = axes('Position',obj.axratio.xz.Position - [0 0.5 0 0]);
                obj.axred.zy = axes('Position',obj.axratio.zy.Position - [0 0.5 0 0]);

                obj.axgreen.xy = axes('Position',obj.axred.xy.Position - [0.5 0 0 0]);
                obj.axgreen.xz = axes('Position',obj.axred.xz.Position - [0.5 0 0 0]);
                obj.axgreen.zy = axes('Position',obj.axred.zy.Position - [0.5 0 0 0]);
                if obj.op.ratio_white
                    ratcolor = 'k';
                else
                    ratcolor = 'w';
                end
                annotation("textbox",[.5+xw 1-h zw zh],'String',obj.op.rattext,'Color',ratcolor,'EdgeColor','none','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
                annotation("textbox",[xw .5-h zw zh],'String',obj.op.gtext,'Color','g','EdgeColor','none','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
                annotation("textbox",[.5+xw .5-h zw zh],'String',obj.op.rtext,'Color','r','EdgeColor','none','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');

                annotation("line", [.5 .5],[0 1], 'Color','w','LineWidth',2);
                annotation("line", [0 1],[.5 .5], 'Color','w','LineWidth',2);
            else
                vw = obj.op.figwidth/2;
                vh = obj.op.figheight/2;

                obj.axbehav = axes('Position',[0 0 .5 1]); %left side

                [~,x,y] = obj.get_projection(0,'r',1);
                [~,z] = obj.get_projection(0,'r',3);

                xr = diff(x([1 end]));
                yr = diff(y([1 end]));
                zr = diff(z([1 end]));

                ar = (xr + zr)/(yr + zr);

                h = min(vw/ar, vh)/obj.op.figheight;
                w = h*ar*obj.op.figheight/obj.op.figwidth;
                a = 0.005; %small spacing adjustment
                if (h > 0.45 && w > 0.45)
                    h = 0.45;
                    w = h*ar*obj.op.figheight/obj.op.figwidth;
                end


                if (vw - w*obj.op.figwidth) > (vh - h*obj.op.figheight)
                    %put color bar to the side
                    sbw = 0.25*(0.5 - w);
                    sbh = .5 - 2*a;
                    obj.gcp = [1-a-sbw .5+a sbw sbh];
                    obj.gco = 'eastoutside';
                else 
                    %put color bar below 
                    sbh = 0.25*(0.5 - h);
                    sbw = .4;% - 2*a;
                    obj.gcp = [.505 .5+a sbw sbh];
                    obj.gco = 'southoutside';
                end
                obj.rcp = obj.gcp + [0 -0.5 0 0];
                obj.rco = obj.gco;

                xw = w*xr/(xr + zr);
                zw = w*zr/(xr + zr);
                yh = h*yr/(yr + zr);
                zh = h*zr/(yr + zr);
         

                obj.axgreen.xy = axes('Position',[.5 1-yh xw yh]);
                obj.axgreen.xz = axes('Position',[.5 1-h xw zh]);
                obj.axgreen.zy = axes('Position',[.5+xw 1-yh zw yh]);

                obj.axred.xy = axes('Position',obj.axgreen.xy.Position - [0 0.5 0 0]);
                obj.axred.xz = axes('Position',obj.axgreen.xz.Position - [0 0.5 0 0]);
                obj.axred.zy = axes('Position',obj.axgreen.zy.Position - [0 0.5 0 0]);

                annotation("textbox",[.5+xw 1-h zw zh],'String',obj.op.gtext,'Color','g','EdgeColor','none','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
                annotation("textbox",[.5+xw .5-h zw zh],'String',obj.op.rtext,'Color','r','EdgeColor','none','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');

                annotation("line", [.5 .5],[0 1], 'Color','w','LineWidth',2); %[x0 x1], [y0 y1]
                annotation("line", [.5 1],[.5 .5], 'Color','w','LineWidth',2);

            end

            
                
            
        end


        function obj = set_limits(obj, tx)
            if (nargin < 2)
                tx = obj.get_default_time_axis;
            end


            [im,~,xaxis,yaxis] = obj.get_behavior_image(tx(1));

            vxr = abs(diff(xaxis([1 end])));
            vyr = abs(diff(yaxis([1 end])));

            sl = interp1(obj.tsr.tx, obj.tsr.tracker.stageloc(1:2,:)', tx)';
            
            xl = [min(sl(1,:)) max(sl(1,:))] + vxr*[-.5 .5];
            yl = [min(sl(2,:)) max(sl(2,:))] + vyr*[-.5 .5];

            try
                pcolor(obj.axbehav, xaxis, yaxis, im); 
            catch
                obj = obj.setup_figure();
                pcolor(obj.axbehav, xaxis, yaxis, im);
            end
            
            shading(obj.axbehav, 'flat'); colormap(obj.axbehav, gray(256)); obj.axbehav.Color = 'k'; obj.axbehav.XTick = []; obj.axbehav.YTick = [];
            obj.axbehav.XLim = xl;
            obj.axbehav.YLim = yl;
            axis(obj.axbehav, 'equal');
            abhw = diff(obj.axbehav.XLim);
            abhh = diff(obj.axbehav.YLim);
            
            %abxyrat = abhh/abhw;
            
            obj.axbehav.XLim = mean(xl) + abhw*[-.5 .5];
            obj.axbehav.YLim = mean(yl) + abhh*[-.5 .5];
            % axis(axbehav, 'equal');
            
            labeledBar(obj.axbehav, [.95 .05]*obj.axbehav.XLim'+[0 obj.op.trajectory_scale_bar*1000], [.05 .95]*obj.axbehav.YLim' + [0 0], [num2str(obj.op.trajectory_scale_bar) ' mm'],'below','w',{'LineWidth',6},{'FontSize',14});

            
            obj.behav_xl = obj.axbehav.XLim;
            obj.behav_yl = obj.axbehav.YLim;
            obj.tx = tx;
            obj.traj = interp1(obj.tsr.tx, obj.tsr.neuron(1).loc(1:2,:)', tx)';
            if (isempty(obj.rlim))
                rpct = percentile(obj.get_projection(tx, 'r', 1), [.25 .85 .99]);
                obj.rlim = [0 rpct(3)];
                obj.alpharlim = rpct([1 2]);
            end
            if (isempty(obj.glim))
                obj.glim = [0 percentile(obj.get_projection(tx, 'g', 1), .99)];
            end
            if (iscell(obj.glim))
                gl = obj.glim{1};
            else
                gl = obj.glim;
            end
            if (iscell(obj.rlim))
                rl = obj.rlim{1};
            else
                rl = obj.rlim;
            end

            if (obj.op.log_ratio)
                
                obj.ratiolim = max(0,ceil(log10(gl(2)./rl(2)))) + [-1 0];
            else
                obj.ratiolim = [0 gl(2)./rl(2)];
            end
            

            
        end
        
        
        function plot_behavior(obj, t)
            [im,~, xaxis, yaxis] = obj.get_behavior_image(t);
            pcolor(obj.axbehav, xaxis, yaxis, im); shading(obj.axbehav, 'flat'); colormap(obj.axbehav, gray(256)); obj.axbehav.Color = 'k'; obj.axbehav.XTick = []; obj.axbehav.YTick = [];
            obj.axbehav.XLim = obj.behav_xl;
            obj.axbehav.YLim = obj.behav_yl;
            if (obj.op.show_trajectory)
                
                hold (obj.axbehav, 'on');
                ind = interp1(obj.tx, 1:length(obj.tx), t, 'nearest', 'extrap');
                plot(obj.axbehav, obj.traj(1,1:ind), obj.traj(2,1:ind), 'm.', 'MarkerSize',7);
                hold (obj.axbehav,'off');
            end
        end
        
        function obj = plot_activity(obj,t)
            
            p = {'xy', 'xz', 'zy'};
            for j = 1:3
                [imr,uaxis,vaxis] = obj.get_projection(t,'r',j);
                img = obj.get_projection(t,'g',j);
                
                pcolor(obj.axred.(p{j}), uaxis, vaxis, imr); 
                shading (obj.axred.(p{j}), 'flat');
                axis(obj.axred.(p{j}), 'off', 'tight');
                colormap(obj.axred.(p{j}), obj.op.cmap_r);
                if (iscell(obj.rlim))
                    obj.axred.(p{j}).CLim = obj.rlim{j};
                else
                    obj.axred.(p{j}).CLim = obj.rlim;
                end
                if (obj.proj_overlay.do_overlay)
                    hold (obj.axred.(p{j}), 'all');
                    plot(obj.axred.(p{j}), obj.proj_overlay.box{j}(1,:),  obj.proj_overlay.box{j}(2,:), 'w--');
                    hold (obj.axred.(p{j}), 'off');
                end
            
                             
                pcolor(obj.axgreen.(p{j}), uaxis, vaxis, img);
                shading (obj.axgreen.(p{j}), 'flat');
                axis(obj.axgreen.(p{j}), 'off', 'tight');
                colormap(obj.axgreen.(p{j}), obj.op.cmap_g);
                if (iscell(obj.glim))
                    obj.axgreen.(p{j}).CLim = obj.glim{j};
                else
                    obj.axgreen.(p{j}).CLim = obj.glim;
                end
                if (obj.proj_overlay.do_overlay)
                    hold (obj.axgreen.(p{j}), 'all');
                    plot(obj.axgreen.(p{j}), obj.proj_overlay.box{j}(1,:),  obj.proj_overlay.box{j}(2,:), 'w--');
                    hold (obj.axgreen.(p{j}), 'off');
                end

                if (obj.op.show_ratio)
                    if (obj.op.log_ratio)
                        imratio = log10(img)-log10(imr);
                    else
                        imratio = img./imr;
                    end
                    srf = pcolor(obj.axratio.(p{j}), uaxis, vaxis, imratio);
                    shading (obj.axratio.(p{j}), 'flat');
                    srf.AlphaData = max(min(1,(imr - obj.alpharlim(1))/diff(obj.alpharlim)), 0);
                    srf.FaceAlpha = 'flat';

                    axis(obj.axratio.(p{j}), 'off', 'tight');
                    colormap(obj.axratio.(p{j}), obj.op.cmap_ratio);
                    obj.axratio.(p{j}).CLim = obj.ratiolim;
                    if (obj.proj_overlay.do_overlay)
                        hold (obj.axratio.(p{j}), 'all');
                        plot(obj.axratio.(p{j}), obj.proj_overlay.box{j}(1,:),  obj.proj_overlay.box{j}(2,:), 'w--');
                        hold (obj.axratio.(p{j}), 'off');
                    end
                end
                


                if (isempty(obj.firstrun) || obj.firstrun)
                    axis(obj.axred.(p{j}), 'equal');
                    axis(obj.axgreen.(p{j}), 'equal');
                    if (obj.op.show_ratio)
                        axis(obj.axratio.(p{j}), 'equal');

                    end

                end
                
            end
            if (isempty(obj.firstrun) || obj.firstrun)
                if (obj.op.show_ratio)
                    if (obj.gco(1) == 'e') %colorbar to the right
                        sbax = obj.axratio.zy;
                    else
                        sbax = obj.axratio.xy;
                    end
                else
                    if (obj.gco(1) == 'e') %colorbar to the right
                        sbax = obj.axgreen.zy;
                    else
                        sbax = obj.axgreen.xy;
                    end
                end
                if (obj.gco(1) == 'e')
                     if (obj.op.reverse_z)
                        sbx = sbax.XLim(1) + [-obj.op.im_scale_bar 0] - 5;
                     else
                         sbx = sbax.XLim(2) + 5+ [0 obj.op.im_scale_bar];
                     end
                     sby = sbax.YLim(2) - 5 + [0 0];
                else
                    sbx = sbax.XLim(1) + 5 + [0 obj.op.im_scale_bar];
                    sby = sbax.YLim(2) +10 + [0 0];
                end
                labeledBar(sbax, sbx,sby, sprintf('%d {\\mu}m',obj.op.im_scale_bar), 'below', 'w',{'LineWidth',6}, {'Interpreter', 'Tex','FontSize', 14});
            end


            if (obj.op.reverse_z)
                obj.axred.zy.XDir = 'reverse';
                obj.axgreen.zy.XDir = 'reverse';
                if (obj.op.show_ratio)
                    obj.axratio.zy.XDir = 'reverse';
                end
            else
                obj.axred.xz.YDir = 'reverse';
                obj.axgreen.xz.YDir = 'reverse';
                if (obj.op.show_ratio)
                    obj.axratio.xz.YDir = 'reverse';
                end
            end
                
            obj.firstrun = false;

        end
        
        function obj = plot_one_frame(obj,t)
            if (obj.op.ratio_white)
                ratcolor = 'k';
            else
                ratcolor = 'w';
            end
            obj.plot_behavior(t);
            obj = obj.plot_activity(t);
            if (~isempty(obj.op.tfwd) && any(t >= obj.op.tfwd(:,1) & t <= obj.op.tfwd(:,2)) )
                text([.05 .95]*obj.axbehav.XLim', [.05 .95]*obj.axbehav.YLim',{'forward','crawl'},'Color',obj.op.fwdcolor, 'Parent', obj.axbehav, 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);
            end
            if (~isempty(obj.op.tbck) && any(t >= obj.op.tbck(:,1) & t <= obj.op.tbck(:,2)) )
                text([.05 .95]*obj.axbehav.XLim', [.05 .95]*obj.axbehav.YLim',{'backward','crawl'},'Color',obj.op.bckcolor, 'Parent', obj.axbehav, 'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);
            end
            if (~isempty(obj.op.genotype))
                text([.95 .05]*obj.axbehav.XLim', [.95 .05]*obj.axbehav.YLim',obj.op.genotype,'Color','w', 'Parent', obj.axbehav, 'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
            end
            text([.05 .95]*obj.axbehav.XLim', [.95 .05]*obj.axbehav.YLim',sprintf('t = %.1f s', t - obj.tref),'Color','w', 'Parent', obj.axbehav, 'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',12);

            if (iscell(obj.glim))
                gl = obj.glim{1};
            else
                gl = obj.glim;
            end
            tcks = linspace(gl(1),gl(2),3);
            clear tlbl;
            for j = 1:length(tcks)
                tlbl{j} = sprintf('%.1f', tcks(j)/1e6); %#ok<AGROW> 
            end
            tlbl{end} = [tlbl{end} ' MHz'];


            colorbar(obj.axgreen.xy,obj.gco,'Position',obj.gcp,'Color','w','AxisLocation','in','FontSize',8,'Ticks',tcks, 'TickLabels', tlbl);
            if (iscell(obj.rlim))
                rl = obj.rlim{1};
            else
                rl = obj.rlim;
            end
            tcks = linspace(rl(1),rl(2),3);            
            clear tlbl;
            for j = 1:length(tcks)
                tlbl{j} = sprintf('%.1f', tcks(j)/1e6); %#ok<AGROW> 
            end
            tlbl{end} = [tlbl{end} ' MHz'];

            colorbar(obj.axred.xy,obj.rco,'Position',obj.rcp,'Color','w','AxisLocation','in','FontSize',8,'Ticks',tcks, 'TickLabels', tlbl);
            if (obj.op.log_ratio)
                tcks = log10(logspace(obj.ratiolim(1), obj.ratiolim(2), 3));
                clear tlbl;

                for j = 1:length(tcks)
                    tlbl{j} = sprintf('%.1f', 10.^tcks(j)); %#ok<AGROW> 
                end
            else
                tcks = linspace(obj.ratiolim(1), obj.ratiolim(2), 3);
                clear tlbl;
                for j = 1:length(tcks)
                    tlbl{j} = sprintf('%.1f', tcks(j)); %#ok<AGROW> 
                end
            end
            if (obj.op.show_ratio)
                colorbar(obj.axratio.xy,obj.rco,'Position',obj.rcp+[0 0.5 0 0],'Color',ratcolor,'AxisLocation','in','FontSize',8,'Ticks',tcks, 'TickLabels', tlbl);
            end
        end
        
        function obj = make_movie(obj, tx)
            obj = obj.setup_figure();
            obj = obj.set_limits(tx); %consider separating out
            obj = obj.set_projection_box();
            obj.play_movie(tx);
        end

        function play_movie(obj, tx)
            
            
            writeFrameOrPause(); %reset internal storage
            vidobj = [];
            if (~isempty(obj.outputname))
                try
                    vidobj = VideoWriter(obj.outputname, 'MPEG-4');
                    vidobj.open();
                    [dstdir, fname] = fileparts(obj.outputname);
                    options = obj.op;
                    save(fullfile(dstdir, [fname '_options.mat']), 'options');
                catch
                    vidobj = [];
                    disp ('did not create video writer')
                end
            end
            for j = 1:length(tx)
                obj = obj.plot_one_frame(tx(j));
                writeFrameOrPause(vidobj);
            end
            writeFrameOrPause(vidobj,[],[],[], true);
        end

        %function obj = set_output_name(obj, fpath) - use fname
        %function obj = set_output_name(obj, dstdir) - use default name but place in
        %dstdir
        %function obj = set_output_name(obj) - use default name and dstdir
        function obj = set_output_name(obj, fname)
            if (nargin < 2 || isempty(fname))
                fname = obj.get_default_fname();
            else
                if (exist(fname, 'dir'))
                    dstdir = fname;
                    fname = obj.get_default_fname();
                    [~,fname] = fileparts(fname);
                    fname = fullfile(dstdir, [fname '.mp4']);
                end
            end
                
            obj.outputname = fname;
        end
      


    end
end