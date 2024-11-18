function [ah,op,annots] = activity_summary_figure_vnc (figData, opp)
    op.cmap_template = hot(256);
    op.cmap_r = inferno(256);%inferno(256);
    op.cmap_g = viridis(256);%inferno(256);

    op.cmap_rat = parula(256);%matlab changed default colormap from jet to parula in 2014! I learned that today 8/15/23! %op.cmap_r(:,[2 1 3]);
    op.rclim = [];
    op.gclim = [];

    op.bar_width = 25; %scale bar width
    op.mar_width = 0.25; % unit: ratio of panel width
    op.sbl = 'right';

    op.fwdcolor = [0 .8 .6];
    op.bckcolor = [.8 0 0.2];
    op.othercolor = [0 0 0];

    op.trajectory_scale_bar = 500; %um
    op.fignum = 1;
    op.activity_overlap = 0.5;
    op.activity_range = [0.5 4];
    op.norm_rat_range = [0.5 4];

    op.activity_linewidth_zoom = 1.5;

    op.stride_marker_linewidth = 1.5;
    
    op.data_channel = 'ratio';

    if (nargin > 0)
       op.voi_colors = winter(ceil(1.1*numel(figData.vois)))*0.7;%viridis(numel(figData.vois));
       op.dobck = ~isempty(figData.parameters.bckTrange);
    else
        ah = op;
        return;
    end

    try
        fn = fieldnames(opp);
        for j = 1:length(fn)
            op.(fn{j}) = opp.(fn{j});
        end
    catch
    end


    if (isempty(op.gclim))
        op.gclim = [0 percentile(figData.xtproj.grate(:), .99)];
    end
    if (isempty(op.rclim))
        op.rclim = [0 percentile(figData.xtproj.rrate(:), .99)];
    end

    [ad,po] = blank8x10Figure(op.fignum, 'nrows',5,'fontsize',12);

    annots = [];
    %%stacked activity ratio full time range
    an = 1;


    cbwidth = (ad.dh-ad.h)/2;

    pos(1) = (ad.lx4 + ad.w4 + ad.clx4)/2;%+ad.w4;
    pos(2) = ad.h0-ad.h;
    pos(3) = ad.rx4+ad.w4-pos(1) - 2*cbwidth;
    pos(4) = ad.h;
    ah(an).pos = pos;
    ah(an).name = 'activity large time window';
    switch lower(op.data_channel)
        case {'green', 'g', 'grn'}
            act_plot_y = figData.activityPlot.green./figData.activityPlot.green0';
        case {'red', 'r', 'rd'}
            act_plot_y = figData.activityPlot.red./figData.activityPlot.red0';
        case {'red_uc', 'r_uc'}
            act_plot_y = figData.activityPlot.red_uc./figData.activityPlot.red_uc0';
        otherwise
            act_plot_y = figData.activityPlot.y;
    end
    handles = stackedActivityPlots(figData.activityPlot.tx,act_plot_y, op.voi_colors, op.activity_overlap, ah(an).pos, op.activity_range);

    set(handles.p, 'LineWidth',1);
    set(handles.lines, 'LineWidth',1);
    set(handles.text, po.labelOptions{:});
    ah(an).ax = handles.ax;
    ah(an).h = rmfield(handles,'ax');
    act_x_lim = get(ah(an).ax(1), 'XLim');
    for k = 1:length(handles.ax)
        text(ah(an).ax(k), ah(an).ax(k).XLim(2), ah(an).ax(k).YLim(1),num2str(k),'Color', op.voi_colors(k,:), 'FontSize',po.fontsize,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','bottom');
    end


     ygap =  (ad.dh-ad.h)/2;
     ytop = pos(2) - 0.15 * ygap;
     ybottom = pos(2) - 0.85*ygap;
     lw = 0.5;
     for j = 1:size(figData.parameters.tfwd, 1)
         figx = dsxy2figxy_marc(ah(an).ax(1), figData.parameters.tfwd(j,:), ah(an).ax(1).YLim);
         annots = [annots annotation("textbox", [figx(1) ybottom diff(figx) ytop-ybottom], "String",'F', 'BackgroundColor', op.fwdcolor, 'Color','w','EdgeColor',op.fwdcolor','FontSize', po.fontsize-2, ...
             'FontWeight','bold','FontName',po.font,'LineWidth',lw, 'HorizontalAlignment','center','VerticalAlignment','middle')];
     end
     if (op.dobck)
         for j = 1:size(figData.parameters.tbck, 1)
             figx = dsxy2figxy_marc(ah(an).ax(1), figData.parameters.tbck(j,:), ah(an).ax(1).YLim);
             annots = [annots annotation("textbox", [figx(1) ybottom diff(figx) ytop-ybottom], "String",'B', 'BackgroundColor', op.bckcolor, 'Color','w','EdgeColor',op.bckcolor,'FontSize', po.fontsize-2,...
                 'FontWeight','bold','FontName',po.font,'LineWidth',lw, 'HorizontalAlignment','center','VerticalAlignment','middle')];
             dy = ytop-ybottom;
             nx = ceil(diff(figx)/dy);
             dx = diff(figx)/nx;
             for k = 1:nx
                 annots = [annots annotation('line',figx(1)+(k-1)*dx + [0 dx], [ybottom, ytop], 'LineStyle','-','LineWidth',2,'Color','w')];
             end
         end
     end

    %%template with vois

    an = an+1;
    pos = ah(an-1).pos;
    pos(1) = ad.lx5;
    pos(3) = ad.w4;
    pos(2) = pos(2)-ad.h;
    pos(4) = ad.h*2;
    ah(an).pos = pos;
    ah(an).name = 'template with vois';
    ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});


    bwim = Voi.makeLabeledMask(figData.vois); %bwim is shortand for black and white (aka logical) image

    irng = find(any(bwim,[2 3]),1,'first'):find(any(bwim,[2 3]),1,'last');
    jrng = find(any(bwim,[1 3]),1,'first'):find(any(bwim,[1 3]),1,'last');
    krng = find(any(bwim,[1 2]),1,'first'):find(any(bwim,[1 2]),1,'last');
    
    pcolor(ah(an).ax, figData.templatey, figData.templatex,  max(figData.templateim(:,:,krng),[],3)); shading interp; colormap (ah(an).ax, op.cmap_template); hold on; axis equal; axis off;
    textx =  figData.templatey(end); %size(figData.templateim,2) + 1;
    %textx = max(jrng) + 7
    for k = 1:length(figData.vois)
        B = bwboundaries(imdilate(any(figData.vois{k}.makeMask,3), ones(3)));                
        
        for q = 1:length(B)
            plot (ah(an).ax,  figData.templatey(B{q}(:,2)), figData.templatex(B{q}(:,1)), 'Color', 'w', 'LineWidth',2);
            plot (ah(an).ax,  figData.templatey(B{q}(:,2)), figData.templatex(B{q}(:,1)), 'Color', op.voi_colors(k,:), 'LineWidth',1.5);
        end
        texty = mean(figData.templatex(B{1}(:,1)));
        text(ah(an).ax, textx, texty,num2str(k),'Color', op.voi_colors(k,:), 'FontSize',po.fontsize+2,'FontWeight','bold','HorizontalAlignment','left');
    end

    plot(ah(an).ax, figData.xtproj.boxregion(2,:), figData.xtproj.boxregion(1,:),'w:','LineWidth',2);



    %%xz projection
    pos = dsxy2figxy_marc(ah(an).ax(1), [figData.templatey(1) figData.templatex(1) diff(figData.templatey([1 end])) diff(figData.templatex([1 end]))]);
    h = diff(figData.templatez([end, 1]))/diff(figData.templatex([end,1])) * pos(3);
    pos(2) = pos(2)-h; pos(4) = h;
    ah(an).ax(2) = axes('Position',pos, po.axesopts{:});

    pcolor(ah(an).ax(2), figData.templatey, figData.templatez,  squeeze(max(figData.templateim,[],1))'); shading interp; colormap (ah(an).ax(2), op.cmap_template); hold(ah(an).ax(2), 'on');  axis off; %axis equal;
    plot(ah(an).ax(2), min(figData.xtproj.boxregion(2,:))*[1 1 0 0 1] +  max(figData.xtproj.boxregion(2,:))*[0 0 1 1 0], min(figData.xtproj.zrange)*[1 0 0 1 1] + max(figData.xtproj.zrange)*[0 1 1 0 0],'w:','LineWidth',2);

    for k = 1:length(figData.vois)
        B = bwboundaries(imdilate(squeeze(any(figData.vois{k}.makeMask,1))', ones(3)));
        for q = 1:length(B)
            plot (ah(an).ax(2),  figData.templatey(B{q}(:,2)), figData.templatez(B{q}(:,1)), 'Color', 'w', 'LineWidth',2);
            plot (ah(an).ax(2),  figData.templatey(B{q}(:,2)), figData.templatez(B{q}(:,1)), 'Color', op.voi_colors(k,:), 'LineWidth',1);
        end
    end

    labeledBar(ah(an).ax(1), ah(an).ax(1).XLim(1)+[0 op.bar_width], ah(an).ax(1).YLim(2) + [5 5], sprintf('%d {\\mu}m', op.bar_width), 'above', 'k', {'LineWidth',4},{'Interpreter','tex'});


    %%tx projection
    an = an+1;
    pos = ah(an-2).pos - [0 (ad.dh+ad.h)/2 0 0];
    ah(an).pos = pos;
    ah(an).name = 'tx projection';
    ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});

    norm_rat_proj = figData.xtproj.ratio./median(figData.xtproj.ratio,'all','omitnan');

    switch lower(op.data_channel)
        case {'green', 'g', 'grn'}
            xt_proj_im = figData.xtproj.grate;
            xt_proj_lim = op.gclim;    
            xt_proj_map = op.cmap_g;
        case {'red', 'r', 'rd'}
            xt_proj_im = figData.xtproj.rrate;
            xt_proj_lim = op.rclim;  
            xt_proj_map = op.cmap_r;
        case {'red_uc', 'r_uc'}
            xt_proj_im = figData.xtproj.rrate; %TODO need rrate_uc
            xt_proj_lim = op.rclim;  
            xt_proj_map = op.cmap_r;

        otherwise
            xt_proj_im = norm_rat_proj;
            xt_proj_lim = op.norm_rat_range;  
            xt_proj_map = op.cmap_rat;

    end
    


    pcolor(ah(an).ax, figData.xtproj.tx, figData.xtproj.px, xt_proj_im);
    shading(ah(an).ax, 'flat');
    set(ah(an).ax,'Clim',xt_proj_lim,'XLim',act_x_lim);    
    axis(ah(an).ax, 'off');

    colormap(ah(an).ax, xt_proj_map);

    alpha = 0;

    lw = 5;
    [fwdarrow_start_x, fwdarrow_start_y] = dsxy2figxy_marc(ah(an).ax(1), figData.parameters.fwdTrange, (1+alpha)*min(ah(an).ax(1).YLim) - alpha*max(ah(an).ax(1).YLim) + [0 0]);
    fwdarrow_start_y = fwdarrow_start_y-lw/(11*72*2);
    annots = [annots annotation("line",fwdarrow_start_x, fwdarrow_start_y, 'Color', op.fwdcolor,'LineWidth',lw)];

    if (op.dobck)
        [bckarrow_start_x, bckarrow_start_y] = dsxy2figxy_marc(ah(an).ax(1), figData.parameters.bckTrange, (1+alpha)*min(ah(an).ax(1).YLim) - alpha*max(ah(an).ax(1).YLim) + [0 0]);
        bckarrow_start_y =  bckarrow_start_y-lw/(11*72*2);
        annots = [annots annotation("line",bckarrow_start_x, bckarrow_start_y, 'Color', op.bckcolor,'LineWidth',lw)];
    end

    cbpos = pos;
    cbpos(3) = cbwidth;
    cbpos(1) = pos(1) + pos(3) + cbwidth;
    ah(an).ax(2) = colorbar(ah(an).ax, 'manual', 'Position', cbpos, 'FontSize', po.fontsize, 'FontName', po.font);


 
    %%forward crawling exceprt


    right_axes_location = ad.rx2; %0.6 * (ad.lx2 + ad.w2) + 0.4 * ad.rx2;

    

    an = an+1;
    pos = ah(an-1).pos;
    if (~op.dobck || isempty(figData.parameters.bckTrange) || (figData.parameters.fwdTrange(1) < figData.parameters.bckTrange(1)))
        pos(1) = ad.lx2;
    else
        pos(1) = right_axes_location;
    end
    pos(2) = pos(2) - 0.5*(ad.h+ad.dh);
    pos(3) = ad.w2;
    ah(an).pos = pos;
    ah(an).name = 'fwd crawl';

    
    
    inds = figData.activityPlot.tx >= figData.parameters.fwdTrange(1) & figData.activityPlot.tx <= figData.parameters.fwdTrange(2);
    handles = stackedActivityPlots(figData.activityPlot.tx(inds),act_plot_y(:,inds), op.voi_colors, op.activity_overlap, ah(an).pos,op.activity_range);
    set(handles.lines, 'LineWidth',1);

    set(handles.p, 'LineWidth',op.activity_linewidth_zoom);

    set(handles.text, po.labelOptions{:});
    ah(an).ax = handles.ax;
    ah(an).h = rmfield(handles,'ax');
    
    %dashed verticle lines at t=0 points
    tc = figData.alignedActivity.tcross(figData.alignedActivity.tcross > figData.parameters.fwdTrange(1) & figData.alignedActivity.tcross < figData.parameters.fwdTrange(2));
    for j = 1:length(ah(an).ax)
        hold(ah(an).ax(j),'on');            
        for k = 1:length(tc)
            plot(ah(an).ax(j),tc(k) + [0 0], ah(an).ax(j).YLim, 'k:', 'LineWidth', op.stride_marker_linewidth);
        end
    end
            
    [arr_x, arr_y] = dsxy2figxy_marc(ah(an).ax(end), figData.parameters.fwdTrange, ah(an).ax(end).YLim(2)+[0 0]);

    for j = 1:2
        ah(an).h.arrow(j) = annotation("arrow", [fwdarrow_start_x(j) arr_x(j)], [fwdarrow_start_y(j) arr_y(j)],'Color',op.fwdcolor);
    end

    %%fwd crawling exceprt projection
    an = an+1;
    pos = ah(an-1).pos - [0 ad.h 0 0];
    ah(an).pos = pos;
    ah(an).name = 'fwd tx projection';
    ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
    pcolor(ah(an).ax, figData.xtproj.tx, figData.xtproj.px, xt_proj_im);
    shading(ah(an).ax, 'flat');
    %     pcolor(ah(an).ax, figData.xtproj.tx, figData.xtproj.xaxis, lowpass1D(figData.xtproj.ratio,1)./median(figData.xtproj.ratio,'omitnan'));
%     shading(ah(an).ax, 'flat')
%     set(ah(an).ax,'Clim',[0.5 3]);
%     
    set(ah(an).ax,'XLim',figData.parameters.fwdTrange,'Clim',xt_proj_lim);
    axis(ah(an).ax, 'off');

    colormap(ah(an).ax, xt_proj_map);


    

    tc = figData.alignedActivity.tcross(figData.alignedActivity.tcross > figData.parameters.fwdTrange(1) & figData.alignedActivity.tcross < figData.parameters.fwdTrange(2));
    
    hold(ah(an).ax,'on');
    for k = 1:length(tc)
        plot(ah(an).ax,tc(k) + [0 0], ah(an).ax.YLim, 'w:', 'LineWidth', op.stride_marker_linewidth);
    end
    
   % ah(an).ax(2) = colorbar(ah(an).ax, 'manual', 'Position', cbpos);

    if (op.dobck && numel(figData.parameters.bckTrange) >= 2)

        %backward crawling excerpt
    
        an = an+1;
        pos = ah(an-2).pos;
        if (isempty(figData.parameters.bckTrange) || figData.parameters.fwdTrange(1) < figData.parameters.bckTrange(1))
            pos(1) = right_axes_location;
        else
            pos(1) = ad.lx2;
        end
        ah(an).pos = pos;
        ah(an).name = 'bck crawl';

        inds = figData.activityPlot.tx >= figData.parameters.bckTrange(1) & figData.activityPlot.tx <= figData.parameters.bckTrange(2);
        if (any(inds))
            handles = stackedActivityPlots(figData.activityPlot.tx(inds),act_plot_y(:,inds), op.voi_colors, op.activity_overlap, ah(an).pos,op.activity_range);
            delete(handles.lines);
            delete(handles.text);

            set(handles.p, 'LineWidth',op.activity_linewidth_zoom);

            ah(an).ax = handles.ax;
            for j = 1:length(ah(an).ax)
                ah(an).ax(j).YLim = ah(an-2).ax(j).YLim;
            end

            ah(an).h = rmfield(handles,{'ax','lines','text'});

            %dashed verticle lines at t=0 points
            tc = figData.alignedActivity.tcross(figData.alignedActivity.tcross > figData.parameters.bckTrange(1) & figData.alignedActivity.tcross < figData.parameters.bckTrange(2));
            for j = 1:length(ah(an).ax)
                hold(ah(an).ax(j),'on');
                for k = 1:length(tc)
                    plot(ah(an).ax(j),tc(k) + [0 0], ah(an).ax(j).YLim, 'k:', 'LineWidth', op.stride_marker_linewidth);
                end
            end

            [arr_x, arr_y] = dsxy2figxy_marc(ah(an).ax(end), figData.parameters.bckTrange, ah(an).ax(end).YLim(2)+[0 0]);

            for j = 1:2
                ah(an).h.arrow(j) = annotation("arrow", [bckarrow_start_x(j) arr_x(j)], [bckarrow_start_y(j) arr_y(j)],'Color',op.bckcolor);
            end
        end

        %bck tx projection
        an = an+1;
        pos = ah(an-1).pos - [0 ad.h 0 0];
        ah(an).pos = pos;
        ah(an).name = 'bck tx projection';
        ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
        pcolor(ah(an).ax, figData.xtproj.tx, figData.xtproj.px, xt_proj_im);
        shading(ah(an).ax, 'flat')
        set(ah(an).ax,'XLim',figData.parameters.bckTrange,'Clim',xt_proj_lim);
        axis(ah(an).ax, 'off');

        colormap(ah(an).ax, xt_proj_map);
    
            
        tc = figData.alignedActivity.tcross(figData.alignedActivity.tcross > figData.parameters.bckTrange(1) & figData.alignedActivity.tcross < figData.parameters.bckTrange(2));
        
        hold(ah(an).ax,'on');
        for k = 1:length(tc)
            plot(ah(an).ax,tc(k) + [0 0], ah(an).ax.YLim, 'w:','LineWidth',  op.stride_marker_linewidth);
        end
    end
        

         switch lower(op.data_channel)
            case {'green', 'g', 'grn'}
                fld = 'fwd_green';
            case {'red', 'r', 'rd'}
                fld = 'fwd_red';
            case {'red_uc', 'r_uc'}
                fld = 'fwd_red_uc';
                
            otherwise
                fld = 'fwd';
        end
    
        % aligned activity to motion
        an = an+1;
        pos = ah(an-1).pos;
        if (~op.dobck || figData.parameters.fwdTrange(1) < figData.parameters.bckTrange(1))
            pos(1) = ad.lx2;
        else
            pos(1) = ad.rx2;
        end
        pos(2) = pos(2) - ad.dh;
        pos(3) = ad.w2;
        ah(an).pos = pos;
        ah(an).name = 'fwd avg';
        ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
        for i = 1:size(figData.alignedActivity.(fld),1)
            plot(ah(an).ax, rad2deg(figData.alignedActivity.phase_axis), figData.alignedActivity.(fld)(i,:), 'Color', op.voi_colors(i,:)); hold(ah(an).ax, 'on');
            [mv,I] = max(figData.alignedActivity.(fld)(i,:));
            text(rad2deg(figData.alignedActivity.phase_axis(I)),mv+.05, num2str(i),"FontSize",po.fontsize,'Color',op.voi_colors(i,:),'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontName',po.font);
        end
        xlabel(ah(an).ax, '\phi','interpreter','tex');
        ylabel(ah(an).ax,'Ravg/minRavg');
        set(ah(an).ax, po.axesopts{:}, 'box','off','XLim',[0 360],'XTick',0:90:360);
        mb = min(figData.alignedActivity.bck,[],'all','omitnan');
        if (isempty(mb))
            mb = 10000;
        end
        minrat = min(min(figData.alignedActivity.(fld),[],'all','omitnan'), mb,'omitnan');
        if (~isfinite(minrat))
            minrat = 0;
        end
        mb = max(figData.alignedActivity.bck,[],'all','omitnan');
        if (isempty(mb))
            mb = 0;
        end
        maxrat = max(max(figData.alignedActivity.(fld),[],'all','omitnan'), mb, 'omitnan');
        maxrat = max(maxrat,minrat + 0.1,'omitnan');

        ah(an).ax.YLim = [minrat maxrat];
        if (~op.dobck || figData.parameters.fwdTrange(1) < figData.parameters.bckTrange(1))
            text(ah(an).ax.XLim(1), ah(an).ax.YLim(2), sprintf('%d forward strides',figData.alignedActivity.nfwd), "FontSize",po.fontsize,'Color','k','HorizontalAlignment','left','VerticalAlignment','bottom', 'FontName',po.font);
        else
            text(ah(an).ax.XLim(2), ah(an).ax.YLim(2), sprintf('%d forward strides',figData.alignedActivity.nfwd), "FontSize",po.fontsize,'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom', 'FontName',po.font);
        end

        if (op.dobck)
            
            switch lower(op.data_channel)
                case {'green', 'g', 'grn'}
                    fld = 'bck_green';
                case {'red', 'r', 'rd'}
                    fld = 'bck_red';
                case {'red_uc', 'r_uc'}
                    fld = 'bck_red_uc';
                otherwise
                    fld = 'bck';
            end
            an = an+1;
            pos = ah(an-1).pos;
            if (figData.parameters.fwdTrange(1) < figData.parameters.bckTrange(1))
                pos(1) = ad.rx2;
            else
                pos(1) = ad.lx2;
            end
            ah(an).pos = pos;
            ah(an).name = 'bck avg';
            ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
            for i = 1:size(figData.alignedActivity.fwd,1)
                plot(ah(an).ax, rad2deg(figData.alignedActivity.phase_axis), figData.alignedActivity.(fld)(i,:), 'Color', op.voi_colors(i,:)); hold(ah(an).ax, 'on');
                [mv,I] = max(figData.alignedActivity.(fld)(i,:));
                text(rad2deg(figData.alignedActivity.phase_axis(I)),mv+.05, num2str(i),"FontSize",po.fontsize,'Color',op.voi_colors(i,:),'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontName',po.font);
            end
            xlabel(ah(an).ax, '\phi','Interpreter','tex');
            %  ylabel(ah(an).ax,'Ravg/minRavg')
            set(ah(an).ax, po.axesopts{:}, 'box','off', 'YLim', ah(an-1).ax.YLim,'XLim',[0 360],'XTick',0:90:360);
            if (figData.parameters.fwdTrange(1) < figData.parameters.bckTrange(1))
                text(ah(an).ax.XLim(2), ah(an).ax.YLim(2), sprintf('%d backward strides',figData.alignedActivity.nbck), "FontSize",po.fontsize,'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom', 'FontName',po.font);
            else
                text(ah(an).ax.XLim(1), ah(an).ax.YLim(2), sprintf('%d backward strides',figData.alignedActivity.nbck), "FontSize",po.fontsize,'Color','k','HorizontalAlignment','left','VerticalAlignment','bottom', 'FontName',po.font);
            end
        end

    end

    %     ah(an).name = 'trajectory';
%     ah(an).pos = [ad.rx4 ad.h0-ad.h ad.w4, ad.h];
%     ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
%     clear h;
%     for j = 1:length(figData.labeledTrajectory.fwd)
%         h(j) = plot(ah(an).ax,  figData.labeledTrajectory.fwd{j}(1,:), figData.labeledTrajectory.fwd{j}(2,:), 'Color', op.fwdcolor); %#ok<*AGROW>   
%         hold(ah(an).ax, 'on');
%         %text(mean(figData.labeledTrajectory.fwd{j}(1,:)), mean(figData.labeledTrajectory.fwd{j}(2,:)),falpha(j));
%     end
%     set(h(1),  'DisplayName', 'Forward');
%     set(h(2:end), 'HandleVisibility', 'off');
% 
%     ah(an).h = h;
% 
%     clear h;
%     for j = 1:length(figData.labeledTrajectory.trn)
%         h(j) = plot(ah(an).ax,  figData.labeledTrajectory.trn{j}(1,:), figData.labeledTrajectory.trn{j}(2,:), 'Color', op.othercolor); 
%     end
%     set(h(1),  'DisplayName', 'Sweep/Pause');
%     set(h(2:end), 'HandleVisibility', 'off');
% 
%     ah(an).h = [ah(an).h h];
% 
%     
%     clear h;
%     try
%         for j = 1:length(figData.labeledTrajectory.bck)
%             h(j) = plot(ah(an).ax,  figData.labeledTrajectory.bck{j}(1,:), figData.labeledTrajectory.bck{j}(2,:), 'Color', op.bckcolor); 
%         end
%         set(h(1),  'DisplayName', 'Backward');
%         set(h(2:end), 'HandleVisibility', 'off');
% 
%         ah(an).h = [ah(an).h h];
%     catch
%     end
% 
% 
% 
%     set(ah(an).ax, 'YTickLabel', [],'XTickLabel', [], 'XTick', [], 'YTick', [], 'box','on');
%     axis(ah(an).ax, 'equal');
% %     set(legend(ah(an).ax),'Location','BestOutside');
%     
%     yl = ah(an).ax.YLim;
%     yl = mean(yl) + diff(yl)*[-.55 .55];
%     xl = ah(an).ax.XLim;
%     xl = mean(xl) + diff(xl)*[-.55 .55];
%     ah(an).ax.XLim = xl;
%     ah(an).ax.YLim = yl;
%     axis(ah(an).ax, 'equal');
% 
%     yl = ah(an).ax.YLim;
%     xl =  ah(an).ax.XLim;
% 
%     labeledBar(ah(an).ax, dot(xl,[0.1 0.9])+[-op.trajectory_scale_bar 0], dot(yl,[0.1 0.9]) + [0 0], sprintf('%d um', op.trajectory_scale_bar), 'below', 'k',{'LineWidth',3}, po.labelOptions);
% 
% 
%     an = an+1;
%    pos = ah(an-1).pos;
