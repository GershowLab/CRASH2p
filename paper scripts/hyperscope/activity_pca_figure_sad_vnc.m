function [ah,op] = activity_pca_figure_sad_vnc (sad, opp)
%combines pca and phase sequence
op.cmap_r = inferno(256);
op.cmap_g = viridis(256);%
op.glim = [];%[0 1.5e6];
op.rlim = [];%[0 3e6];

op.fwdcolor = [0 .8 .6];
op.bckcolor = [.8 0 0.2];
op.othercolor = [0 0 0];

op.trajectory_scale_bar = 500; %um
op.fignum = 2;
op.activity_overlap = 0.5;
op.activity_range = [0.5 4];
op.norm_rat_range = [0.5 4];
op.color = 1;
op.activity_linewidth_zoom = 1.5;
op.nrows = 1; %needs revision to work with nrows = 2 or more
op.dovois = false;
op.voi_plot_ratio = false;
op.dobck = [];
op.dofwd = [];
op.restrict_to_pca_inds = true;

if (nargin <= 0)
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

if (isempty(op.dofwd))
    op.dofwd = sad.nfwd > 4;
end
if (isempty(op.dobck))
    op.dobck = sad.nbck > 4;
end

[ad,po] = blank8x10Figure(op.fignum, 'nrows',2*op.nrows + 1,'fontsize',8);

an = 0;
if (op.dofwd)
    f = {'fwdaligned'};
end
if (op.dobck)
    f = {'bckaligned'};
end
if (op.dofwd && op.dobck)
    f = {'fwdaligned', 'bckaligned'};
end

if (op.restrict_to_pca_inds && isfield(sad.fwdaligned, 'pca_xinds'))
    xi = sad.fwdaligned.pca_xinds;
    yi = sad.fwdaligned.pca_yinds;
else
    xi = 1:size(sad.fwdaligned.green,1);
    yi = 1:size(sad.fwdaligned.green,2);
end
    

if(isempty(op.glim))
    op.glim = [0 percentile(abs(sad.(f{1}).green_pca_wave_z(xi,yi)) + mean(sad.(f{1}).green(xi,yi,:),3), .95)];
end

if(isempty(op.rlim))
    op.rlim = [0 percentile(abs(sad.(f{1}).red_pca_wave_z(xi,yi)) + mean(sad.(f{1}).red(xi,yi,:),3), .95)];
end



c = {'green', 'red'};
pleft = [ad.lx3, ad.cx3, ad.rx3];
cm = {op.cmap_g, op.cmap_r};
ttla = {'forward', {'pca amplitude','backward'},'fused'};
ttlp = {'forward', 'backward'};

ncolumns = ceil(numel(sad.fwdaligned.phase_axis)/op.nrows);


dh = ad.dh;
blockdh = op.nrows*ad.dh;

if (op.restrict_to_pca_inds && isfield(sad.fwdaligned, 'pca_xaxis'))
    xaxis = sad.fwdaligned.pca_xaxis;
    yaxis = sad.fwdaligned.pca_yaxis;
else
    xaxis = sad.fwdaligned.xaxis;
    yaxis = sad.fwdaligned.yaxis;
end
if (op.dofwd)
    rm = sad.fwdaligned.red_mask;
else
    rm = sad.bckaligned.red_mask;
end
if (op.dobck)
    rm = rm & sad.bckaligned.red_mask;
end
rm = imclose(rm,ones(7));
rmnan = double(rm); rmnan(~rm) = NaN;

j = op.color;
for i = 1:length(f)
    for k = 1:numel(sad.fwdaligned.phase_axis)
        if (mod(k,ncolumns) == 1)
            pos = [ad.leftmargin ad.h0-ad.h-((i-1)*blockdh + floor(k/ncolumns)*dh)  (1-ad.leftmargin-ad.rightmargin)/ncolumns ad.h];
        end
        an = an + 1;
        ah(an).name = sprintf('%s %s phase %.1f', f{i}, c{j}, rad2deg(sad.(f{i}).phase_axis(k)));
        ah(an).pos = pos;
        ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
        im = sad.(f{i}).(c{j})(:,:,k);
        pcolor(ah(an).ax, sad.(f{i}).yaxis, sad.(f{i}).xaxis, im); shading flat; axis equal; colormap(ah(an).ax, cm{j} );
        ah(an).ax.CLim = op.([c{j}(1) 'lim']);

        axis(ah(an).ax, 'tight'); 
   %     set(ah(an).ax, 'XTick',[],'YTick',[])
        axis(ah(an).ax, 'off')
 %       axis (ah(an).ax, 'off');
        if (i == 1 && mod(k,2) == 1)
            tthndl = title(sprintf(' %.1f',rad2deg(sad.(f{i}).phase_axis(k))),'Interpreter','tex');
            ttpos = dsxy2figxy_marc(tthndl.Parent, tthndl.Extent);
        end
        [x,y] = dsxy2figxy_marc(ah(an).ax(1), sad.(f{i}).yaxis([1 end]), sad.(f{i}).xaxis([1 end]));
        if (k == 1)
            hy = diff(y);
            wx = diff(x);
            annotation("textbox",[x(1)-wx/10,y(1),hy, wx/4], "String", sprintf('\\langle{%s}_{%s}(x,y,\\phi)\\rangle',upper(c{j}(1)), upper(f{i}(1))), 'interpreter','tex','FontSize',po.fontsize, 'EdgeColor', 'none',...
                'FontName', po.font, 'Rotation',90,'HorizontalAlignment','left','VerticalAlignment','middle','FitBoxToText','off');
                     %   ylabel(ah(an).ax, sprintf('\\langle{%s}_{%s}(x,y,\\phi)\\rangle',upper(c{j}(1)), upper(f{i}(1))));

        end
        ah(an).pos = [x(1) y(1) wx hy];
        set (ah(an).ax(1),'Position',ah(an).pos, po.axesopts{:});
        pos = ah(an).pos;
        pos(1) = pos(1) + 1.05*pos(3);
        
        dh = 1.05*ah(an).pos(4);%0.5*(ad.dh+ah(an).pos(4));
        dhbar = dh;
        if (i == 1 && op.nrows > 1) 
            dh = ttpos(2)+ttpos(4) - pos(2);
            dhbar = 0.5*dhbar + 0.5*dh;
        end
        blockdh = op.nrows*dhbar;
    end
end

pos = [ad.leftmargin pos(2)-ad.dh  ad.h ad.h];

for i = 1:(3-2*mod(length(f),2)) %1 if length(f) is 1, 3 if length(f) is 2
    an = an + 1;
    if (i < 3)
        ah(an).name = [f{i} '_' c{j} '_pca_ampl'];
    else
        ah(an).name = ['fused_' c{j} '_pca_ampl'];
    end
    ah(an).pos = pos;
    ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
    if (i < 3)
        im = abs(sad.(f{i}).([c{j} '_pca_z']));
        pcolor(ah(an).ax, yaxis, xaxis, im(xi,yi).*rmnan(xi,yi)); shading flat; axis equal; colormap(ah(an).ax, cm{j} );
        if (op.dovois)
            for k = 1:length(sad.vois)
                B = bwboundaries(imdilate(any(sad.vois{k}.makeMask,3), ones(3)));       
                hold(ah(an).ax, 'on')
                
           %     plot (ah(an).ax,  sad.templatey(B{1}(:,2)), sad.templatex(B{1}(:,1)), 'Color', 'w', 'LineWidth',2);
                plot (ah(an).ax,  sad.fwdaligned.yaxis(B{1}(:,2)), sad.fwdaligned.xaxis(B{1}(:,1)), 'Color', op.voi_colors(k,:), 'LineWidth',1.5,'LineStyle','-');
            end
        end
        ah(an).ax.CLim = op.([c{j}(1) 'lim']);
    else
        im = imfuse(abs(sad.(f{1}).([c{j} '_pca_z']))/ op.([c{j}(1) 'lim'])(2), abs(sad.(f{2}).([c{j} '_pca_z']))/ op.([c{j}(1) 'lim'])(2), 'Scaling','none');
        im(repmat(~rm,[1 1 3])) = 255;
        image(ah(an).ax, yaxis, xaxis, im); shading flat; axis xy;
        if (op.dovois)
            for k = 1:length(sad.vois)
                B = bwboundaries(imdilate(any(sad.vois{k}.makeMask,3), ones(3)));       
                hold(ah(an).ax, 'on')
                
                %plot (ah(an).ax,  sad.templatey(B{1}(:,2)), sad.templatex(B{1}(:,1)), 'Color', 'w', 'LineWidth',2);
                plot (ah(an).ax,  sad.fwdaligned.yaxis(B{1}(:,2)), sad.fwdaligned.xaxis(B{1}(:,1)), 'Color', op.voi_colors(k,:), 'LineStyle', '-' , 'LineWidth',1.5);
%                 texty = mean(sad.templatex(B{1}(:,1)));
%                 text(ah(an).ax, textx, texty,num2str(k),'Color', op.voi_colors(k,:), 'FontSize',po.fontsize+2,'FontWeight','bold','HorizontalAlignment','left');
            end
        end
    end
    axis(ah(an).ax, 'tight'); axis (ah(an).ax, 'off');
    title(ah(an).ax, ttla{i});



    [x,y] = dsxy2figxy_marc(ah(an).ax(1), yaxis([1 end]), xaxis([1 end]));
    ah(an).pos = [x(1) y(1) diff(x) diff(y)];
    set (ah(an).ax(1),'Position',ah(an).pos, po.axesopts{:});

    if (i == 1)
        cbw = (ah(an).pos(1)-pos(1))/3;
        ah(an).pos(1) = ah(an).pos(1) - (ah(an).pos(1)-pos(1))/2;

        ah(an).ax(1).Position = ah(an).pos;
        cbpos = pos;
        cbpos(3) = cbw;
        ah(an).ax(2) = colorbar(ah(an).ax, 'manual', 'Position', cbpos, 'FontSize', po.fontsize, 'FontName', po.font);
    end
    pos = ah(an).pos;
    pos(1) = pos(1) + 1.05*pos(3);
end



for i = 1:length(f)



    an = an + 1;
    ah(an).name = [f{i} '_' c{j} '_pca_phase'];
    ah(an).pos =pos;
    ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
    %im = rad2deg(angle(imblur(sad.(f{i}).([c{j} '_pca_z']),.5)));
    im = rad2deg(angle(sad.(f{i}).([c{j} '_pca_z'])));
    hsvrev = hsv(256);
    hsvrev = hsvrev(end:-1:1,:);
    pcolor(ah(an).ax, yaxis, xaxis, im(xi,yi).*rmnan(xi,yi)); shading flat; axis equal; colormap(ah(an).ax, hsvrev);
    axis(ah(an).ax, 'tight'); axis off;
    title(ah(an).ax, ttlp{i});
    [x,y] = dsxy2figxy_marc(ah(an).ax, yaxis([1 end]), xaxis([1 end]));
    ah(an).pos = [x(1) y(1) diff(x) diff(y)];
    set (ah(an).ax,'Position',ah(an).pos, po.axesopts{:},'CLim',[-180 180]);
    pos = ah(an).pos;
    pos(1) = pos(1) + 1.05*pos(3);

    if (i == length(f))
        cbpos = pos;
        cbpos(3) = cbw;
        ah(an).ax(2) = colorbar(ah(an).ax, 'manual', 'Position', cbpos, 'FontSize', po.fontsize, 'FontName', po.font, 'Ticks',-180:90:180);
    end

    pos = ah(an).pos;
    pos(1) = pos(1) + 1.05*pos(3);
end
an = an+1;
ah(an).name = 'principal componenents';
ah(an).pos = cbpos;
ah(an).pos(1) = ah(an).pos(1)+3*ah(an).pos(3);
ah(an).pos(3) = max(ad.w5, min(ad.w3, 1-ad.rightmargin-ah(an).pos(1)));
ah(an).pos(1) = 1-ad.rightmargin-ah(an).pos(3);
ah(an).pos(2) = ah(an).pos(2) + ah(an).pos(4)/4;
ah(an).pos(4) = 0.75*ah(an).pos(4);
ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
if (op.dofwd)
    hh = plot(ah(an).ax, rad2deg(sad.fwdaligned.phase_axis), sad.fwdaligned.([c{j} '_pca_coeff']), 'Color',op.fwdcolor, 'LineWidth',op.activity_linewidth_zoom);
    hold(ah(an).ax, 'on')
    set(hh(2),'LineStyle','--');
end
if (op.dobck)
    hh = plot(ah(an).ax, rad2deg(sad.bckaligned.phase_axis), sad.bckaligned.([c{j} '_pca_coeff']), 'Color',op.bckcolor, 'LineWidth',op.activity_linewidth_zoom);
    set(hh(2),'LineStyle','--');
%legend('fwd1','fwd2','bck1','bck2','Location','northoutside');
end
xlabel(ah(an).ax, '\phi','Interpreter','tex')
set(ah(an).ax, 'XTick',0:90:360,'XLim',[0 360],'YAxisLocation', 'Right',po.axesopts{:}, 'box','on','FontSize',12)

if (~op.dovois)
    return;
end

an = an+1;
pos = ah(an-2).pos;
if (~op.dobck || isempty(sad.parameters.bckTrange) || (sad.parameters.fwdTrange(1) < sad.parameters.bckTrange(1)))
    pos(1) = ad.lx2;
else
    pos(1) = ad.rx2;
end
pos(2) = pos(2) - ad.dh;
pos(3) = ad.w2;
ah(an).pos = pos;
ah(an).name = 'fwd crawl';

inds = sad.activityPlot.tx >= sad.parameters.fwdTrange(1) & sad.activityPlot.tx <= sad.parameters.fwdTrange(2);
% ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
% ah(an).h = plot(sad.activityPlot.tx(inds),sad.activityPlot.y(:,inds));
% for j = 1:length(ah(an).h)
%     set(ah(an).h(j),'Color',op.voi_colors(j,:));
% end



if (op.voi_plot_ratio)
    handles = stackedActivityPlots(sad.activityPlot.tx(inds),sad.activityPlot.y(:,inds), op.voi_colors, op.activity_overlap, ah(an).pos,op.activity_range);
    set(handles.lines, 'LineWidth',1);
    set(handles.p, 'LineWidth',op.activity_linewidth_zoom);
    set(handles.text, po.labelOptions{:});
    ah(an).ax = handles.ax;
    ah(an).h = rmfield(handles,'ax');
else
    yr =[0 max(sad.activityPlot.(c{op.color}),[],"all")]/1e6;% [min(sad.activityPlot.(c{op.color})(:,inds),[],'all') max(sad.activityPlot.(c{op.color})(:,inds),[],"all")];
    ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
    ah(an).h = plot(sad.activityPlot.tx(inds)-sad.activityPlot.tx(find(inds,1,'first')),sad.activityPlot.(c{op.color})(:,inds)/1e6,'LineWidth',op.activity_linewidth_zoom );
    for j = 1:length(ah(an).h)
        set(ah(an).h(j),'Color',op.voi_colors(j,:));
    end
    xlabel(ah(an).ax, 't(s)');
    ylabel(ah(an).ax, 'count rate (MHz)');
    set(ah(an).ax, 'YLim',yr, po.axesopts{:}, 'box','on','FontSize',12)
  %  handles = stackedActivityPlots(sad.activityPlot.tx(inds),sad.activityPlot.(c{op.color})(:,inds), op.voi_colors, op.activity_overlap, ah(an).pos,yr);
end


an = an+1;

if (isempty(sad.parameters.bckTrange) || sad.parameters.fwdTrange(1) < sad.parameters.bckTrange(1))
    pos(1) = ad.rx2;
else
    pos(1) = ad.lx2;
end
ah(an).pos = pos;
ah(an).name = 'bck crawl';

if (op.dobck && numel(sad.parameters.bckTrange) >= 2)

    inds = sad.activityPlot.tx >= sad.parameters.bckTrange(1) & sad.activityPlot.tx <= sad.parameters.bckTrange(2);
    if (any(inds))
      %  handles = stackedActivityPlots(sad.activityPlot.tx(inds),sad.activityPlot.y(:,inds), op.voi_colors, op.activity_overlap, ah(an).pos,op.activity_range);
      if (op.voi_plot_ratio)
          handles = stackedActivityPlots(sad.activityPlot.tx(inds),sad.activityPlot.y(:,inds), op.voi_colors, op.activity_overlap, ah(an).pos,op.activity_range);
          delete(handles.lines);
          delete(handles.text);

          set(handles.p, 'LineWidth',op.activity_linewidth_zoom);

          ah(an).ax = handles.ax;
          for j = 1:length(ah(an).ax)
              ah(an).ax(j).YLim = ah(an-1).ax(j).YLim;
          end

          ah(an).h = rmfield(handles,{'ax','lines','text'});
      else
%          handles = stackedActivityPlots(sad.activityPlot.tx(inds),sad.activityPlot.(c{op.color})(:,inds), op.voi_colors, op.activity_overlap, ah(an).pos,yr);
           % yr =[0 max(sad.activityPlot.(c{op.color})(:,inds),[],"all")];% [min(sad.activityPlot.(c{op.color})(:,inds),[],'all') max(sad.activityPlot.(c{op.color})(:,inds),[],"all")];
            ah(an).ax = axes('Position',ah(an).pos, po.axesopts{:});
            ah(an).h = plot(sad.activityPlot.tx(inds)-sad.activityPlot.tx(find(inds,1,'first')),sad.activityPlot.(c{op.color})(:,inds)/1e6,'LineWidth',op.activity_linewidth_zoom);
            for j = 1:length(ah(an).h)
                set(ah(an).h(j),'Color',op.voi_colors(j,:));
            end
            xlabel(ah(an).ax, 't(s)');
            set(ah(an).ax, 'YLim', yr,  po.axesopts{:}, 'box','on','FontSize',12, 'YTickLabel',{})
      end
      
    end
end

end

