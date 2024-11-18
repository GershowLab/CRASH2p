function handles = stackedActivityPlots(xdata,ydata, colors, overlap, pos, yrange)

if (size(xdata,1)~=1)
    xdata = xdata';
end
if (size(ydata,2) ~= size(xdata,2))
    ydata = ydata';
end

ntrace = size(ydata,1);
if ~existsAndDefault('pos') || (isempty(pos))
    axtemp = axes();
    axp = get(axtemp,'pos');
    delete(axtemp);
else
    axp = pos;
end
h = axp(4)/((1-overlap)*(ntrace-1) + 1);
dh = (1-overlap)*h;

ymax = min(max(ydata,[],'all','omitnan'),100);
ymin = min(ydata,[],'all','omitnan');
xl = [min(xdata) max(xdata)];
existsAndDefault('yrange', [ymin ymax]);
yl = yrange;

for i = 1:ntrace

    pos = [axp(1) axp(2)+dh*(i-1) axp(3) h];
    
    ax(i) = axes('Position',pos, 'Color','none','box','off'); %#ok<*AGROW> %box on is for debugging overlap - 
    p(i) = plot(ax(i), xdata, ydata(i,:), 'DisplayName', sprintf('voi %x',i));
    p(i).Color = colors(i,:);%[colors(i,:) 0.2];
    axis(ax(i),'off');
    ax(i).Color = 'none';
    
    xlim(ax(i), xl)
    ylim(ax(i), yl);

end

handles.ax = ax;
handles.p = p;


dx = 10^(floor(log10(diff(xl)/10))); %largest power of 10 no greater than 10% of the axis length


dy = max(1,10^(floor(log10(diff(yl)/4))));
[handles.lines,handles.text] = labeledBar(ax(end), xl(2) + dx + [0 dx], yl(2)+[-dy -dy], [num2str(dx) 's'], 'below', 'k', {'LineWidth',2}, {'FontSize',12} );

[lll,ttt] =  labeledBar(ax(end), xl(2) + dx + [0 0], yl(2)+[-dy 0], [num2str(dy) 'x'], 'left', 'k', {'LineWidth',2}, {'FontSize',12} ); %may not be x, but keeping for consistency
handles.lines = [handles.lines lll];
handles.text = [handles.text ttt];
%[handles.lines(2),handles.text(2)] = labeledBar(ax(end), xl(2) + dx + [0 0], yl(2)+[-1 0], ['1x'], 'left', 'k', {'LineWidth',2}, {'FontSize',12} );



