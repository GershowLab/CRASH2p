function h = shadedRegions(ax, xd, logicalState, color)
%function h = shadedRegions(ax, xd, logicalState, color)
%undershades region of ax indicated by logical state

if (length(ax) > 1)
    for j = 1:length(ax)
        h{j} = shadedRegions(ax(j), xd, logicalState, color);
    end
    return;
end

if (~any(logicalState))
    h = [];
    return;
end

xl = ax.XLim;
yl = ax.YLim;
logicalState(~isfinite(logicalState)) = 0;
maxsize = min(length(xd), length(logicalState));
xd = xd(1:maxsize);
logicalState = logicalState(1:maxsize);
if (~any(logicalState))
    h = [];
    return;
end
ls = xd((diff(logicalState) > 0));
le = xd((diff(logicalState) < 0));
% if (isempty(ls) && isempty(ld))
%     if any(logicalState)
%         ls = xl(1); le = xl(2);
%     else
hc = ax.Children;
ls = ls(ls > xl(1) & ls < xl(2));
le = le(le > xl(1) & le < xl(2));


if (interp1(xd, double(logicalState), max(min(xd), xl(1)), 'nearest', 'extrap'))
    ls = [xl(1) ls];
end
if (interp1(xd, double(logicalState), min(max(xd), xl(2)), 'nearest', 'extrap'))
    le = [le xl(2)];
end
if (isempty(ls))
    h = [];
    return;
end

np = ax.NextPlot;
ax.NextPlot = 'add';
for j = 1:length(ls)
    h(j) = patch([ls(j) ls(j) le(j) le(j)], yl([1 2 2 1]), color, 'FaceColor', color', 'EdgeColor', 'none', 'Parent', ax);
end
ax.Children = [hc' h];
ax.NextPlot = np;