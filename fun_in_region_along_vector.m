function [proj,px,vpar,vper,boxregion, dpar] = fun_in_region_along_vector (xaxis, yaxis, data, pt1, pt2, dist, fun)
%function fun_in_region_along_vector (xaxis, yaxis, data, pt1, pt2, dist, fun)
%data is ordered x then y then t
%pt1, pt2 

existsAndDefault('fun', @(x) mean(x,'all','omitnan'));

[xx,yy] = ndgrid(xaxis, yaxis);

vpar = (pt2-pt1); vpar = vpar./sqrt(sum(vpar.^2));
lmax = dot(pt2-pt1, vpar);
vper = [-vpar(2) vpar(1)];

xx = xx - pt1(1);
yy = yy - pt1(2);

dpar = vpar(1)*xx + vpar(2)*yy;
dper = vper(1)*xx + vper(2)*yy;

valid = (dpar >= 0 & dpar <= lmax) & abs(dper) <= dist;

dx = median(diff(xaxis));
nedges = ceil(lmax/dx + 1);
peredges = linspace(-dx/2,lmax+dx/2,nedges);

proj = zeros(nedges-1, size(data,3));

[~,~,bin] = histcounts(dpar(valid), peredges);
px = 0.5*(peredges(1:end-1)+peredges(2:end));

boxregion = [pt1(:)-vper(:)*dist pt1(:)+vper(:)*dist pt2(:)+vper(:)*dist pt2(:)-vper(:)*dist  pt1(:)-vper(:)*dist];

if (isempty(data))
    return;
end

for j = 1:size(data,3)
    f = data(:,:,j);
%     clf;
%     pcolor(yaxis, xaxis, f); shading flat; axis equal; colormap gray; 
%     ax2 = axes('position', get(gca,'Position'),'color','none');
%     f2 = NaN(size(f)); f2(valid) = f(valid);
%     pcolor(ax2,yaxis, xaxis, f2); shading flat; axis equal; colormap(ax2,'hot');
%     pause(0.1);

    proj(:,j) = accumarray(bin, f(valid), [nedges-1 1], fun);
end
