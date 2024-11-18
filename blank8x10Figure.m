function [axesdims, po] = blank8x10Figure (fignum, varargin)
%function [axesdims, po] = blank8x10Figure (fignum, varargin)


hspace = .6/11;
nrows = 4;
leftmargin = .5/8;
rightmargin = .5/8;
topmargin = 1/11;
bottommargin = 1/11;
fontsize = 7;
wspace_scale = 1;
varargin = assignApplicable(varargin);

axesdims.leftmargin = leftmargin;
axesdims.rightmargin = rightmargin;
axesdims.hspace = hspace;

allaxeswidth = 1 - leftmargin - rightmargin;
if (existsAndDefault('fignum', []))
    f = figure(fignum);
else
    f = figure();
end
clf(f);
axesdims.fignum = f;
ss = get(0, 'ScreenSize');

screenwidth = ss(3); screenheight = ss(4);
figratio = 8/10;
r = 0.8;
figheight = screenheight*r;
figwidth = figheight*figratio;
figpos = round([(screenwidth - figwidth)/2, (screenheight-figheight)/2, figwidth, figheight]);
    
set(f, 'Position', figpos);
get(f, 'Position');  %calling get f,position here prevents a weird bug where text is incorrectly spaced later

set(f, 'PaperType', 'usletter', 'PaperPosition', [0.25 0.25 8 10], 'PaperOrientation', 'portrait', 'color', 'w', 'inverthardcopy', 'off');


axesdims.h0 = 1-topmargin;
allaxesheight = axesdims.h0-bottommargin;

wspace2 = wspace_scale*0.08 * allaxeswidth;
wspace3 = wspace_scale*0.08 * allaxeswidth;
wspace4 = wspace_scale*0.08 * allaxeswidth;
wspace5 = wspace_scale*0.05 * allaxeswidth;

axesdims.wspace = wspace2;
axesdims.wspace5 = wspace5;

axesdims.h = (allaxesheight - (nrows-1)*hspace)/nrows;
axesdims.dh = axesdims.h + hspace;

axesdims.w1 = allaxeswidth;
axesdims.w2 = min((allaxeswidth - wspace2)/2);%, h*1.61803399);
axesdims.w3 = (allaxeswidth - 2*wspace3)/3;
axesdims.w4 = (allaxeswidth - 3*wspace4)/4;
axesdims.w5 = (allaxeswidth - 4*wspace5)/5;

centerx = 0.5*(1+leftmargin -rightmargin);

axesdims.lx1 = centerx - axesdims.w1/2;
axesdims.lx2 = centerx - wspace2/2 - axesdims.w2;
axesdims.rx2 = centerx + (wspace2)/2;

axesdims.lx3 = centerx - wspace3 - 3*axesdims.w3/2;
axesdims.cx3 = centerx - axesdims.w3/2;
axesdims.rx3 = centerx + axesdims.w3/2 + wspace3;


axesdims.lx4 = centerx - 3*wspace4/2 - 2*axesdims.w4;
axesdims.clx4 = centerx - wspace4/2 - axesdims.w4;
axesdims.crx4 = centerx + (wspace4)/2;
axesdims.rx4 = centerx + 3*(wspace4)/2+axesdims.w4;

axesdims.lx5 = centerx - 2*wspace5 - 2.5*axesdims.w5;
axesdims.clx5 = centerx - wspace5 - 1.5*axesdims.w5;
axesdims.cx5 = centerx - axesdims.w5/2;
axesdims.crx5 = centerx + axesdims.w5/2 + (wspace5);
axesdims.rx5 = centerx + 1.5*axesdims.w5 + 2*(wspace5);

if nargout > 1
    po.lineWidth = 1;
    po.font = 'Arial';
    po.fontsize = fontsize;
    po.bigfontsize = 2*fontsize;
    po.color = 'k';
    po.axesopts = {'FontName', po.font, 'FontSize', po.fontsize, 'LineWidth', po.lineWidth/2, 'box', 'off'};
    po.plotOptions = {'LineWidth', po.lineWidth};
    po.labelOptions = {'Interpreter', 'Tex', 'FontSize', po.fontsize};
    po.redcolor = [1 0 0];
    po.greencolor = [0 1 0.25];
    po.bluecolor = [0 0 1];

end