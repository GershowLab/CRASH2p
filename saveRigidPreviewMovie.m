function saveRigidPreviewMovie(vra,vraopts)

% save preview movie for rigid registration results
% - voxel size = 1x1x2um
% - frame rate ~ 1Hz (but use vra.binFrame to combine multiple frames
% together, so not always exactly 1Hz)
% - blursigma = 1 (use 3D gaussian blur)

% disp('saving rigid registration preview movie...');

vr = vra.vr;

%--------------------------------
% movie playback and saving prep
%--------------------------------

fpath = fullfile(vraopts.dstdir,[vraopts.timestamp '_preview_rigid']);
if isunix % mp4 profile not available on hpc
    fpath = [fpath '.avi'];
    vw = VideoWriter(fpath,'Motion JPEG AVI');
else
    fpath = [fpath '.mp4'];
    vw = VideoWriter(fpath,'MPEG-4');
end
n = round(vr.frame.framerate_Hz); % number of original frames to combine
vw.FrameRate = vr.frame.framerate_Hz/n;
% vw.FrameRate = vr.frame.framerate_Hz;
open(vw);

frameinds = vraopts.alignedInds;
et = vr.frame.edges(frameinds)+vr.frame.frametime/2;
if mod(length(frameinds),n)~=0
    frameinds = [frameinds NaN(1,n-mod(length(frameinds),n))];
    et = [et NaN(1,n-mod(length(et),n))];
end
frameinds = reshape(frameinds,n,[])';
et = mean(reshape(et,n,[]),1,'omitnan'); % only used for behavior and tracker panels

%-----------
% data prep
%-----------

% behavior
behav = vr.behavior;
bvr = VideoReader(behav.videoFileName);

% tracker
nloc = interp1(vr.tsr.tx,lowpass1D(vr.tsr.neuron.loc(1:2,:),median(diff(et))/3)',et)';
vi = interp1(vr.tsr.tx,double(vr.tsr.tracker.tracking),et,'nearest');
vi(isnan(vi)) = 0;
nloc(:,~vi) = NaN;

% fluorescence
vra = vra.createDisplayAssembler([],[1 1 2]);
fa = vra.fadisplay;

%-------------
% figure prep
%-------------

% utility variables
cmap_r = inferno(256);
cmap_g = cmap_r(:,[2 1 3]);
plottitle = fileparts(vr.filename);
maxlen = 67;
if length(plottitle)>maxlen
    plottitle = fullfile(plottitle(1:2),['...' plottitle(end-(maxlen-5):end)]);
end

% set up figure window
set(0,'Units','pixels');
mp = get(0,'MonitorPositions');
u = floor(min(mp(1,3)/20,mp(1,4)/19));
fig = figure('Units','pixels','Position',[mp(1,1)+u/2 mp(1,2)+u/2 19*u 18*u]);

% set up axes
ax.behav = axes(fig,'Units','pixels','Position',[1,11,5,5]*u);
ax.tracker = axes(fig,'Units','pixels','Position',[1,3,5,5]*u);
ax.green.xy = axes(fig,'Units','pixels','Position',[7,11,5,5]*u);
tmp = get(ax.green.xy,'Title');
set(tmp,'Units','pixels');
fluotitlepos = tmp.Position(1:2);
ax.green.xz = axes(fig,'Units','pixels','Position',[7,9,5,2]*u);
ax.green.yz = axes(fig,'Units','pixels','Position',[12,11,2,5]*u);
ax.red.xy = axes(fig,'Units','pixels','Position',[7,3,5,5]*u);
ax.red.xz = axes(fig,'Units','pixels','Position',[7,1,5,2]*u);
ax.red.yz = axes(fig,'Units','pixels','Position',[12,3,2,5]*u);
ax.stats = axes(fig,'Units','pixels','Position',[16,11,2,5]*u);

%-------------
% first frame
%-------------

w = 10; % pixel stats panel bar/whisker size

% get data
% behavior
try
    im_behav = vr.getBehaviorImage(bvr,et(1));
    im_behav = flip(flip(im_behav,1),2);
catch
    im_behav = zeros(bvr.Height,bvr.Width);
end
% fluorescence and voxel stats
fi = frameinds(1,:);
fi(isnan(fi)) = [];
% fi = frameinds(1);
[ct_r,dw_r] = vra.binFrame(fi,fa,[],'r');
rv = imgaussfilt3(ct_r,1)./imgaussfilt3(dw_r,1);
rv(~isfinite(rv)) = NaN;
mmat = imclose(dw_r>0,ones([9 9 9])); rv(~mmat) = NaN;
[ct_g,dw_g] = vra.binFrame(fi,fa,[],'g');
gv = imgaussfilt3(ct_g,1)./imgaussfilt3(dw_g,1);
gv(~isfinite(gv)) = NaN;
mmat = imclose(dw_g>0,ones([9 9 9])); gv(~mmat) = NaN;
sz = size(rv);

% adjust fluorescence panel sizes if image is not square
if sz(1)~=sz(2)
    % adjust xy panel
    if sz(1)>sz(2)
        newwidth = ax.green.xy.Position(3);
        newheight = (sz(2)/sz(1))*newwidth;
        offset = [0 (ax.green.xy.Position(4)-newheight)/2];
        newleft_g = ax.green.xy.Position(1);
        newbottom_g = ax.green.xy.Position(2)+offset(2);
        newleft_r = ax.red.xy.Position(1);
        newbottom_r = ax.red.xy.Position(2)+offset(2);
    else
        newheight = ax.green.xy.Position(4);
        newwidth = (sz(1)/sz(2))*newheight;
        offset = [(ax.green.xy.Position(3)-newwidth)/2 0];
        newbottom_g = ax.green.xy.Position(2);
        newleft_g = ax.green.xy.Position(1)+offset(1);
        newbottom_r = ax.red.xy.Position(2);
        newleft_r = ax.red.xy.Position(1)+offset(1);
    end
    set(ax.green.xy,'Position',[newleft_g newbottom_g newwidth newheight]);
    set(ax.red.xy,'Position',[newleft_r newbottom_r newwidth newheight]);    
    % adjust xz panel
    newwidth = ax.green.xy.Position(3);
    newheight = (2/5)*ax.green.xy.Position(4);
    newleft_g = ax.green.xy.Position(1);
    newbottom_g = ax.green.xy.Position(2)-newheight;
    newleft_r = ax.red.xy.Position(1);
    newbottom_r = ax.red.xy.Position(2)-newheight;
    set(ax.green.xz,'Position',[newleft_g newbottom_g newwidth newheight]);
    set(ax.red.xz,'Position',[newleft_r newbottom_r newwidth newheight]);
    fluotitlepos = fluotitlepos-offset;
    % adjust yz panel
    newheight = ax.green.xy.Position(4);
    newwidth = ax.green.xz.Position(4);
    newleft_g = ax.green.xy.Position(1)+ax.green.xy.Position(3);
    newbottom_g = ax.green.xy.Position(2);
    newleft_r = ax.red.xy.Position(1)+ax.red.xy.Position(3);
    newbottom_r = ax.red.xy.Position(2);
    set(ax.green.yz,'Position',[newleft_g newbottom_g newwidth newheight]);
    set(ax.red.yz,'Position',[newleft_r newbottom_r newwidth newheight]);
end

% get projections (mean)
rI.xy = mean(rv,3,'omitnan');
rI.xz = squeeze(mean(rv,2,'omitnan'));
rI.yz = squeeze(mean(rv,1,'omitnan'));
gI.xy = mean(gv,3,'omitnan');
gI.xz = squeeze(mean(gv,2,'omitnan'));
gI.yz = squeeze(mean(gv,1,'omitnan'));
% use first frame to determine color scale limits
cl_r = prctile(rI.xy(rI.xy>0 & isfinite(rI.xy)),[50 100]);
cl_g = prctile(gI.xy(gI.xy>0 & isfinite(gI.xy)),[50 100]);
yl_r = prctile(rI.xy(rI.xy>cl_r(1) & isfinite(rI.xy)),[99 100])';
yl_g = prctile(gI.xy(gI.xy>cl_g(1) & isfinite(gI.xy)),[99 100])';
yl = [min([yl_r yl_g],[],'all','omitnan') max([yl_r yl_g],[],'all','omitnan')];
yl(1) = yl(1)-0.1*range(yl);
yl(2) = yl(2)+0.1*range(yl);
% get top 1% non-zero xy pixel stats
rmin = prctile(rI.xy(rI.xy>cl_r(1) & isfinite(rI.xy)),99);
rstats = prctile(rI.xy(rI.xy>=rmin & isfinite(rI.xy)),[0 25 75 100]);
gmin = prctile(gI.xy(gI.xy>cl_g(1) & isfinite(gI.xy)),99);
gstats = prctile(gI.xy(gI.xy>=gmin & isfinite(gI.xy)),[0 25 75 100]);

% fill axes
% behavior
h.behav = imagesc(ax.behav,im_behav);
set(ax.behav,'XTick',[],'YTick',[],'ColorMap',gray,'CLim',[0 200]);
str = ['behavior camera (AVI time: '...
    datestr(bvr.CurrentTime/(24*60*60),'MM:SS.FFF') ')'];
ax.behav.Title.String = str;
% tracker
plot(ax.tracker,nloc(1,:),nloc(2,:),'k'); hold(ax.tracker,'on');
h.tracker = plot(ax.tracker,nloc(1,1),nloc(2,1),'ro','LineWidth',1.5);
hold(ax.tracker,'off');
axis(ax.tracker,'equal','padded');
title(ax.tracker,'tracked neuron trajectory');
% fluorescence
% 1) green
% 1.a) xy projection
h.green.xy = imagesc(ax.green.xy,gI.xy');
set(ax.green.xy,'XTick',[],'YTick',[],'YDir','normal','ColorMap',cmap_g,'CLim',cl_g);
title(ax.green.xy,'green fluorescence (photon rates)','Units','pixels','Position',fluotitlepos);
% 1.b) xz projection
h.green.xz = imagesc(ax.green.xz,gI.xz');
set(ax.green.xz,'XTick',[],'YTick',[],'YDir','reverse','ColorMap',cmap_g,'CLim',cl_g);
% 1.c) yz projection
h.green.yz = imagesc(ax.green.yz,gI.yz);
set(ax.green.yz,'XTick',[],'YTick',[],'YDir','normal','ColorMap',cmap_g,'CLim',cl_g);
pos = [14.25 11 0.5 5]*u;
colorbar(ax.green.xy,'Units','pixels','Position',pos);
% 2) red
% 2.a) xy projection
h.red.xy = imagesc(ax.red.xy,rI.xy');
set(ax.red.xy,'XTick',[],'YTick',[],'YDir','normal','ColorMap',cmap_r,'CLim',cl_r);
title(ax.red.xy,'red fluorescence (photon rates)','Units','pixels','Position',fluotitlepos);
% 2.b) xz projection
h.red.xz = imagesc(ax.red.xz,rI.xz');
set(ax.red.xz,'XTick',[],'YTick',[],'YDir','reverse','ColorMap',cmap_r,'CLim',cl_r);
% 2.c) yz projection
h.red.yz = imagesc(ax.red.yz,rI.yz);
set(ax.red.yz,'XTick',[],'YTick',[],'YDir','normal','ColorMap',cmap_r,'CLim',cl_r);
pos = [14.25 3 0.5 5]*u;
colorbar(ax.red.xy,'Units','pixels','Position',pos);
% pixel stats
hold(ax.stats,'on');
h.stats.g1 = plot(ax.stats,[1 1],gstats([1 4]),'Color',cmap_g(127,:),'LineWidth',1,'Marker','_','MarkerSize',w);
h.stats.g2 = plot(ax.stats,[1 1],gstats([2 3]),'Color',cmap_g(127,:),'LineWidth',w);
h.stats.r1 = plot(ax.stats,[1 1]+1,rstats([1 4]),'Color',cmap_r(127,:),'LineWidth',1,'Marker','_','MarkerSize',w);
h.stats.r2 = plot(ax.stats,[1 1]+1,rstats([2 3]),'Color',cmap_r(127,:),'LineWidth',w);
hold(ax.stats,'off');
set(ax.stats,'XLim',[0 3],'XTick',[],'YLim',yl,'YAxisLocation','right','Box','on');
title(ax.stats,'xy pixel stats');
xlabel(ax.stats,{'box = 25-75%','of top 1% pixel'},'Interpreter','none');
% metadata
str = {plottitle;...
    sprintf('time: %.2f-%.2fs, frames %d-%d',...
        vr.frame.edges(fi(1)),vr.frame.edges(fi(end)),fi(1),fi(end))};
% str = {plottitle;...
%         sprintf('time: %.2fs, frame %d',et(1),fi)};
h.plottitle = text(ax.behav,'String',str,...
    'Units','pixels','Position',[0.5 6.25]*u,...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold','Interpreter','none');

%----------------------
% play and save movies
%----------------------

disp('saving rigid preview movie...');

for i=1:size(frameinds,1)
% for i=1:length(frameinds)
    % 1) get data
    % behavior
    try
        im_behav = vr.getBehaviorImage(bvr,et(i));
        im_behav = flip(flip(im_behav,1),2);
    catch
        im_behav = zeros(bvr.Height,bvr.Width);
    end
    % fluorescence
    fi = frameinds(i,:);
    fi(isnan(fi)) = [];
%     fi = frameinds(i);
    [ct_r,dw_r] = vra.binFrame(fi,fa,[],'r');
    rv = imgaussfilt3(ct_r,1)./imgaussfilt3(dw_r,1);
    rv(~isfinite(rv)) = NaN;
    mmat = imclose(dw_r>0,ones([9 9 9])); rv(~mmat) = NaN;
    [ct_g,dw_g] = vra.binFrame(fi,fa,[],'g');
    gv = imgaussfilt3(ct_g,1)./imgaussfilt3(dw_g,1);
    gv(~isfinite(gv)) = NaN;
    mmat = imclose(dw_g>0,ones([9 9 9])); gv(~mmat) = NaN;
    % get projections (mean)
    rI.xy = mean(rv,3,'omitnan');
    rI.xz = squeeze(mean(rv,2,'omitnan'));
    rI.yz = squeeze(mean(rv,1,'omitnan'));
    gI.xy = mean(gv,3,'omitnan');
    gI.xz = squeeze(mean(gv,2,'omitnan'));
    gI.yz = squeeze(mean(gv,1,'omitnan'));
    % get top 1% non-zero xy pixel stats
    rmin = prctile(rI.xy(rI.xy>cl_r(1) & isfinite(rI.xy)),99);
    rstats = prctile(rI.xy(rI.xy>=rmin & isfinite(rI.xy)),[0 25 75 100]);
    gmin = prctile(gI.xy(gI.xy>cl_g(1) & isfinite(gI.xy)),99);
    gstats = prctile(gI.xy(gI.xy>=gmin & isfinite(gI.xy)),[0 25 75 100]);
    % metadata
    str = {plottitle;...
        sprintf('time: %.2f-%.2fs, frames %d-%d',...
        vr.frame.edges(fi(1)),vr.frame.edges(fi(end)),fi(1),fi(end))};
%     str = {plottitle;...
%         sprintf('time: %.2fs, frame %d',et(i),fi)};
    str_behav = ['behavior camera (AVI time: '...
        datestr(bvr.CurrentTime/(24*60*60),'MM:SS.FFF') ')'];
    
    % 2) update axes
    % behavior
    h.behav.CData = im_behav;
    ax.behav.Title.String = str_behav;
    % tracker
    h.tracker.XData = nloc(1,i);
    h.tracker.YData = nloc(2,i);
    % fluorescence
    h.green.xy.CData = gI.xy';
    h.green.xz.CData = gI.xz';
    h.green.yz.CData = gI.yz;
    h.red.xy.CData = rI.xy';
    h.red.xz.CData = rI.xz';
    h.red.yz.CData = rI.yz;
    h.stats.g1.YData = gstats([1 4]);
    h.stats.g2.YData = gstats([2 3]);
    h.stats.r1.YData = rstats([1 4]);
    h.stats.r2.YData = rstats([2 3]);
    % metadata
    h.plottitle.String = str;
    
    drawnow;
    
    writeVideo(vw,getframe(fig));
end

close(vw);

disp('...done!');

end