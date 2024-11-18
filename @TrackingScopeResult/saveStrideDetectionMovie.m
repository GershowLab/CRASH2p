function saveStrideDetectionMovie(tsr,varargin)
%tsr.saveStrideDetectionMovie;
% saves movie of behavior image + vfwd with stride and behavior type labels
% to behavior_video/stride_alignment/ if dstdir is not specified

dstdir = [];
tw = 20; % unit: s; vfwd plot window width
trange = tsr.behavior.et([1 end]);
frame_rate = tsr.behavior.framerate_Hz;

assignApplicable(varargin);

if range(trange)>20
    x = input(sprintf('saving this movie will take ~%s; proceed? y/n\n',estetime(range(trange)*3)),'s');
    if strcmpi(x,'n')
        return;
    end
end

if isempty(dstdir)
    dstdir = fullfile(fileparts(tsr.behavior.videoFileName),'stride_alignment');
end

behav_str = {'','forward','backward'};
behav_text_colors = [0 0 0;0 1 0;1 0 0];
behav_shade_colors = min(behav_text_colors+0.75,ones(size(behav_text_colors)));

% get everything from tsr
t_behav = tsr.behavior.et;
X = tsr.com.vfwd-tsr.com.vfwd_smooth;
X_behav = interp1(tsr.tx,X,t_behav);
yl = max(abs(X))*[-1 1];
nstrides = length(tsr.com.strides);
behav_labels_t = zeros(size(t_behav));
for j=1:nstrides
    tedges = tsr.com.strides(j).t([1 end]);
    behav_labels_t(t_behav>=tedges(1) & t_behav<=tedges(2)) = tsr.com.strides(j).behav_label;
end

vw = VideoWriter(fullfile(dstdir,'stride-detection.mp4'),'MPEG-4');
vw.FrameRate = frame_rate;

% first frame
fig = figure('Unit','pixels','Position',[50 100 480 720]);
% 1) behavior image
subplot(3,1,[1 2]);
bvr = VideoReader(tsr.behavior.videoFileName);
im = tsr.getBehaviorImage(bvr,t_behav(1));
h_behav = imagesc(im); axis image off; caxis([0 200]); colormap(gray(256));
k = behav_labels_t(1)+1;
h_label = text(size(im,2)/2,size(im,1)-10,behav_str{k},'Color',behav_text_colors(k,:),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
% 2) vfwd
ax_vfwd = subplot(3,1,3);
% plot vfwd
plot(t_behav,X_behav); ylim(yl*1.1); hold on;
% shade by forward/backward
shadedRegions(ax_vfwd,t_behav,behav_labels_t==1,behav_shade_colors(2,:));
shadedRegions(ax_vfwd,t_behav,behav_labels_t==2,behav_shade_colors(3,:));
% plot strides
for j=1:nstrides
    xline(tsr.com.strides(j).t([1 end]),'k:');
    text(tsr.com.strides(j).t(1),yl(2)*0.9,num2str(j));
end
% TODO fix this: text goes outside axes box
% plot time pointer
h_vfwd = xline(t_behav(1),'r'); hold off;
xlim(t_behav(1)+[-1 1]*tw); xticks(t_behav(1)+[-1 0 1]*tw);
xticklabels({sprintf('-%ds',tw),sprintf('t=%.1fs',t_behav(1)),sprintf('+%ds',tw)});
ylim(yl);
title('vfwd');

validinds = find(t_behav>=trange(1) & t_behav<=trange(2));
open(vw);
for i=validinds
    if ~ishandle(fig)
        return;
    end
    [im,frameind] = tsr.getBehaviorImage(bvr,t_behav(i));
    h_behav.CData = im;
    k = behav_labels_t(i)+1;
    h_label.String = behav_str{k};
    h_label.Color = behav_text_colors(k,:);
    h_vfwd.Value = t_behav(i);
    xlim(ax_vfwd,t_behav(i)+[-1 1]*tw); xticks(ax_vfwd,t_behav(i)+[-1 0 1]*tw);
    xticklabels(ax_vfwd,{sprintf('-%ds',tw),sprintf('t=%.1fs',t_behav(i)),sprintf('+%ds',tw)});
    sgtitle(sprintf('t=%.2fs, behavior frame %d',t_behav(i),frameind));
    drawnow;
    writeVideo(vw,getframe(fig));
end

close(vw);
close(fig);

disp('...done!');

end