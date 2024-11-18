function tsr = markStridesWithBehaviorLabels(tsr,varargin)
%tsr = tsr.markStridesWithBehaviorLabels('valid_trange',trange);
% trange: dataset-specific time range where larva is uncompressed
% use Hilbert transform on centered vfwd (vfwd-vfwd_smooth) to detect strides
% forward strides start/end at vfwd peaks
% backward strides start/end at vfwd troughs

% default parameters
vmin = 10; % unit: um/s; tracker speed threshold for behavior labeling
valid_trange = tsr.tx([1 end]); % larva uncompressed and sleap-annotated time range
behav_colors = [0.75 1 0.75; 1 0.75 0.75]; % behavior shading in verbose plot: (1,:)=light green=forward, (2,:)=light pink=backward
nphasebins = 36;
verbose = 0;

assignApplicable(varargin);

% detect forward and backward time ranges using fwd_dp and vfwd_smooth
behavlabels = zeros(size(tsr.com.fwd_dp));
behavlabels(tsr.com.fwd_dp>cosd(45) & abs(tsr.com.vfwd_smooth)>=vmin) = 1;
behavlabels(tsr.com.fwd_dp<-cosd(45) & abs(tsr.com.vfwd_smooth)>=vmin) = 2;
validinds = tsr.tx>=valid_trange(1) & tsr.tx<=valid_trange(2);
behavlabels(~validinds) = 0;

% get vfwd-vfwd_smooth and make sure there's no NaN (hilbert() doesn't like
% any NaN)
X = tsr.com.vfwd-tsr.com.vfwd_smooth; % tracker speed // to spot-gut axis
X = fillmissing(X,'linear','EndValues',0); % set ends to 0 and interpolate middles

% detect forward stride by vfwd peaks
sig = hilbert(-X'); % do Hilbert transform on -vfwd so that phase rollover occurs at peaks
phi_fwd = angle(sig);
dphi_fwd = diff(phi_fwd);
dphi_fwd = cat(1,NaN,dphi_fwd);
stridestartinds = find(dphi_fwd<=-2*pi*0.9);
nstrides = length(stridestartinds)-1;
S_fwd = repmat(struct('txinds',[],'t',[],'vfwd',[],'tt',[],'vfwd_tt',[],'behav_label',0),nstrides,1);
validstrides = false(nstrides,1);
for j=1:nstrides
    inds = stridestartinds(j):(stridestartinds(j+1)-1);
    % only include strides entirely within forward time range
    if any(behavlabels(inds)~=1)
        % reset forwards to pauses in this stride for visual (won't disturb backward labels)
        tmp = behavlabels(inds);
        tmp(tmp==1) = 0;
        behavlabels(inds) = tmp;
        continue;
    end
    S_fwd(j).txinds = inds;
    S_fwd(j).vfwd = X(inds);
    t = tsr.tx(inds);
    S_fwd(j).t = t; % original tracker time
    S_fwd(j).tt = linspace(t(1),t(end),nphasebins); % tracker time resampled to get 1-to-1 map to phase bins
    S_fwd(j).vfwd_tt = interp1(S_fwd(j).t,S_fwd(j).vfwd,S_fwd(j).tt);
    S_fwd(j).behav_label = 1;
    validstrides(j) = true;
end
S_fwd = S_fwd(validstrides); % remove empty ones
fprintf('%d forward strides detected\n',sum(validstrides));

% detect backward stride by vfwd troughs
sig = hilbert(X'); % do Hilbert transform on +vfwd so that phase rollover occurs at troughs
phi_bck = angle(sig);
dphi_bck = diff(phi_bck);
dphi_bck = cat(1,NaN,dphi_bck);
stridestartinds = find(dphi_bck<=-2*pi*0.9);
nstrides = length(stridestartinds)-1;
S_bck = repmat(struct('txinds',[],'t',[],'vfwd',[],'tt',[],'vfwd_tt',[],'behav_label',0),nstrides,1);
validstrides = false(nstrides,1);
for j=1:nstrides
    inds = stridestartinds(j):(stridestartinds(j+1)-1);
    % only include strides entirely within backward time range
    if any(behavlabels(inds)~=2)
        % reset backwards to pauses in this stride for visual (won't disturb forward labels)
        tmp = behavlabels(inds);
        tmp(tmp==2) = 0;
        behavlabels(inds) = tmp;
        continue;
    end
    S_bck(j).txinds = inds;
    S_bck(j).vfwd = X(inds);
    t = tsr.tx(inds);
    S_bck(j).t = t; % original tracker time
    S_bck(j).tt = linspace(t(1),t(end),nphasebins); % tracker time resampled to get 1-to-1 map to phase bins
    S_bck(j).vfwd_tt = interp1(S_bck(j).t,S_bck(j).vfwd,S_bck(j).tt);
    S_bck(j).behav_label = 2;
    validstrides(j) = true;
end
S_bck = S_bck(validstrides); % remove empty ones
fprintf('%d backward strides detected\n',sum(validstrides));

% combine and sort by stride start time
S = cat(1,S_fwd,S_bck);
nstrides = length(S);
icross = NaN(nstrides,1);
for j=1:nstrides
    icross(j) = S(j).txinds(1);
end
[~,sortorder] = sort(icross,'ascend');
S = S(sortorder);

% detect crawling bouts, merge too-short strides in each bout, then update tcross and icross
S = tsr.detectAndMergeShortStrides(S);
nstrides = length(S);
tcross = NaN(nstrides,1);
icross = NaN(nstrides,1);
for j=1:nstrides
    tcross(j) = S(j).t(1);
    icross(j) = S(j).txinds(1);
end

% save into tsr.com
tsr.com.tcross = tcross; % tracker time (tsr.tx) of stride starts
tsr.com.icross = icross; % tsr.tx ind of stride starts
tsr.com.strides = S;
tsr.com.phasebins = linspace(0,2*pi,nphasebins);

% plot if verbose
if verbose>0
    phi = NaN(size(X));
    phi(behavlabels==1) = phi_fwd(behavlabels==1);
    phi(behavlabels==2) = phi_bck(behavlabels==2);
    dphi = NaN(size(X));
    dphi(behavlabels==1) = dphi_fwd(behavlabels==1);
    dphi(behavlabels==2) = dphi_bck(behavlabels==2);
    figure;
    ax(1) = subplot(3,1,1); % vfwd
    plot(tsr.tx,X); hold on;
    for j=1:nstrides
        xline(S(j).t([1 end]),'r--');
    end
    hold off;
    title('data');
    ax(2) = subplot(3,1,2); % instantaneous phase
    plot(tsr.tx,phi); hold on;
    for j=1:nstrides
        xline(S(j).t([1 end]),'r--');
    end
    hold off;
    title('\phi'); ylim([-1 1]*1.5*pi); yticks([-1 0 1]*pi); yticklabels({'-\pi','0','\pi'});
    ax(3) = subplot(3,1,3); % diff(instantaneous phase)
    plot(tsr.tx,dphi); hold on;
    for j=1:nstrides
        xline(S(j).t([1 end]),'r--');
    end
    hold off;
    title('\Delta\phi'); ylim([-1 1]*2*pi); yticks([-1 0 1]*2*pi); yticklabels({'-2\pi','0','2\pi'});
    shadedRegions(ax,tsr.tx,behavlabels==1,behav_colors(1,:));
    shadedRegions(ax,tsr.tx,behavlabels==2,behav_colors(2,:));
    linkaxes(ax,'x'); xlim(tsr.tx([1 end]));
end

end