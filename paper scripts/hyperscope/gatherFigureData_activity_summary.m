function figData = gatherFigureData_activity_summary (data, parameters)



figData.parameters = parameters;
vs = parameters.vs;
traceinds = parameters.traceinds;

if (~isfield(parameters, 'domovie'))
    parameters.domovie = false;
end


alpha = 0.1;
ma = motion_alignmentHS;

tsr = data.tsr;

%tsr = tsr.findCounterMovements();%'minPeakSep',1, 'vel_sigma', 0.05, 'lp_sigma', 4);

%%trajectory, broken down into forward, backward, and turn

for i = 1:size(parameters.tfwd, 1)
    tidx = ma.get_time_index_tsr(tsr, parameters.tfwd(i,1), parameters.tfwd(i,2));
    if (~isempty(tidx(1):tidx(2)))
        figData.labeledTrajectory.fwd{i} =  tsr.neuron.loc(1:2, tidx(1):tidx(2));
    end
%     figData.labeledTrajectory.fwd_perp{i} =  [-1 1].*diff(tsr.neuron.loc(2:1,[tidx(1) tidx(2)]),[],2);
%     figData.labeledTrajectory.fwd_perp{i} = figData.labeledTrajectory.fwd_perp{i}/sqrt(sum(figData.labeledTrajectory.fwd_perp{i}.^2));
end

for i = 1:size(parameters.ttrn, 1)
    tidx = ma.get_time_index_tsr(tsr, parameters.ttrn(i,1), parameters.ttrn(i,2));
    if (~isempty(tidx(1):tidx(2)))
        figData.labeledTrajectory.trn{i} =  tsr.neuron.loc(1:2, tidx(1):tidx(2));
    end
end

for i = 1:size(parameters.tbck, 1)
    tidx = ma.get_time_index_tsr(tsr, parameters.tbck(i,1), parameters.tbck(i,2));
    if (~isempty(tidx(1):tidx(2)))
        figData.labeledTrajectory.bck{i} =  tsr.neuron.loc(1:2, tidx(1):tidx(2));
    end
end


activity = data.activity;
tfirst_act = min([parameters.tbck; parameters.tfwd; parameters.ttrn],[],'all');
tlast_act = max([parameters.tbck; parameters.tfwd; parameters.ttrn],[],'all');

% tfirst_act = activity.txfine(1); %1994; truncate data a bit
% tlast_act = activity.txfine(end); %2320;


t_idx_act = ma.get_time_index_act(activity, tfirst_act, tlast_act);
ti_act = t_idx_act(1); tf_act = t_idx_act(2);
figData.tr_act = activity.txfine(ti_act:tf_act);

t_idx_tsr = ma.get_time_index_tsr(tsr, activity.txfine(ti_act), activity.txfine(tf_act));
ti_tsr = t_idx_tsr(1); tf_tsr = t_idx_tsr(2);
figData.tr_tsr = tsr.tx(ti_tsr:tf_tsr);

existsAndDefault('traceinds', 1:length(activity.voisets(vs).vois));


figData.timing.tfirst_act = tfirst_act;
figData.timing.tlast_act = tlast_act;

%%activity vs. time in vois
figData.activityPlot.tx = figData.tr_act;
y = (lowpass1D(activity.voisets(vs).lpr(traceinds,:)',2))'./activity.voisets(vs).r0(traceinds); 
figData.activityPlot.y = y(:,ti_act:tf_act);
figData.vois = activity.voisets(vs).vois(traceinds);

try
    gc = (lowpass1D(activity.voisets(vs).gcfine(traceinds,:)',2))';
    gd = (lowpass1D(activity.voisets(vs).gdfine_ic(traceinds,:)',2))';
catch
    gc = (lowpass1D(activity.voisets(vs).gcfine_ic(traceinds,:)',2))';
    gd = (lowpass1D(activity.voisets(vs).gdfine(traceinds,:)',2))';
end

figData.activityPlot.green = gc(:,ti_act:tf_act)./gd(:,ti_act:tf_act);
for k = 1:size(figData.activityPlot.green, 1)
    figData.activityPlot.green0(k) = meanAndStdAssymOutlier(figData.activityPlot.green(k,:));
end

try
    rc = (lowpass1D(activity.voisets(vs).rcfine(traceinds,:)',2))';
    rd = (lowpass1D(activity.voisets(vs).rdfine_ic(traceinds,:)',2))';
    rd_uc = (lowpass1D(activity.voisets(vs).rdfine(traceinds,:)',2))';
    rc_uc = rc;
    
catch
    rc = (lowpass1D(activity.voisets(vs).rcfine_ic(traceinds,:)',2))';
    rd = (lowpass1D(activity.voisets(vs).rdfine(traceinds,:)',2))';
    rc_uc = (lowpass1D(activity.voisets(vs).rcfine(traceinds,:)',2))';
    rd_uc = rd;
end

figData.activityPlot.red = rc(:,ti_act:tf_act)./rd(:,ti_act:tf_act);
figData.activityPlot.red_uc = rc_uc(:,ti_act:tf_act)./rd_uc(:,ti_act:tf_act);

for k = 1:size(figData.activityPlot.red, 1)
    figData.activityPlot.red0(k) = meanAndStdAssymOutlier(figData.activityPlot.red(k,:));
    figData.activityPlot.red_uc0(k) = meanAndStdAssymOutlier(figData.activityPlot.red_uc(k,:));

end

figData.activityPlot.speed = interp1(tsr.tx, tsr.com.vfwd, figData.activityPlot.tx);

figData.templateim = data.final_template.F;
figData.templatex = data.final_template.u;
figData.templatey = data.final_template.v;
figData.templatez = data.final_template.w;
try
    figData.alignedActivity.phase_axis = activity.phase_axis;
catch
    figData.alignedActivity.phase_axis = deg2rad(0:10:360);
end
figData.alignedActivity.tcross = tsr.tx(tsr.com.counterInd);
[y,ally] = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).lpr(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
y = lowpass1D(y',1,'padType','circular')';
activity_denom = min(y,[],2);
figData.alignedActivity.fwd = y./activity_denom;
figData.alignedActivity.nfwd = size(ally,2);
figData.alignedActivity.allfwd = ally./activity_denom;

[~,figData.alignedActivity.allfwd_speed] = tsr.alignDataToMotionPhase(tsr.tx, tsr.com.vfwd, figData.alignedActivity.phase_axis, parameters.tfwd);


if (isfield(activity.voisets(vs), 'gdfine_ic'))
    gc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gcfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    gd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gdfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
else
    gc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gcfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    gd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gdfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
end
figData.alignedActivity.fwd_green = lowpass1D((gc./gd)',1,'padtype', 'circular')';

green_denom = min(figData.alignedActivity.fwd_green,[],2);

if (isfield(activity.voisets(vs), 'rdfine_ic'))
    rc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rcfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    rd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rdfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    rd_uc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rdfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    rc_uc = rc;
else
    rc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rcfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    rd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rdfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    rc_uc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rcfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tfwd);
    rd_uc = rd;
end
figData.alignedActivity.fwd_red = lowpass1D((rc./rd)',1,'padtype', 'circular')';
figData.alignedActivity.fwd_red_uc = lowpass1D((rc_uc./rd_uc)',1,'padtype', 'circular')';

red_denom = min(figData.alignedActivity.fwd_red,[],2);
red_denom_uc = min(figData.alignedActivity.fwd_red_uc,[],2);


if (~isempty(parameters.tbck))
    [y,ally] = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).lpr(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck, 'ignore_valid',true);
    y = lowpass1D(y',1,'padType','circular')';
    figData.alignedActivity.bck = y./activity_denom;
    figData.alignedActivity.nbck = size(ally,2);
    figData.alignedActivity.allbck = ally./activity_denom;
    
    [~,figData.alignedActivity.allbck_speed] = tsr.alignDataToMotionPhase(tsr.tx, tsr.com.vfwd, figData.alignedActivity.phase_axis, parameters.tbck, 'ignore_valid',true);

    if (isfield(activity.voisets(vs), 'gdfine_ic'))
        gc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gcfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        gd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gdfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
    else
        gc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gcfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        gd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).gdfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
    end
    figData.alignedActivity.bck_green = lowpass1D((gc./gd)',1,'padtype', 'circular')';
    %green_denom = min(green_denom,
    %min(figData.alignedActivity.bck_green,[],2)); %for consistency with
    %activity_denom that doesn't take into account back
    figData.alignedActivity.bck_green =  figData.alignedActivity.bck_green./green_denom;
    
    if (isfield(activity.voisets(vs), 'rdfine_ic'))
        rc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rcfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        rd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rdfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        rd_uc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rdfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        rc_uc = rc;
    else
        rc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rcfine_ic(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        rd = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rdfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        rc_uc = tsr.alignDataToMotionPhase(activity.txfine, activity.voisets(vs).rcfine(traceinds,:), figData.alignedActivity.phase_axis, parameters.tbck);
        rd_uc = rd;
    end
    
    figData.alignedActivity.bck_red = lowpass1D((rc./rd)',1,'padtype', 'circular')';
   % red_denom = min(red_denom, min(figData.alignedActivity.bck_red,[],2));%for consistency with
    %activity_denom that doesn't take into account back
    figData.alignedActivity.bck_red = figData.alignedActivity.bck_red./red_denom;
    
    figData.alignedActivity.bck_red_uc = lowpass1D((rc_uc./rd_uc)',1,'padtype', 'circular')';
   % red_denom_uc = min(red_denom_uc, min(figData.alignedActivity.bck_red_uc,[],2));
    figData.alignedActivity.bck_red_uc = figData.alignedActivity.bck_red_uc./red_denom_uc;
    
    
    
else
    figData.alignedActivity.bck = NaN(size(figData.alignedActivity.fwd));
    figData.alignedActivity.bck_green = figData.alignedActivity.bck;
    figData.alignedActivity.bck_red = figData.alignedActivity.bck;
    figData.alignedActivity.nbck = 0;
    figData.alignedActivity.allbck = NaN(size(figData.alignedActivity.allfwd).*[1 0 1]);
end

figData.alignedActivity.fwd_green =figData.alignedActivity.fwd_green./green_denom;
figData.alignedActivity.fwd_red = figData.alignedActivity.fwd_red./red_denom;
figData.alignedActivity.fwd_red_uc = figData.alignedActivity.fwd_red_uc./red_denom_uc;


if (~isfield(parameters, 'proj_pt1') || isempty(parameters.proj_pt1))
    bw = any(logical(Voi.makeLabeledMask(data.voisets(vs).vois)),3);
    [xx,yy] = ndgrid(data.final_template.u, data.final_template.v);
    x = xx(bw);
    y = yy(bw);
    p = polyfit(x,y,1);
    cp = [mean(x) polyval(p,mean(x))];
    vpar = [1 p(1)]/sqrt(1 + p(1)^2); %runs from low x to high x
    vper = vpar([2 1]).*[-1 1];
    [~,I] = min(vpar(1).*x + vpar(2).*y);    
    parameters.proj_pt1 = [x(I) y(I)] + vper*(dot(vper, cp-[x(I) y(I)]));
    [~,I] = max(vpar(1).*x + vpar(2).*y);
    parameters.proj_pt2 = [x(I) y(I)] + vper*(dot(vper, cp-[x(I) y(I)]));
    parameters.proj_dist = max(abs(sum(vper.*[x-cp(1) y-cp(2)],2)));
end


sigma = .5;
[gc,px, vpar, vper,boxregion] = fun_in_region_along_vector(data.final_template.u, data.final_template.v, imblur(data.projections.gc,sigma), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, @sum);
gd = fun_in_region_along_vector(data.final_template.u, data.final_template.v, imblur(data.projections.gd,sigma), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, @sum);
rc = fun_in_region_along_vector(data.final_template.u, data.final_template.v, imblur(data.projections.rc,sigma), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, @sum);
rd = fun_in_region_along_vector(data.final_template.u, data.final_template.v, imblur(data.projections.rd,sigma), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, @sum);

try
    rc_uc = fun_in_region_along_vector(data.final_template.u, data.final_template.v, imblur(data.projections.rc_uc,sigma), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, @sum);
    rd_uc = fun_in_region_along_vector(data.final_template.u, data.final_template.v, imblur(data.projections.rd_uc,sigma), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, @sum);
    figData.xtproj.rrate_uc = rc_uc./rd_uc;
catch
end

figData.xtproj.tx = data.projections.tx;
figData.xtproj.px = px;
figData.xtproj.gcount = gc;
figData.xtproj.rcount = rc;
figData.xtproj.grate = gc./gd;
figData.xtproj.rrate = rc./rd;

figData.xtproj.ratio = (gc.*rd)./(rc.*gd); %included dwell time with intensity corrections also
figData.xtproj.boxregion = boxregion;
figData.xtproj.zrange = data.projections.zrange;
figData.xtproj.pt1 = parameters.proj_pt1;
figData.xtproj.pt2 = parameters.proj_pt2;
figData.xtproj.dist = parameters.proj_dist;
figData.xtproj.vpar = vpar;
figData.xtproj.vper = vper;

if (parameters.domovie && isfield(tsr.behavior, 'at_microns2px'))
    [figData.fwdMovie.imstack,figData.fwdMovie.frametime] = tsr.getBehaviorSequence(parameters.fwdTrange);
    [figData.bckMovie.imstack,figData.bckMovie.frametime] = tsr.getBehaviorSequence(parameters.bckTrange);
end

sad = stride_analysis(data,parameters);

%intersect(fieldnames(figData), fieldnames(sad))

fn = fieldnames(sad);
for j = 1:length(fn)
    figData.(fn{j}) = sad.(fn{j});
end

figData.parameters = parameters;

return
%%

gc2 = blockproc(data.projections.gc, [2 2], @(block_struct) sum(block_struct.data,[1 2]));
gd2 = blockproc(data.projections.gd, [2 2], @(block_struct) sum(block_struct.data,[1 2]));
rc2 = blockproc(data.projections.rc, [2 2], @(block_struct) sum(block_struct.data,[1 2]));
rd2 = blockproc(data.projections.rd, [2 2], @(block_struct) sum(block_struct.data,[1 2]));

mask = imopen(imclose(data.final_template.F_2D > median(data.final_template.F_2D,'all'),ones(5)),ones(5));
mask = mask(1:2:end,1:2:end);


f = {'fwdaligned','bckaligned'};
c = {'green', 'red','rat'};
tf = {'tfwd', 'tbck'};
counts = {gc2,rc2};
dwell = {gd2,rd2};
us = 4;
for i = 1:2
    
    figData.(f{i}).xaxis = data.final_template.u(1:2:end);
    figData.(f{i}).yaxis = data.final_template.v(1:2:end);
    figData.(f{i}).mask = mask;
    figData.(f{i}).phase_axis = deg2rad(7.5:15:360);    
    for j = 1:2
        [~,ca] =  data.tsr.alignDataToMotionPhase(data.projections.tx, counts{j},figData.(f{i}).phase_axis, parameters.(tf{i}));
        [~,da] =  data.tsr.alignDataToMotionPhase(data.projections.tx, dwell{j},figData.(f{i}).phase_axis, parameters.(tf{i}));
        ca(da <= 0) = 0;
        figData.(f{i}).(c{j}) = squeeze(sum(ca,3,'omitnan'))./squeeze(sum(da,3, 'omitnan'));
        figData.(f{i}).(c{j})(~isfinite(figData.(f{i}).(c{j}))) = 0;
        figData.(f{i}).([c{j} '_stride_with_nan']) = ca./da;
        figData.(f{i}).([c{j} '_stride']) = nan_infill_nd(ca./da);
    end
    if ~any(figData.(f{i}).green_stride_with_nan > 0,'all')
        continue; %in case there is no stride data in one direction
    end
    gmin = percentile(figData.(f{i}).green_stride_with_nan(mask & (figData.(f{i}).green_stride_with_nan > 0)), .01); 
    gmax = percentile(figData.(f{i}).green_stride_with_nan(mask & (figData.(f{i}).green_stride_with_nan > 0)), .999); 
    rmin = percentile(figData.(f{i}).red_stride_with_nan(mask & (figData.(f{i}).red_stride_with_nan > 0)), .01); 
    figData.(f{i}).rat_stride = max(min(figData.(f{i}).green_stride, gmax),gmin)./max(mean(figData.(f{i}).red,3,'omitnan'), rmin);
    figData.(f{i}).rat = max(min(figData.(f{i}).green, gmax),gmin)./max(mean(figData.(f{i}).red,3,'omitnan'), rmin);


    %commented out because slow and being replaced by fourier approach
%     for j = 1:3
%         s = max(figData.(f{i}).(c{j}),[],3);
%         lim =  [0 percentile(s(mask & isfinite(s)), .97)];
%         [figData.(f{i}).([c{j} '_signal']), figData.(f{i}).([c{j} '_phi']), figData.(f{i}).([c{j} '_amplitude']),~,~, figData.(f{i}).([c{j} '_proj_signal'])] = align_activity_nnmf(mask.*figData.(f{i}).(c{j}), us, lim);
%     end    
end

for i = 1:2
    for j = 1:2
        ft = fftn(figData.(f{i}).(c{j})); %fourier transform of mean rate
        ft_tr = fftn(figData.(f{i}).(c{j})(:,:,end:-1:1)); %fourier transrom of mean time - reverse rate
        ac = ifftn(conj(ft).*ft); ac = ac/max(ac,[],'all');
        ac_tr = ifftn(conj(ft_tr).*ft_tr); ac_tr = ac_tr/max(ac_tr,[],'all');
        xc = ifftn(conj(fftn(ac - ac_tr)).*ft);
        xc_norm = ifftn(conj(fftn(abs(ac - ac_tr))).*ft);

        figData.(f{i}).([c{j} '_mot_corr']) = sum(figData.(f{i}).(c{j}).*xc,3)./(sum(xc_norm,3,'omitnan').*mean(figData.(f{i}).(c{j}),'all'));
        %figData.(f{i}).([c{j} '_mot_corr_alt_norm']) = sum(figData.(f{i}).(c{j}).*xc,3)./(sum(figData.(f{i}).(c{j}).*xc_norm,3,'omitnan'));
        
    end
end

mean_no_nan = @(x) mean(x,'all','omitnan');

if (~isfield(parameters,'pca_xrange') || isempty(parameters.pca_xrange))
    parameters.pca_xrange = [-Inf Inf];
end
if (~isfield(parameters,'pca_yrange') || isempty(parameters.pca_yrange))
    parameters.pca_yrange = [-Inf Inf];
end


for i = 1:2
    xi = figData.(f{i}).xaxis >= min(parameters.pca_xrange) & figData.(f{i}).xaxis <= max(parameters.pca_xrange);
    yi = figData.(f{i}).yaxis >= min(parameters.pca_yrange) & figData.(f{i}).yaxis <= max(parameters.pca_yrange);
    figData.(f{i}).pca_xinds = xi;
    figData.(f{i}).pca_yinds = yi;
    figData.(f{i}).pca_xaxis = figData.(f{i}).xaxis(xi);
    figData.(f{i}).pca_yaxis = figData.(f{i}).yaxis(yi);
    
    for j = 1:2
        
        y = figData.(f{i}).(c{j})(xi,yi,:);
        my = mean(y,3);
        y = y - my;
        [coeff,score] = pca(reshape(y,[],size(y,3)), 'NumComponents',2,'Economy',true,'Centered',false);
        xc = xcorr(coeff(:,1), coeff(:,2));
        if (diff(xc(size(y,3) + [0 1])) > 0)
            coeff = coeff(:,[2 1]);
            score = score(:,[2 1]);
        end
        figData.(f{i}).([c{j} '_pca_coeff']) = coeff;
        figData.(f{i}).([c{j} '_pca_z']) = reshape(score(:,1) + 1i*score(:,2), size(y,[1 2]));
        figData.(f{i}).([c{j} '_pca_reconstr']) = reshape(score(:,1), size(y,[1 2])).*reshape(coeff(:,1),1,1,[]) +  reshape(score(:,2), size(y,[1 2])).*reshape(coeff(:,2),1,1,[]) + my;
        [k,~,zfit] = line_fit_phase_z(figData.(f{i}).yaxis(yi), figData.(f{i}).xaxis(xi), figData.(f{i}).([c{j} '_pca_z']));
        figData.(f{i}).([c{j} '_pca_kvector']) = k([2 1]); %note reversal to make kx,ky instead of y,x 
        figData.(f{i}).([c{j} '_pca_wave_z']) = zfit;
        figData.(f{i}).([c{j} '_pca_wave_reconstr']) = real(zfit).*reshape(coeff(:,1),1,1,[]) + imag(zfit).*reshape(coeff(:,2),1,1,[]) + my;
        figData.(f{i}).([c{j} '_var_expl_pca']) = 1 - var(figData.(f{i}).(c{j})(xi,yi,:)-figData.(f{i}).([c{j} '_pca_reconstr']),1,'all')./var(y,1,'all');
        figData.(f{i}).([c{j} '_var_expl_pca_wave']) = 1 - var(figData.(f{i}).(c{j})(xi,yi,:)-figData.(f{i}).([c{j} '_pca_wave_reconstr']),1,'all')./var(y,1,'all');     
        [temp,figData.(f{i}).proj_axis,~,~,~,dproj] = fun_in_region_along_vector(figData.(f{i}).xaxis(xi), figData.(f{i}).yaxis(yi), figData.(f{i}).([c{j} '_pca_z']), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, mean_no_nan);
        figData.(f{i}).([c{j} '_phase_vs_proj']) = angle(mean(temp,2,'omitnan'));
    end
end


return
c = {'green', 'red','rat'};
mean_no_nan = @(x) mean(x,'all','omitnan');
deg = 1;
for i = 1:2
     if ~any(figData.(f{i}).green_stride_with_nan > 0,'all')
        continue; %in case there is no stride data in one direction
    end
    for j = 1:3
        [temp,figData.(f{i}).proj_axis,~,~,~,dproj] = fun_in_region_along_vector(figData.(f{i}).xaxis, figData.(f{i}).yaxis, mask.*exp(1i*figData.(f{i}).([c{j} '_phi'])).*figData.(f{i}).([c{j} '_amplitude']), parameters.proj_pt1, parameters.proj_pt2, parameters.proj_dist, mean_no_nan);
        figData.(f{i}).([c{j} '_phase_vs_proj']) = angle(mean(temp,2,'omitnan'));
        if (j == 3)
            m = double(mask).*mean(figData.(f{i}).red_proj_signal, 3, 'omitnan');
        else
            m = double(mask);
        end
        [figData.(f{i}).([c{j} '_polyfit']),figData.(f{i}).([c{j} '_polyfit_projection'])] = fit_signal_polynomial(deg, figData.(f{i}).phase_axis, upsample_circular(figData.(f{i}).([c{j} '_signal']), 1/us), dproj, figData.(f{i}).([c{j} '_amplitude']), figData.(f{i}).(c{j})-min(figData.(f{i}).(c{j}),[],3), m);
        figData.(f{i}).([c{j} '_polyfit_projection']) = figData.(f{i}).([c{j} '_polyfit_projection']) + min(figData.(f{i}).(c{j}),[],3);
        %        [figData.(f{i}).([c{j} '_polyfit']), figData.(f{i}).([c{j} '_polyfit_hessian'])] =  fitLineAngle(dproj, figData.(f{i}).([c{j} '_phi']),  mask.*figData.(f{i}).([c{j} '_amplitude']));
        temp = unwrap([polyval(figData.(f{i}).([c{j} '_polyfit']), figData.(f{i}).proj_axis)' figData.(f{i}).([c{j} '_phase_vs_proj'])],[],2);
        figData.(f{i}).([c{j} '_phase_vs_proj']) = temp(:,2);
        figData.(f{i}).([c{j} '_var_explained']) = 1-sum(mask.*(figData.(f{i}).(c{j})-figData.(f{i}).([c{j} '_polyfit_projection'])).^2,'all','omitnan')/ sum(mask.*(figData.(f{i}).(c{j})-mean(figData.(f{i}).(c{j}),3,'omitnan')).^2,'all','omitnan');
        dphase = polyval(figData.(f{i}).([c{j} '_polyfit']),0);
        figData.(f{i}).([c{j} '_phase_vs_proj']) = figData.(f{i}).([c{j} '_phase_vs_proj']) - dphase;
        figData.(f{i}).([c{j} '_polyfit'])(end) = figData.(f{i}).([c{j} '_polyfit'])(end) - dphase;
        figData.(f{i}).([c{j} '_phi']) = figData.(f{i}).([c{j} '_phi']) - dphase;
    end
end

return;


