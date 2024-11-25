% CRASH2p image registration and data analysis pipeline
% start: .bin, .txt, and .json files created by the LabVIEW data
%        acquisition software
% end: VOI and PCA analysis figures as shown in the paper

setupDirectories; % add entire code repo to matlab path; only need to do once per MATLAB session

%% 1) parse binary hyperscope data and save as .mat

% if srcdir is not specified in base workspace, look for the dataset folder
% named with the 14-digit TIMESTAMP in pwd/demo_data/
if ~exist('srcdir','var')
    try
        d = dir(fullfile(pwd,'demo_data'));
        d(~[d.isdir]) = []; % filter for folders only
        validinds = cellfun(@(x,y) ~isempty(regexp(x,y,'once')),{d.name},repmat({'\d{14}'},size({d.name})));
        d = d(validinds); % filter for 14-digit folder names
        srcdir = fullfile(d(1).folder,d(1).name); % use first dataset folder if multiple
        fprintf('found demo dataset %s\n',srcdir);
    catch
        warning('ATTENTION: no demo dataset found in current working folder; please specify full dataset folder path as srcdir in base workspace');
        return;
    end
else
    fprintf('using user-specified demo dataset %s\n',srcdir);
end

% load and parse vr (photon arrival and tracker data)
% try loading vr from .mat file
[~,timestamp,~] = fileparts(srcdir);
fpath = fullfile(srcdir,sprintf('%s_vr.mat',timestamp));
if isfile(fpath)
    load(fpath);
    vr = vr.changeFilename(srcdir);
    vr.vrfa.vrfc.unload;
else
    % parse vr from .bin file and make necessary calculations
    vr = parseHyperscopeDataFromBinary(srcdir);
end

% will (might) add the following new files to srcdir/ if running for the first time:
%   TIMESTAMP_vr.mat
%   (TIMESTAMP_parse_error_message.txt)

%% 1a) optional preview movie

saveMovie(srcdir,[],'preview');

% will add the following new files to srcdir/:
%   TIMESTAMP_preview_movie.mp4
%   TIMESTAMP_preview_movie_options.mat

%% 2) rigid registration

% determine time ranges manually
regset.rigid_trange = [3999 4255]; % determine from preview movie, exclude when larva is compressed
regset.template_trange = [4103 4104]; % determine from preview movie, choose 1-s interval
[~,timestamp,~] = fileparts(srcdir);
if ~exist(fullfile(srcdir,'registered'),'dir' )
    mkdir(fullfile(srcdir,'registered'));
end
fpath = fullfile(srcdir,'registered',sprintf('%s_registration_settings.mat',timestamp));
save(fpath,'regset');

pipeline_rigid_registration(srcdir);

% will add the following new files to srcdir/registered/:
%   TIMESTAMP_initial2Dalignment.mat
%   TIMESTAMP_registration_settings.mat
%   TIMESTAMP_rigid_registration_log.txt
%   TIMESTAMP_template.mat
%   TIMESTAMP_vra.mat
%   TIMESTAMP_vra_options.mat

%% 2a) optional rigid-registered movie

saveMovie(srcdir,'registered','rigid');

% will add the following new files to srcdir/registered/:
%   TIMESTAMP_rigid_movie.mp4
%   TIMESTAMP_rigid_movie_options.mat

%% 3) MANUAL STEP: finalize template for non-rigid registration and intensity correction

% follow script below and run each section of code separately
open pipeline_after_rigid.m

% will add the following new files to srcdir/registered/:
%   TIMESTAMP_final_template_options.mat
% will (might) update the following files in srcdir/registered/:
%   TIMESTAMP_registration_settings.mat
%   TIMESTAMP_template.mat
%   (TIMESTAMP_vra.mat)

%% 4) non-rigid registration and intensity correction

% default: use same rigid time ranges for non-rigid
[~,timestamp,~] = fileparts(srcdir);
fpath = fullfile(srcdir,'registered',sprintf('%s_registration_settings.mat',timestamp));
load(fpath,'regset');
regset.non_rigid_trange = regset.rigid_trange;
save(fpath,'regset');

pipeline_non_rigid_registration(srcdir);

% will add the following new files to srcdir/registered/:
%   TIMESTAMP_cdmffd_options.mat
%   TIMESTAMP_metric.mat
%   TIMESTAMP_movie_projections_full.mat
%   TIMESTAMP_non_rigid_registration_log.txt
% will add the follwing new files to srcdir/registered/registered_frames/:
%   TIMESTAMP_frame*.mat
% will modify the follwing files in srcdir/registered/:
%   TIMESTAMP_registration_settings.mat
%   TIMESTAMP_vra.mat

%% 4a) optional final movie

saveMovie(srcdir,'registered','final');

% will add the following new files to srcdir/registered/:
%   TIMESTAMP_movie_projections_full_movie.mp4
%   TIMESTAMP_movie_projections_full_movie_options.mp4

%% 5.1) MANUAL STEP: draw VOIs using VOI selector GUI

% use this GUI to select spherical or cylindrical VOIs for trace extraction
% after non-rigid
% - spherical VOIs: move crosshair to center of sphere then click "select"
% - cylindrical VOIs: move crosshair to center of top surface of cylinder
%                     then click "select"; move crosshair to center of
%                     bottom surface of cylinder then click "select" again

disp('VOISelector app is running...');

[~,timestamp,~] = fileparts(srcdir);
templatefpath = fullfile(srcdir,'registered',sprintf('%s_template.mat',timestamp));
vs = VOISelectorObject(templatefpath);
VOISelector(vs);

disp('...done!');

% will add the following new file to srcdir/registered/ (file name
% modifiable in GUI):
%   voi_selections.mat

%% 5.2) extract VOI activity

% extract activity from saved voisets and frames, then save activity to
% disk
vsfpath = fullfile(srcdir,'registered','voi_selections.mat');
activity = addActivityInVOIsToDirectory(fullfile(srcdir,'registered'),vsfpath);

% will add the following new file to srcdir/registered/:
%   voi_activity.mat

%% 6.1) prep for figures 

% ATTENTION: the following 3 sections reproduce Figures 2, S1, 3, and S2
% from the CRASH2p paper using the same A27h Gcamp6f dataset; in addition
% to files generated by the above steps, SLEAP labels are required with the
% following file path pattern for auto-detection:
%   srcdir/behavior_video/*.analysis.h5
% or specify expt.sleap_path below to supply file location

% prep manual annotations for sample dataset
expt.srcdir = srcdir;
expt.register = 'registered';
expt.sleap_path = []; % auto-detect
expt.tfwd = [3995 4013;
             4022 4027;
             4028 4068;
             4069 4110;
             4153 4162;
             4166 4183;
             4186 4188;
             4240 4244;
             4248 4254];
expt.ttrn = [4013 4022;
             4027 4028;
             4068 4069;
             4110 4115;
             4141 4153;
             4162 4166;
             4183 4186;
             4232 4240;
             4244 4248];
expt.tbck = [4115 4141;
             4188 4232];
expt.bckTrange = [4190 4230];
expt.fwdTrange = [4028 4068];
expt.vs = 6;
expt.traceinds = [];

% save xyt projections for analysis
add_xytprojections_to_dir(expt.srcdir,expt.register);

% load sample dataset
data = loadFigureData_activity_summary(expt.srcdir,expt.register,expt.sleap_path);
figData = gatherFigureData_activity_summary(data,expt);
figData.templatez = data.final_template.w;
figData.xtproj.zrange = data.projections.zrange;

% will add the following new file to srcdir/registered/:
%   TIMESTAMP_xytprojections.mat

%% 6.2) VOI figures

% prep common figure options
figOp_voi = activity_summary_figure_vnc;
figOp_voi.activity_range = [0.4 5];
figOp_voi.norm_rat_range = [.5 3.5];
figOp_voi.activity_linewidth_zoom = 1.5;
figOp_voi.voi_colors = winter(ceil(1.1*numel(figData.vois)))*1;

% paper Figure 2
figOp_voi.fignum = 1;
ah_voi = activity_summary_figure_vnc(figData,figOp_voi);

% paper Supplemental Figure S1
figOp_voi.fignum = 2;
figOp_voi.data_channel = 'red_uc';
ah_voi_red = activity_summary_figure_vnc(figData,figOp_voi);
set([ah_voi_red([end-1, end]).ax], 'YLim', ah_voi(end).ax.YLim);

%% 6.3) PCA figures

% prep common figure options
figOp_pca = activity_pca_figure_sad_vnc();
figOp_pca.glim = [0 1.5e6];
figOp_pca.rlim = [0 5e6];
figOp_pca.dovois = true;
figOp_pca.voi_colors = [240 150 100; 0 255 255]/255;
figOp_pca.activity_range = [0 6];
figOp_pca.activity_overlap = 0.2;
fd = figData;
fd.activityPlot.y = (lowpass1D(data.activity.voisets(5).lpr(10,:)',2))'./data.activity.voisets(5).r0(10);
fd.activityPlot.y(2,:) = (lowpass1D(data.activity.voisets(6).lpr(4,:)',2))'./data.activity.voisets(6).r0(4);

% paper Figure 3
% prep data
gc = (lowpass1D(data.activity.voisets(5).gcfine(10,:)',2))';
gd = (lowpass1D(data.activity.voisets(5).gdfine_ic(10,:)',2))';
fd.activityPlot.green = gc./gd;
gc = (lowpass1D(data.activity.voisets(6).gcfine(4,:)',2))';
gd = (lowpass1D(data.activity.voisets(6).gdfine_ic(4,:)',2))';
fd.activityPlot.green(2,:) = gc./gd;
fd.vois = {data.voisets(5).vois{10}, data.voisets(6).vois{4}};
% plot figure
figOp_pca.color = 1;
figOp_pca.fignum = 3;
figOp_pca.voi_plot_ratio = true;
ah_pca = activity_pca_figure_sad_vnc(fd,figOp_pca);

% paper Supplemental Figure S2
% prep data
rc = (lowpass1D(data.activity.voisets(5).rcfine(10,:)',2))';
rd = (lowpass1D(data.activity.voisets(5).rdfine_ic(10,:)',2))';
fd.activityPlot.red = rc./rd;
rc = (lowpass1D(data.activity.voisets(6).rcfine(4,:)',2))';
rd = (lowpass1D(data.activity.voisets(6).rdfine_ic(4,:)',2))';
fd.activityPlot.red(2,:) = rc./rd;
% plot figure
figOp_pca.color = 2;
figOp_pca.fignum = 4;
figOp_pca.voi_plot_ratio = false;
ah_pca_red = activity_pca_figure_sad_vnc(fd,figOp_pca);