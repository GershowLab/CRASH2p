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
