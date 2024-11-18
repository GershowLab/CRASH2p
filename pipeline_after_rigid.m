% image registration pipeline: evaluate rigid registration results and make
% template for non-rigid registration

% ATTENTION: run each section of code separately

%% ==========BELOW: for each data set (or each track within a data set)==========
%  set srcdir:
%  e.g. srcdir = 'R:\hyperscope registered data\A27h_gCamp7f_mCherry\tracking\20220204113228';
%  set tracknum accordingly if got subfolder registered_1/2/... in srcdir
%  e.g. tracknum = 1;

%% load rigid registration result

if ~exist('tracknum','var')
    tracknum = [];
end

% get some additional variables make finding files easier
[~,timestamp,~] = fileparts(srcdir);
if isempty(tracknum)
    subdir = 'registered';
else
    subdir = sprintf('registered_%d',tracknum);
end

% load vr
fpath = fullfile(srcdir,sprintf('%s_vr.mat',timestamp));
load(fpath,'vr');
vr = vr.changeFilename(srcdir);
vr.vrfa.vrfc.unload;
disp('vr loaded from mat file');

% load vra
load(fullfile(srcdir,subdir,[timestamp '_vra.mat']),'vra');
vra.vr = vr;
disp('vra loaded from mat file');

disp('...done!');

%% run volume rotation GUI if necessary

disp('VolumeRotator app is running...');

vrot = VolumeRotator(vra);
uiwait(vrot.RotateVolumeUIFigure);

% apply to vra
vra.additional_rotation = rotation_angle_deg;

fprintf('...done!\nAdditional rotations:\n\tX: %.1f degrees\n\tY: %.1f degrees\n\tZ: %.1f degrees\n',rotation_angle_deg);

%% run template setup GUI

% use this GUI to select binning parameters for template and individual
% frames to be used in hpc pipeline

disp('TemplateCreator app is running...');

[templateopts,tm] = setupTemplate(vra);

% click on "Done" to quit the GUI and output params to workspace

disp('...done!');

%% update saved files and get ready for data transfer to hpc

% save template settings and TemplateMaker object to disk
fpath = fullfile(srcdir,subdir,[timestamp '_final_template_options.mat']);
save(fpath,'templateopts');
disp('template and binning options saved');
% fpath = fullfile(srcdir,fstub,[timestamp '_tm.mat']);
% save(fpath,'tm');
% disp('tm saved');

% update and save vra.mat
vra.accumulated_images = [];
vra.vr = [];
fpath = fullfile(srcdir,subdir,[timestamp '_vra.mat']);
save(fpath,'vra','-v7.3');
disp('vra updated and saved');

% update and save template.mat
final_template = templateopts.template;
fpath = fullfile(srcdir,subdir,[timestamp '_template.mat']);
save(fpath,'final_template','-append');
disp('final template saved');

% update and save registration_settings.mat
fpath = fullfile(srcdir,subdir,[timestamp '_registration_settings.mat']);
load(fpath,'regset');
regset.voxelsize = [templateopts.xybin templateopts.xybin templateopts.zbin];
regset.smooth_s = templateopts.smooth_s;
save(fpath,'regset');
disp('registration settings updated and saved');

disp('...all done! ready for non-rigid registration');
