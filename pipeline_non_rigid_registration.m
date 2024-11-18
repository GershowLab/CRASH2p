function pipeline_non_rigid_registration(srcdir,istest,skipnr)
% non-rigid registration and intensity correction

setupDirectories;

existsAndDefault('istest',false);
existsAndDefault('skipnr',false);

if istest
    disp('this is a test');
end

ws = warning; % save warning settings to restore later
warning('off','verbose');
warning('off','backtrace');

% detect multiple tracks
d = dir(fullfile(srcdir,'registered*'));
if istest
    ntracks = 1;
else
    ntracks = numel(d);
end
for i=1:ntracks
    [~,timestamp,~] = fileparts(srcdir);
    load(fullfile(d(i).folder,d(i).name,[timestamp '_registration_settings.mat']),'regset');
    if ~isfield(regset,'non_rigid_trange')
        fprintf('track %d/%d (%s/) not properly prepped for non-rigid, skipped\n',i,ntracks,d(i).name);
        continue;
    end
    if any(isnan(regset.non_rigid_trange))
        fprintf('track %d/%d (%s/) is marked as invalid, skipped\n',i,ntracks,d(i).name);
        continue;
    else
        fprintf('registering track %d/%d (%s/)\n',i,ntracks,d(i).name);
        subdir = d(i).name;
        doNonRigidRegistration(srcdir,subdir,timestamp,regset,istest,skipnr);
    end
end

warning(ws); % restore previous warning settings

end

function doNonRigidRegistration(srcdir,subdir,timestamp,regset,istest,skipnr)

%-------------
% parse input
%-------------

if istest
    non_rigid_trange = regset.non_rigid_trange(1)+[0 5];
else
    non_rigid_trange = regset.non_rigid_trange;
end

if isunix
    dstdir = strrep(srcdir,'data_in','data_out');
    if ~isfolder(dstdir)
        mkdir(dstdir);
    end
elseif ispc
    dstdir = srcdir;
end
if ~isfolder(fullfile(dstdir,subdir))
    mkdir(dstdir,subdir);
end

if ispc
    fpath = fullfile(dstdir,subdir,[timestamp '_non_rigid_registration_log.txt']);
    if isfile(fpath)
        delete(fpath);
    end
    diary(fpath);
end

tic;
if skipnr
    fprintf('%s beginning intensity correction on data set %s...\n',datestr(now),timestamp);    
else
    fprintf('%s beginning non-rigid registration and intensity correction on data set %s...\n',datestr(now),timestamp);
end

%-----------
% prep work
%-----------

% load vra
load(fullfile(srcdir,[timestamp '_vr.mat']),'vr');
vr = vr.changeFilename(srcdir);
vr.vrfa.vrfc.unload;
% for backward compatibility
vr = vr.applyScanOffset;
load(fullfile(srcdir,subdir,[timestamp '_vra.mat']),'vra');
vra.vr = vr;
load(fullfile(srcdir,subdir,[timestamp '_vra_options.mat']),'vraopts');
% for backward compatibility
if isempty(vra.valid)
    vra.valid = false(1,length(vra.vr.frame.edges)-1);
end
if all(~vra.valid)
    vra.valid(vraopts.alignedInds) = true;
end

% load final template and binning options
load(fullfile(srcdir,subdir,[timestamp '_final_template_options.mat']),'templateopts');
template = templateopts.template.F;
smooth_s = templateopts.smooth_s;
rmin = templateopts.analysis_mask.threshold;
bw = templateopts.analysis_mask.mask;
fa = templateopts.fa;
[xc,yc,zc] = fa.getBinCenters;

disp('rigid registration results loaded');

% trim frame range according to non_rigid_trange
validinds = [];
for i=1:size(non_rigid_trange,1)
    validinds = [validinds find(vr.frame.edges>=non_rigid_trange(i,1) & vr.frame.edges<=non_rigid_trange(i,2))]; %#ok<*AGROW>
end
frameinds = intersect(find(vra.valid),validinds);
framestart = min(frameinds);
framestop = max(frameinds);
% (this should also trim non_rigid_trange to be within rigid_trange)

% for backward compatibility: align subframe if needed
if any(~vra.hassubframe(frameinds))
    fprintf([datestr(now) ' subframes rigid registration not done, doing now...\n']);
    for i=frameinds
        if ~vra.hassubframe(i)
            vra = vra.fineAlignSubFrame(i);
        end
    end
    fprintf([datestr(now) ' ...done!\n']);
end

% set up frame saving output path
framedstdir = fullfile(dstdir,subdir,'registered_frames');
if ~isfolder(framedstdir)
    mkdir(framedstdir);
end
framefstub = fullfile(framedstdir,sprintf('%s_frame',timestamp));

% save rigid frames
fprintf([datestr(now) ' saving rigid registered frames...\n']);
parfor ii=1:length(frameinds)
    i = frameinds(ii);
    [rc,rd] = vra.binSubFrame(i,fa,'r'); %#ok<PFBNS>
    [gc,gd] = vra.binSubFrame(i,fa,'g');
    te = vra.subframe_time_edges(i,:);
    if vra.valid(i)
        fpath = sprintf('%s%d.mat',framefstub,i);
        a = matfile(fpath,'Writable',true);
        a.rigid_red_counts = rc;
        a.rigid_red_dwell = rd;
        a.rigid_green_counts = gc;
        a.rigid_green_dwell = gd;
        a.time_edges = te;
        a.xaxis = xc;
        a.yaxis = yc;
        a.zaxis = zc;
    end
    rc = []; rd = []; gc = []; gd = []; %#ok<NASGU> % release memory; see Rui's Evernote journal 20230405
end
fprintf([datestr(now) ' ...done!\n']);

% save updated vra and options
vra.accumulated_images = [];
vra.vr = [];
save(fullfile(dstdir,subdir,[timestamp '_vra.mat']),'vra','-v7.3');
disp('vra updated and saved');
vra.vr = vr;

%------------------------
% non-rigid registration
%------------------------

% set up cdmffd params
cdmffdopts.nLevel = 3;
cdmffdopts.maxIter = 100;
cdmffdopts.TOL = 0.05;
cdmffdopts.useGaussian = false;
if ~isfield(regset,'reg_gamma') || (isfield(regset,'reg_gamma') && isempty(regset.reg_gamma))
    cdmffdopts.regGamma = 1;
else
    cdmffdopts.regGamma = regset.reg_gamma;
end

% prep for saving image metrics
nframe = length(vr.frame.edges)-1;
m0 = NaN(1,nframe);
m1 = NaN(1,nframe);

fprintf([datestr(now) ' beginning non-rigid registration...\n']);
parfor i=framestart:framestop
    % load vra frame
    fpath = sprintf('%s%d.mat',framefstub,i);
    try
        if ~isfile(fpath)
            fprintf('skipped invalid frame %d',i);
            continue;
        end
        a = matfile(fpath,'Writable',true);
        if skipnr
            a.non_rigid_Tx = zeros(size(template));
            a.non_rigid_Ty = zeros(size(template));
            a.non_rigid_Tz = zeros(size(template));
            a.non_rigid_red_counts = a.rigid_red_counts;
            a.non_rigid_red_dwell = a.rigid_red_dwell;
            a.non_rigid_green_counts = a.rigid_green_counts;
            a.non_rigid_green_dwell = a.rigid_green_dwell;
        else
            ct0_r = squeeze(sum(a.rigid_red_counts,4));
            dw0_r = squeeze(sum(a.rigid_red_dwell,4));
            mmat = imclose(dw0_r>0,ones([9 9 9]));
            rv0 = l1spline(ct0_r./dw0_r, smooth_s, 1, 100, 1, 1e-5);
            rv0 = mmat.*rv0;
            % get vra image metric
            m0(i) = multissim3(rv0,template);
            % do cdmffd
            rv0_rmin = rv0;
            rv0_rmin(~bw) = rmin;
            rv0_rmin(rv0_rmin<rmin) = rmin;
            mmat = imclose(dw0_r>0,ones([9 9 9]));
            [~,Vx,Vy,Vz] = MultiresolutionRegistration3D_silent(rv0_rmin,mmat.*template,...
                cdmffdopts.nLevel,cdmffdopts.maxIter,false,false,false,false,...
                cdmffdopts.TOL,cdmffdopts.useGaussian,cdmffdopts.regGamma); %#ok<PFBNS>
            % save cdmffd results
            a.non_rigid_Tx = Vx;
            a.non_rigid_Ty = Vy;
            a.non_rigid_Tz = Vz;
            % save count/dwell with subframes
            [ct1_r,dw1_r] = getCDMFFDSubFrame(vra,fa,i,struct('Tx',Vx,'Ty',Vy,'Tz',Vz),'r');
            [ct1_g,dw1_g] = getCDMFFDSubFrame(vra,fa,i,struct('Tx',Vx,'Ty',Vy,'Tz',Vz),'g');
            a.non_rigid_red_counts = ct1_r;
            a.non_rigid_red_dwell = dw1_r;
            a.non_rigid_green_counts = ct1_g;
            a.non_rigid_green_dwell = dw1_g;
            mmat = imclose(sum(dw1_r,4)>0,ones([9 9 9]));
            rv1 = l1spline(sum(ct1_r,4)./sum(dw1_r,4), smooth_s, 1, 100, 1, 1e-5);
            rv1 = mmat.*rv1;
            % get cdmffd metric
            m1(i) = multissim3(rv1,template);
            fclose all;
        end
    catch me
        fprintf('ERROR at %s',fpath);
        fclose all;
        rethrow(me);
    end
end
fprintf([datestr(now) ' ...done!\n']);

if ~skipnr
    metric.rigid = m0;
    metric.non_rigid = m1;
    save(fullfile(dstdir,subdir,[timestamp '_metric.mat']),'metric');
    disp('image metrics saved');
    save(fullfile(dstdir,subdir,[timestamp '_cdmffd_options.mat']),'cdmffdopts');
    disp('non-rigid registration options saved');
end

regset.saved_frame_inds = frameinds;
regset.voxel_size = [templateopts.xybin templateopts.xybin templateopts.zbin];
regset.smooth_s = smooth_s;
regset.r_min = rmin;
if ~skipnr
    regset.reg_gamma = cdmffdopts.regGamma;
end
save(fullfile(dstdir,subdir,sprintf('%s_registration_settings.mat',timestamp)),'regset');
disp('registration settings updated and saved');

%----------------------
% intensity correction
%----------------------

fprintf([datestr(now) ' beginning intensity correction...\n']);
if isunix
    d = dir(fullfile(srcdir,subdir,'*_template.mat'));
    fpath1 = fullfile(srcdir,subdir,d(1).name);
    fpath2 = fullfile(dstdir,subdir,d(1).name);
    copyfile(fpath1,fpath2);
end
IntensityCorrectorBSpline.correct_directory(dstdir,subdir,true);
fprintf([datestr(now) ' ...done!\n']);

%---------------------
% save 2D projections
%---------------------

fprintf([datestr(now) ' saving 2D projections for movie...\n']);
projopts.sigma = [1.3 1.3 1.3 1.3];
projopts.savefile = true;
projopts.calc_mip = false;
projopts.savename = sprintf('%s_movie_projections_full.mat',timestamp);
get_projections_for_movie(dstdir,subdir,[],[],[],projopts);
fprintf([datestr(now) ' ...done!\n']);

fprintf([datestr(now) ' ...all done!\n']);
toc;

if ispc
    diary off;
end

end