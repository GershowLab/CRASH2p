function pipeline_rigid_registration(srcdir,regset,istest)
% rigid registration

setupDirectories;

existsAndDefault('regset',[]);
existsAndDefault('istest',false);

ws = warning; % save warning settings to restore later
warning('off','verbose');
warning('off','backtrace');

% set up default registration settings: entire tracking duration as one track,
% auto-detect template range
defaultregset.rigid_trange = [];
defaultregset.template_trange = [];
[~,timestamp,~] = fileparts(srcdir);
fstub = 'registered'; % name of default or user-supplied dstdir

% detect multiple tracks
if ~isempty(regset)
    disp('using user-supplied registration settings');
    fn = fieldnames(defaultregset);
    for i=1:length(fn)
        if ~isfield(regset,fn{i}) || (isfield(regset,fn{i}) && isempty(regset.(fn{i})))
            regset.(fn{i}) = defaultregset.(fn{i});
        end
    end
    doRigidRegistration(srcdir,fstub,timestamp,regset,istest);
else
    d = dir(fullfile(srcdir,'registered*'));
    if isempty(d)
        disp('using default registration settings');
        regset = defaultregset;
        doRigidRegistration(srcdir,fstub,timestamp,regset,istest);
    else
        disp('using pre-saved registration settings');
        ntracks = numel(d);
        for i=1:ntracks
            fprintf('registering track %d/%d (%s/)\n',i,ntracks,d(i).name);
            load(fullfile(d(i).folder,d(i).name,[timestamp '_registration_settings.mat']),'regset');
            fstub = d(i).name;
            doRigidRegistration(srcdir,fstub,timestamp,regset,istest);
        end
    end
end

warning(ws); % restore previous warning settings

end
    
function doRigidRegistration(srcdir,fstub,timestamp,regset,istest)
    
    %-------------
    % parse input
    %-------------

    if istest
        regset.rigid_trange = regset.template_trange(1)+[-2.5 2.5];
    end
    
    dstdir = fullfile(srcdir,fstub);
    if isunix
        dstdir = strrep(dstdir,'data_in','data_out');
    end
    if ~isfolder(dstdir)
        mkdir(dstdir);
    end
    
    if ispc
        fpath = fullfile(dstdir,[timestamp '_rigid_registration_log.txt']);
        if isfile(fpath)
            delete(fpath);
        end
        diary(fpath);
    end
        
    ts = tic;
    fprintf('%s beginning rigid registration on data set %s...\n',datestr(now),timestamp);
    
    %-----------
    % prep work
    %-----------
    
    % set up vraligner params
    fprintf([datestr(now) ' setting up parameters for rigid registration...\n']);
    vraopts = VRAligner.processOptions;
    if ~isempty(regset.rigid_trange)
        vraopts.trange_alignment = regset.rigid_trange;
    end
    if ~isempty(regset.template_trange)
        vraopts.template_timerange = regset.template_trange;
    end
    vraopts.dstdir = dstdir;
    vraopts.save_frames = false;
    
    %---------------------
    % rigid registration
    %---------------------
    
    % do rigid registration, save vra.mat and other relavent data
    fprintf([datestr(now) ' beginning rigid registration and saving results...\n']);
    VRAligner.processDirectory(srcdir,vraopts);
    
    %-----------
    % post work
    %-----------
    
    % load registration results
    load(fullfile(srcdir,[timestamp '_vr.mat']),'vr');
    vr = vr.changeFilename(srcdir);
    vr.vrfa.vrfc.unload;
    load(fullfile(dstdir,[timestamp '_vra.mat']),'vra');
    vra.vr = vr;
    load(fullfile(dstdir,[timestamp '_vra_options.mat']),'vraopts');
    
    % update and save registration settings
    if isempty(regset.rigid_trange)
        regset.rigid_trange = vraopts.trange_alignment;
    end
    if isempty(regset.template_trange)
        regset.template_trange = vraopts.template_timerange;
    end
    fpath = fullfile(dstdir,sprintf('%s_registration_settings.mat',timestamp));
    save(fpath,'regset');
    
    fprintf([datestr(now) ' ...all done!\n']);
    toc(ts);
    
    if ispc
        diary off;
    end
        
end