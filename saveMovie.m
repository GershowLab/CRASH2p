function saveMovie(srcdir,subdir,movietype,istest)
%saveMovie saves 4-panel registered movie (behavior image with tracker
%trajectory overlay, green fluorescence, red fluorescence, and green/red
%ratiometrics) from a saved projections .mat file. Use HyperMovieMaker
%directly for more customizable options.
%
%   Must-have inputs:
%       srcdir:     absolute path ending in 14-digit numerical TIMESTAMP)
%       subdir:     'registered' (or 'registered_*' if srcdir contains multiple tracks)
%
%   Optional inputs:
%       projname:	default = 'TIMESTAMP_movie_projections_full.mat'
%       istest:     default = false, will only save a short 100-frame movie if true

existsAndDefault('movietype','final');
existsAndDefault('istest',false);

tic;

% parse inputs
[~,timestamp,~] = fileparts(srcdir);
if strcmp(movietype(end-3:end),'.mat')
    projfname = movietype;
    movietype = 'final';
else
    projfname = [];
end

% load tsr
fpath = fullfile(srcdir,sprintf('%s_vr.mat',timestamp));
load(fpath,'vr');
vr = vr.changeFilename(srcdir);
vr.vrfa.vrfc.unload;
tsr = vr.tsr;
if isempty(tsr.behavior)
    tsr = tsr.addBehaviorVideo;
end
if isempty(tsr.com)
    tsr = tsr.conditionData;
end
try
    if ~isfield(tsr.behavior,'sleap')
        tsr = tsr.addSleapResult;
    end
    if ~isfield(tsr.com,'periphase')
        tsr = tsr.findCounterMovements;
    end
catch
    disp('failed to load SLEAP behavior labels, traker overlay in behavior panel will have a small offset');
end
vr.tsr = tsr;
% Q: should I save this updated vr to disk?

% load vra/projections depending on movietype
switch movietype
    case 'rigid'
        fpath = fullfile(srcdir,subdir,sprintf('%s_vra.mat',timestamp));
        load(fpath,'vra');
        vra.vr = vr;
        if isempty(vra.valid)
            fpath = fullfile(srcdir,subdir,sprintf('%s_vra_options.mat',timestamp));
            load(fpath,'vraopts');
            % for backward compatibility
            if isempty(vra.valid)
                vra.valid = false(1,length(vra.vr.frame.edges)-1);
            end
            if all(~vra.valid)
                vra.valid(vraopts.alignedInds) = true;
            end
        end
        disp('rigid results loaded');
    case 'final'
        if isempty(projfname)
            projfname = sprintf('%s_movie_projections_full.mat',timestamp);
        end
        fpath = fullfile(srcdir,subdir,projfname);
        try
            proj = load(fpath);
        catch
            warning('Error loading projections file %s, please check file path!',projfname);
            return;
        end
        proj.fname = fpath;
        disp('final projections loaded');
end

% save movie
switch movietype
    case 'preview'
        hmm = HyperMovieMaker_VR(vr,tsr);
        hmm = hmm.set_output_name(fullfile(srcdir,sprintf('%s_preview_movie.mp4',timestamp)));
    case 'rigid'
        hmm = HyperMovieMaker_VRA(vra,tsr);
        hmm = hmm.set_output_name(fullfile(srcdir,subdir,sprintf('%s_rigid_movie.mp4',timestamp)));
    case 'final'
        hmm = HyperMovieMaker_Projections(proj,tsr);
        hmm = hmm.set_output_name(fullfile(srcdir,subdir,strrep(projfname,'.mat','_movie.mp4')));
end
tx = hmm.get_default_time_axis;
if istest
    tx = tx(1:100);
end
hmm = hmm.setup_figure;
if length(tx)<=100
    hmm = hmm.set_limits(tx);
else
    hmm = hmm.set_limits(tx(unique(round(linspace(1,length(tx),100)))));
end
fprintf('saving movie to %s...\n',hmm.outputname);
hmm.play_movie(tx);
fprintf('...done! Took %.1fs\n',toc);

end