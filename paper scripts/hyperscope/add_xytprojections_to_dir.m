function projections = add_xytprojections_to_dir (srcdir, register, redo, zi, appendz)
%srcdir = 'R:\hyperscope registered data\A27h_gCamp7f_mCherry\tracking\20201204151516';
%register = 'registered';
[~,timestamp,~] = fileparts(srcdir);

existsAndDefault('redo',false);
existsAndDefault('zi',[]);
existsAndDefault('appendz', ~isempty(zi));

%don't redo

 if (appendz && ~isempty(zi))
    fname = fullfile(srcdir, register, [timestamp sprintf('_z%dto%d_xytprojections.mat', find(zi, 1 ),find(zi, 1, 'last' ))]);
 else
    fname = fullfile(srcdir, register, [timestamp '_xytprojections.mat']);
 end

if (~redo && exist(fname, 'file'))
    if (nargout > 0)
        load(fname, 'projections');
    end
    return;
end

if (exist(fname, 'file'))
    [dd,ff] = fileparts(fname);
    movefile(fname, fullfile(dd, [ff '_backup' num2str(redo) '.mat']));
end

try
    if (isempty(zi))
        data = load(fullfile(srcdir,register,'voi_selections.mat'));
        v = [data.voisets.vois];
        bwim = v{1}.makeMask();
        for j = 1:length(v)
            bwim = bwim | v{j}.makeMask();
        end
        zi = squeeze(any(bwim,[1 2]));
    end
    existsAndDefault('template_name', [timestamp '_template.mat']);
    
    data = load(fullfile(srcdir, register, template_name));
    w = data.final_template.w;
    zrange = [min(w(zi)) max(w(zi))];
catch
    zrange = [-Inf Inf];
end
% 
% min_ctrl_pts = 3;
% ctrl_pt_spacing = 20;
% 
% load(fullfile(srcdir, register, [timestamp '_template.mat']), 'final_template');
% 
% npts = max(ceil([diff(final_template.u([1 end])) diff(final_template.v([1 end])) diff(final_template.w([1 end]))]/ctrl_pt_spacing),min_ctrl_pts)
% tic
% ic = IntensityCorrectorBSpline(final_template.F,npts);
% toc

if ~exist(fullfile(srcdir, register, [timestamp '_intensity_corrections.mat']), 'file')
    corrections = IntensityCorrectorBSpline.correct_directory (srcdir, register, saveresults);
else
    corrections = [];
end


projections = xytproject(fullfile(srcdir, register), [],  zrange,corrections, true, appendz); %corrections = [] ->load corrections from disk

