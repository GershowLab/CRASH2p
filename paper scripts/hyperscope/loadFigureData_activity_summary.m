function [data,hassleap] = loadFigureData_activity_summary (srcdir, register, sleap_path, abortifnosleap) %fourier_map_name, vr_name, template_name)
%srcdir = 'R:\hyperscope registered data\A27h_gCamp7f_mCherry\tracking\20201204151516';
%register = 'registered';
[~,timestamp,~] = fileparts(srcdir);

existsAndDefault('vr_name', [timestamp '_vr.mat']);
existsAndDefault('abortifnosleap', false);
data = load(fullfile(srcdir,vr_name));
data.vr = data.vr.changeFilename(srcdir);
data.tsr = data.vr.tsr.conditionData();
data.tsr = data.tsr.addBehaviorVideo();
try
    data.tsr = data.tsr.addSleapResult(sleap_path);
    hassleap = true;
catch
     disp('did not load sleap results');
     hassleap = false;
     if (abortifnosleap)
         return;
     end
end
data.tsr = data.tsr.findCounterMovements();


data.vr.vrfa.vrfc.unload;

try
    s = load(fullfile(srcdir,register,'voi_selections.mat'));
    data = mergeStructs(data,s);
    try
        s = load(fullfile(srcdir, register, 'voi_activity.mat'));
    catch
        try
            s = load(fullfile(srcdir, register, 'voi_activity_old.mat'));
        catch
            d = dir(fullfile(srcdir, register, 'voi_activity*.mat'));
            if (length(d) == 1)
                s = load(fullfile(d.folder, d.name));
            else
                error ('can not find unambiguous voi activity file');
            end
        end
    end
    data = mergeStructs(data,s);
catch
end

existsAndDefault('template_name', [timestamp '_template.mat']);

s = load(fullfile(srcdir, register, template_name));
data = mergeStructs(data,s);

try
    s = load(fullfile(srcdir, register, [timestamp '_xytprojections.mat']));
    data = mergeStructs(data,s);
catch
    try
        d = dir(fullfile(srcdir, register, [timestamp '*_xytprojections.mat']));
        s = load(fullfile(d(1).folder, d(1).name));
        data = mergeStructs(data, s);
    catch
        disp ('xyt projections not loaded');
    end
end
try 
    d = dir(fullfile(srcdir, register, [timestamp '*_movie_projections*.mat']));
    for j = 1:length(d)
        s = load(fullfile(d(j).folder, d(j).name));
        fn = fieldnames(s);
        for k = 1:length(fn)
            data.mov_projections(j).(fn{k}) = s.(fn{k});
        end
    end
    for j = 1:length(d)
        data.mov_projections(j).fname = d(j).name;
    end
catch
    disp ('mov projections not loaded');
end


%fourier map and vra not used yet
% try
%     s = load(fullfile(srcdir, register, fourier_map_name));
%     data = mergeStructs(data,s);
% catch
% end
% 
% s = load(fullfile(srcdir, register, sprintf('%s_vra.mat',timestamp)));
% data = mergeStructs(data,s);

%data.activity = activity_vs_phase_in_voi(data.activity, data.tsr, [min(data.fmap.tx) max(data.fmap.tx)]);
try
    data.activity = prep_activity_for_analysis(data.activity);
catch
end
data.templateim = data.final_template.F;

end

function a = mergeStructs(a,b)
    fn = fieldnames(b);
    for i = 1:length(fn)
        a.(fn{i}) = b.(fn{i});
    end
end