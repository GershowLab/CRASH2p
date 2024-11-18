function activity = addActivityInVOIsToDirectory(dirname, voisetfilename, framestub)

existsAndDefault('voisetfilename', fullfile(dirname,'voi_selections.mat'));
existsAndDefault('framestub','');
load(voisetfilename,'voisets');

if (isempty(framestub))
    d = dir(fullfile(dirname, 'registered_frames', '*_frame*.mat'));
else
    d = dir([framestub '*.mat']);
end
ss = strsplit(d(1).name, '_frame');
framestub = fullfile(d(1).folder, [ss{1}, '_frame']);
ind = NaN(size(d));
for j = 1:length(d)
    try
        ss = strsplit(d(j).name,'_frame');
        ind(j) = sscanf(ss{2}, '%d.mat');
    catch
    end
end
framerange = [min(ind,[],'omitnan') max(ind,[],'omitnan')];

try
    d = dir(fullfile(dirname,'*_intensity_corrections.mat'));
    load(fullfile(d.folder, d.name), 'corrections')
catch
    corrections = [];
end
activity = activity_in_vois(framestub, framerange, voisets,corrections, true);
save (fullfile(dirname, 'voi_activity.mat'), 'activity');

end

