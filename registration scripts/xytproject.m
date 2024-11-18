function projections = xytproject(framestub, frame_range, zrange, corrections, saveInDir, appendzrange)
%projections = xytproject(framestub, frame_range, zrange, corrections, saveInDir, appendzrange)

existsAndDefault('saveInDir',true);
existsAndDefault('appendzrange', false);

if (isfolder(framestub))
    dirname = framestub;
    d = dir(fullfile(dirname, 'registered_frames', '*_frame*.mat'));
    if (nargin < 4 || ~isstruct(corrections) || ~isfield(corrections, 'ic') || ~isfield(corrections.p))
        dd = dir(fullfile(dirname, '*_intensity_corrections.mat'));
        if (length(dd) ~= 1)
            error ('correction struct not provided and there is not a unique intensity corrections in the directory');
        end
        load(fullfile(dd(1).folder, dd(1).name), 'corrections');
    end
else
    dirname = '';
    d = dir([framestub '*.mat']);
    saveInDir = false;

    if (nargin < 4 || ~isstruct(corrections) || ~isfield(corrections, 'ic') || ~isfield(corrections.p))
        error('xtyproject needs a corrections struct with fields ic and p');
    end
end

ss = strsplit(d(1).name, '_frame');
framestub = fullfile(d(1).folder, [ss{1}, '_frame']);

timestamp = ss{1};

ind = NaN(size(d));
for j = 1:length(d)
    try
        ss = strsplit(d(j).name,'_frame');
        ind(j) = sscanf(ss{2}, '%d.mat');
    catch
    end
end
fullframerange = [min(ind,[],'omitnan') max(ind,[],'omitnan')];

existsAndDefault('frame_range',fullframerange);





existsAndDefault('zrange', [-Inf Inf]);

frame = load(sprintf('%s%d.mat', framestub, frame_range(1)));
% 
% % helper variables for intensity correction\
sz = size(frame.non_rigid_red_counts,1:3);
% npatch = size(frame.intensity_patch);
% imGrid = {1:sz(1),1:sz(2),1:sz(3)};
% xx = linspace(1,sz(1),npatch(1));
% yy = linspace(1,sz(2),npatch(2));
% zz = linspace(1,sz(3),npatch(3));
% patchGrid = {xx,yy,zz};


projections.frame_range = frame_range;
projections.zrange = zrange;
projections.xaxis = frame.xaxis;
projections.yaxis = frame.yaxis;
nfine = size(frame.non_rigid_green_counts, 4);

projections_tx = zeros([diff(frame_range)+1 nfine]);
projections_rc =  zeros([sz([1 2]) size(projections_tx,1) nfine]);

% projections_tx = zeros([1 nfine*(diff(frame_range)+1)]);
% projections_rc = zeros([sz([1 2]) length(projections.tx)]);
projections_gc = projections_rc;
projections_rd = projections_rc;
projections_gd = projections_rc;

projections_rc_uc = projections_rc;
projections_rd_uc = projections_rc;
zi = frame.zaxis >= zrange(1) & frame.zaxis <= zrange(end);

 if (appendzrange)
    fname = fullfile(dirname, [timestamp sprintf('_z%dto%d_xytprojections.mat', find(zi, 1 ),find(zi, 1, 'last' ))]);
 else
    fname = fullfile(dirname, [timestamp '_xytprojections.mat']);
 end

 if (exist(fname,'file'))
     [d,f,e] = fileparts(fname);
     f = [f '_backup_' datestr(now, 'yymmddHHMMSS')];
     movefile(fname, fullfile(d,[f e]));
 end


frameinds = frame_range(1):frame_range(end);
corrections_inds = interp1(corrections.frame_range(1):corrections.frame_range(end), 1:size(corrections.p,1), frameinds, 'nearest','extrap');
ic = corrections.ic;
corr_p = corrections.p(corrections_inds,:);

ts1 = tic;


D = parallel.pool.DataQueue;
h = waitbar(0, 'Calculating Projections ...');
afterEach(D, @nUpdateWaitbar);
N = length(frameinds); 
p = 1;
    function nUpdateWaitbar(~)
        waitbar(p/N, h);
        p = p + 1;
    end


parfor j = 1:N
   % tind = (1:4) + (j-1)*4;
    frame = matfile(sprintf('%s%d.mat', framestub, frameinds(j)));%, varsToLoad{:});
    projections_tx(j,:) = conv(frame.time_edges, [.5 .5], 'valid');
    alpha = ic.calculateAlpha(corr_p(j,:), 10)

%     %restrict correction to factor of 10
%     alpha(alpha < 0.1) = 0.1;
%     alpha(alpha > 10) = 10;

    %apply correction symmetrically to counts and dwll

    rc = frame.non_rigid_red_counts./sqrt(alpha); %auto expands dimension - thanks Matlab!
    rd = frame.non_rigid_red_dwell.*sqrt(alpha); 
    gc = frame.non_rigid_green_counts./sqrt(alpha);
    gd = frame.non_rigid_green_dwell.*sqrt(alpha);
    
    rc_uc = frame.non_rigid_red_counts;
    rd_uc = frame.non_rigid_red_dwell;
    
    projections_rc(:,:,j,:) = sum(rc(:,:,zi,:),3);
    projections_rd(:,:,j,:) = sum(rd(:,:,zi,:),3);
    projections_gc(:,:,j,:) = sum(gc(:,:,zi,:),3);
    projections_gd(:,:,j,:) = sum(gd(:,:,zi,:),3);
    projections_rc_uc(:,:,j,:) = sum(rc_uc(:,:,zi,:), 3);
    projections_rd_uc(:,:,j,:) = sum(rd_uc(:,:,zi,:), 3);
    

%     if (any(~zi))
%         projections.rc_allz(:,:,j) = sum(rc,3);
%         projections.rd_allz(:,:,j) = sum(rd,3);
%         projections.gc_allz(:,:,j) = sum(gc,3);
%         projections.gd_allz(:,:,j) = sum(gd,3);
%     end
%     if (mod(j,100) == 0)
%         disp([num2str(j) '/' num2str(length(frameinds))]);
%         toc(ts1);
%     end
      send(D, j);
      frame = []; %may help with a memory leak
end    
projections.tx = zeros([1 nfine*(diff(frame_range)+1)]);
projections.rc = zeros([sz([1 2]) length(projections.tx)]);
projections.gc = projections.rc;
projections.rd = projections.rc;
projections.gd = projections.rc;
for j = 1:length(frameinds)
    tind = (1:nfine) + (j-1)*nfine;
    projections.tx(tind) = projections_tx(j,:);
    projections.rc(:,:,tind) = projections_rc(:,:,j,:);
    projections.gc(:,:,tind) = projections_gc(:,:,j,:);
    projections.rd(:,:,tind) = projections_rd(:,:,j,:);
    projections.gd(:,:,tind) = projections_gd(:,:,j,:);
    projections.rc_uc(:,:,tind) = projections_rc_uc(:,:,j,:);
    projections.rd_uc(:,:,tind) = projections_rd_uc(:,:,j,:);  
end

fprintf('saving xyt projections to %s, et %.1fs\n',fname,toc(ts1));

if (saveInDir)
    dt = whos('projections'); 
   
    if dt.bytes < 2e+9
        save (fname, 'projections');
    else
        % added for file sizes greater that 2GB
        save (fname, 'projections', '-v7.3');
    end
end


fprintf('...all done! et %.1fs\n',toc(ts1));


close(h)

end