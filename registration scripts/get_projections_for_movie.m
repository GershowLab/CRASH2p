function projections = get_projections_for_movie(srcdir,register,frameinds, ic, p, opp)
% default projection settings: l1spline smoothing, mean projection

%Q: should there be some kind of mask for mean projection? 

%progress bar - https://www.mathworks.com/help/parallel-computing/parallel.pool.dataqueue.html

ts0 = tic;

op.xi = [];
op.yi = [];
op.zi = [];
op.sigma = [1 1 1 1.3];
op.calc_mip = false;
op.maxcorrection = 10;
op.mindwell = 2.5e-8;
op.savefile = true;
op.savename = '';

op.xytproject = false;
op.nonrigid = true;
op.intensitycorrect = [];

if (nargin <= 0)
    projections = op;
    return;
end
[~,timestamp,~] = fileparts(srcdir);

try
    fn = fieldnames(opp);
    for j = 1:length(fn)
        op.(fn{j}) = opp.(fn{j});
    end
catch
end
if (isempty(op.intensitycorrect))
    op.intensitycorrect = op.nonrigid;
end
if (isempty(op.savefile))
    op.savefile = ~isempty(op.savename);
end
op.xytproject = op.xytproject && op.savefile;

if ~isfolder(fullfile(srcdir,register,'registered_frames'))
    fprintf('no frame saved\n');
    return;
end

if (nargin < 3 || isempty(frameinds))
    framepaths = dir(fullfile(srcdir,register,'registered_frames',[timestamp '_frame*.mat']));
    nframes = length(framepaths);
    frameinds = NaN(1,nframes);
    for i=1:nframes
        frameinds(i) = sscanf(framepaths(i).name,[timestamp '_frame%d.mat']);
    end
    [frameinds,sortinds] = sort(frameinds);
    framepaths = framepaths(sortinds);
else
    for i = 1:length(frameinds)
        framepaths(i) = dir(fullfile(srcdir,register,'registered_frames',sprintf('%s_frame%d.mat',timestamp,frameinds(i))));
    end
    nframes = length(framepaths);
end

try
    startframepath = dir(fullfile(srcdir,register,'registered_frames',sprintf('%s_frame%d.mat',timestamp,frameinds(1)-1)));
    if (isempty(startframepath))
        startframepath = framepaths(1);
    end
catch
    startframepath = framepaths(1);
end

try
    endframepath = dir(fullfile(srcdir,register,'registered_frames',sprintf('%s_frame%d.mat',timestamp,frameinds(end)+1)));
    if (isempty(endframepath))
        endframepath = framepaths(end);
    end
catch
    endframepath = framepaths(end);
end


if (nargin < 5 || isempty(p) && op.intensitycorrect)
    load(fullfile(srcdir,register,sprintf('%s_intensity_corrections.mat',timestamp)));
    ic = corrections.ic;
    p = corrections.p;
    corrections_inds = interp1(corrections.frame_range(1):corrections.frame_range(end), 1:size(corrections.p,1), frameinds, 'nearest','extrap');
    start_correction_ind = interp1(corrections.frame_range(1):corrections.frame_range(end), 1:size(corrections.p,1), frameinds(1)-1, 'nearest','extrap');
    end_correction_ind = interp1(corrections.frame_range(1):corrections.frame_range(end), 1:size(corrections.p,1), frameinds(1)+1, 'nearest','extrap');
    
else
    corrections_inds = 1:length(frameinds);
    start_correction_ind = 1;
    end_correction_ind = length(frameinds);
end
    
ts1 = tic;

frame = matfile(fullfile(startframepath.folder,startframepath.name));

xaxis = frame.xaxis;
yaxis = frame.yaxis;
zaxis = frame.zaxis;


if (op.nonrigid)
    rc_field = 'non_rigid_red_counts';
    rd_field = 'non_rigid_red_dwell';
    gc_field = 'non_rigid_green_counts';
    gd_field = 'non_rigid_green_dwell';
else
    rc_field = 'rigid_red_counts';
    rd_field = 'rigid_red_dwell';
    gc_field = 'rigid_green_counts';
    gd_field = 'rigid_green_dwell';
end

sz = size(frame.(rc_field),1:3);
nfine = size(frame.(rc_field), 4);

if (isempty(op.xi))
    op.xi = true([1 sz(1)]);
end
if (isempty(op.yi))
    op.yi = true([1 sz(2)]);
end
if (isempty(op.zi))
    op.zi = true([1 sz(3)]);
end


frametimes = zeros([diff(frameinds([1 end]))+1 nfine]);
red_xy = zeros([sz(2) sz(1) nframes nfine]); %transposed for you!
red_xz = zeros([sz(3) sz(1) nframes nfine]);
red_zy = zeros([sz(2) sz(3) nframes nfine]);
green_xy = red_xy; green_xz = red_xz; green_zy = red_zy;
if (op.calc_mip)
    green_xy_mip = red_xy; green_xz_mip = red_xz; green_zy_mip = red_zy;
    red_xy_mip = red_xy; red_xz_mip = red_xz; red_zy_mip = red_zy;
end
% if (op.xytproject)
%     projections_rc =  zeros([sz([1 2]) nframes nfine]);
%     projections_gc = projections_rc;
%     projections_rd = projections_rc;
%     projections_gd = projections_rc;
% end

ts1 = tic;


framepaths = [framepaths(:); endframepath];
corrections_inds = [corrections_inds(:); end_correction_ind];
if (op.intensitycorrect)    
    alpha = ic.calculateAlpha(p(start_correction_ind,:), op.maxcorrection);
else
    alpha = 1;
end
prev_rc = frame.(rc_field)./sqrt(alpha);
prev_rd = frame.(rd_field).*sqrt(alpha);
prev_gc = frame.(gc_field)./sqrt(alpha);
prev_gd = frame.(gd_field).*sqrt(alpha);

frame = matfile(fullfile(framepaths(1).folder,framepaths(1).name));
alpha = ic.calculateAlpha(p(corrections_inds(1),:), op.maxcorrection);
this_rc = frame.(rc_field)./sqrt(alpha);
this_rd = frame.(rd_field).*sqrt(alpha);
this_gc = frame.(gc_field)./sqrt(alpha);
this_gd = frame.(gd_field).*sqrt(alpha);


xi = op.xi;
yi = op.yi;
zi = op.zi;

if (isempty(op.savename))
    if (islogical(xi))
        x1 = find(xi,1,"first");
        x2 = find(xi,1,"last");
    else
        x1 = min(xi);
        x2 = max(xi);
    end
    if (islogical(yi))
        y1 = find(yi,1,"first");
        y2 = find(yi,1,"last");
    else
        y1 = min(yi);
        y2 = max(yi);
    end
    if (islogical(zi))
        z1 = find(zi,1,"first");
        z2 = find(zi,1,"last");
    else
        z1 = min(zi);
        z2 = max(zi);
    end
    if ~op.nonrigid
        timestamp = [timestamp '_rigid'];
    end
    op.savename = sprintf('%s_movie_projections_x_%d_%d_y_%d_%d_z_%d_%d.mat',timestamp, x1, x2, y1, y2, z1, z2);
end

for i=1:nframes
    frametimes(i,:) = conv(frame.time_edges, [.5 .5], 'valid');

    frame = matfile(fullfile(framepaths(i+1).folder,framepaths(i+1).name));
    if (op.intensitycorrect)
        alpha = ic.calculateAlpha(p(corrections_inds(i+1),:), op.maxcorrection);
    else
        alpha = 1;
    end
    next_rc = frame.(rc_field)./sqrt(alpha);
    next_rd = frame.(rd_field).*sqrt(alpha);
    next_gc = frame.(gc_field)./sqrt(alpha);
    next_gd = frame.(gd_field).*sqrt(alpha);

    rc = imblur(cat(4, prev_rc, this_rc, next_rc), op.sigma);
    gc = imblur(cat(4, prev_gc, this_gc, next_gc), op.sigma);
    rd = imblur(cat(4, prev_rd, this_rd, next_rd), op.sigma);
    gd = imblur(cat(4, prev_gd, this_gd, next_gd), op.sigma);

    rc = rc(:,:,:,nfine+(1:nfine));
    gc = gc(:,:,:,nfine+(1:nfine));
    rd = rd(:,:,:,nfine+(1:nfine));
    gd = gd(:,:,:,nfine+(1:nfine));

    valid = isfinite(rd) & isfinite(gd) & rd > op.mindwell & gd > op.mindwell;

    rc(~valid) = 0;
    rd(~valid) = 0;
    gc(~valid) = 0;
    gd(~valid) = 0;

    red_xy(:,:,i,:) = permute(sum(rc(:,:,zi,:),3)./ sum(rd(:,:,zi,:),3), [2 1 4 3]);
    red_xz(:,:,i,:) = permute(sum(rc(:,yi,:,:),2)./ sum(rd(:,yi,:,:),2), [3 1 4 2]);
    red_zy(:,:,i,:) = squeeze(sum(rc(xi,:,:,:),1)./ sum(rd(xi,:,:,:),1));
    green_xy(:,:,i,:) = permute(sum(gc(:,:,zi,:),3)./ sum(gd(:,:,zi,:),3), [2 1 4 3]);
    green_xz(:,:,i,:) = permute(sum(gc(:,yi,:,:),2)./ sum(gd(:,yi,:,:),2),[3 1 4 2]);
    green_zy(:,:,i,:) = squeeze(sum(gc(xi,:,:,:),1)./ sum(gd(xi,:,:,:),1));
    if (op.calc_mip)
        rv = rc./rd;
        rv(rd < op.mindwell) = 0;
        red_xy_mip(:,:,i,:) = permute(squeeze(max(rv(:,:,zi,:),[],3,'omitnan')), [2 1 3]);
        red_xz_mip(:,:,i,:) =  permute(squeeze(max(rv(:,yi,:,:),[],2,'omitnan')), [2 1 3]);
        red_zy_mip(:,:,i,:) = squeeze(max(rv(xi,:,:,:),[],1,'omitnan'));

        gv = gc./gd;
        gv(gd < op.mindwell) = 0;
        green_xy_mip(:,:,i,:) = permute(squeeze(max(gv(:,:,zi,:),[],3,'omitnan')),[2 1 3]);
        green_xz_mip(:,:,i,:) =  permute(squeeze(max(gv(:,yi,:,:),[],2,'omitnan')),[2 1 3]);
        green_zy_mip(:,:,i,:) = squeeze(max(gv(xi,:,:,:),[],1,'omitnan'));
    end
   
%     if (op.xytproject)
%         projections_rc(:,:,j,:) = sum(this_rc(:,:,zi,:),3);
%         projections_rd(:,:,j,:) = sum(this_rd(:,:,zi,:),3);
%         projections_gc(:,:,j,:) = sum(this_gc(:,:,zi,:),3);
%         projections_gd(:,:,j,:) = sum(this_gd(:,:,zi,:),3);
%     end

    prev_rc = this_rc;
    prev_rd = this_rd;
    prev_gc = this_gc;
    prev_gd = this_gd;

    this_rc = next_rc;
    this_rd = next_rd; 
    this_gc = next_gc; 
    this_gd = next_gd;
    
end
if (op.calc_mip)
    vars = {'red_xy','red_xz','red_zy','green_xy','green_xz','green_zy',...
        'red_xy_mip','red_xz_mip','red_zy_mip','green_xy_mip','green_xz_mip','green_zy_mip',...
        'frameinds','frametimes','xaxis','yaxis','zaxis','op'};
else
    vars = {'red_xy','red_xz','red_zy','green_xy','green_xz','green_zy',...
        'frameinds','frametimes','xaxis','yaxis','zaxis','op'};
end



if (op.savefile)    
    fpath = fullfile(srcdir,register,op.savename);
    save(fpath,vars{:},'-v7.3');
    projections.fpath = fpath;
end

% if(op.xytproject)
%     projections.frame_range = [min(frameinds) max(frameinds)];
%     projections.zrange = zaxis(zi);
%     projections.xaxis = xaxis;
%     projections.yaxis = yaxis;
%     projections.tx = zeros([1 nfine*nframes]);
%     projections.rc = zeros([sz([1 2]) length(projections.tx)]);
%     projections.gc = projections.rc;
%     projections.rd = projections.rc;
%     projections.gd = projections.rc;
% 
%     
%     for j = 1:length(frameinds)
%         tind = (1:nfine) + (j-1)*nfine;
%         projections.tx(tind) = frametimes(j,:);
%         projections.rc(:,:,tind) = projections_rc(:,:,j,:);
%         projections.gc(:,:,tind) = projections_gc(:,:,j,:);
%         projections.rd(:,:,tind) = projections_rd(:,:,j,:);
%         projections.gd(:,:,tind) = projections_gd(:,:,j,:);
%     end
%     fname = fullfile(srcdir,register, [timestamp sprintf('_z%dto%d_xytprojections.mat', find(zi, 1 ),find(zi, 1, 'last' ))]);
%     
%     dt = whos('projections');
% 
%     if dt.bytes < 2e+9
%         save (fname, 'projections');
%     else
%         % added for file sizes greater that 2GB
%         save (fname, 'projections', '-v7.3');
%     end
%     
%     clear projections projections_rc projections_gc projections_rd projections_gd
% end 


if (nargout > 0)
    for j = 1:length(vars)
        projections.(vars{j}) = eval(vars{j});
    end
end






fprintf('...all done! et %.1fs\n',toc(ts0));

end


