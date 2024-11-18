function activity = activity_in_vois(framestub, frame_range, voisets, corrections, dofine)

use_ic = existsAndDefault('corrections',[]);
existsAndDefault('dofine', true);
maxcorr = 20;
frame = load(sprintf('%s%d.mat', framestub, frame_range(1)));

tx = zeros([1 diff(frame_range)+1]);
activity.totals.rc = tx;
activity.totals.rc_nr = tx;
activity.totals.gc = tx;
activity.totals.gc_nr = tx;
activity.totals.rd = tx;
activity.totals.rd_nr = tx;
activity.totals.gd = tx;
activity.totals.gd_nr = tx;


ntimes = (length(frame.time_edges)-1);
txfine = zeros(1,ntimes*length(tx));

% helper variables for intensity correction
frameinds = frame_range(1):frame_range(end);

if (use_ic)
    pinds = interp1(corrections.frame_range(1):corrections.frame_range(end),1:size(corrections.p,1), frameinds,'nearest','extrap');
else    
    sz = size(frame.non_rigid_red_counts);
    sz = sz(1:3);
    npatch = size(frame.intensity_patch);
    imGrid = {1:sz(1),1:sz(2),1:sz(3)};
    xx = linspace(1,sz(1),npatch(1));
    yy = linspace(1,sz(2),npatch(2));
    zz = linspace(1,sz(3),npatch(3));
    patchGrid = {xx,yy,zz};
end

if (~isfield(frame, 'non_rigid_red_counts'))
    disp('non_rigid_red_counts and related data required; use add_fine_count_and_dwell_to_frames (suggest backing up first)');
    dofine = false;
end

for j = 1:length(voisets)
    voisets(j).red_rate = zeros(length(voisets(j).vois), length(tx));
    voisets(j).green_rate = voisets(j).red_rate;
   
    if (dofine)
        voisets(j).rdfine = zeros(length(voisets(j).vois), length(txfine));
        voisets(j).gdfine = voisets(j).rdfine;
        voisets(j).rcfine = voisets(j).rdfine;
        voisets(j).gcfine = voisets(j).rdfine;
    end
end

%varsToLoad = {'non_rigid_red_counts','non_rigid_red_dwell',...
%              'non_rigid_green_counts','non_rigid_green_dwell',...
%              'intensity_patch','time_edges'};

% missingIntensity = false;
ts1 = tic;

if (use_ic)
    disp('using intensity corrector b spline');
else
    disp('using intensity corrector from frame');
end


for j = 1:length(frameinds)
    frame = load(sprintf('%s%d.mat', framestub, frameinds(j)));%, varsToLoad{:});

   
    


    try
        activity.totals.rc(j) = sum(frame.red_counts,'all');
        activity.totals.gc(j) = sum(frame.green_counts,'all');
        activity.totals.rd(j) = sum(frame.red_dwell,'all');
        activity.totals.gd(j) = sum(frame.green_dwell,'all');
    catch
        % July 2023 frame data format updating: add rigid_ prefix, data
        % unchanged
        activity.totals.rc(j) = sum(frame.rigid_red_counts,'all');
        activity.totals.gc(j) = sum(frame.rigid_green_counts,'all');
        activity.totals.rd(j) = sum(frame.rigid_red_dwell,'all');
        activity.totals.gd(j) = sum(frame.rigid_green_dwell,'all');
    end
    
    %same intensity correction applies to whole frame - calculate it once
    %and move forward
    if (use_ic)
        int_corr = corrections.ic.calculateAlpha(corrections.p(pinds(j),:),maxcorr);
    else
        fi = griddedInterpolant(patchGrid,frame.intensity_patch);
        int_corr = fi(imGrid);
    end
    int_corr4 = repmat(int_corr, [1 1 1 size(frame.non_rigid_red_dwell,4)]);

    

    activity.totals.rc_nr(j) = sum(frame.non_rigid_red_counts,'all');
    activity.totals.gc_nr(j) = sum(frame.non_rigid_green_counts,'all');
    activity.totals.rd_nr(j) = sum(frame.non_rigid_red_dwell,'all');
    activity.totals.gd_nr(j) = sum(frame.non_rigid_green_dwell,'all');
    


    % get non-rigid counts and dwell
    rc = squeeze(sum(frame.non_rigid_red_counts,4));
    rd = squeeze(sum(frame.non_rigid_red_dwell,4));
    gc = squeeze(sum(frame.non_rigid_green_counts,4));
    gd = squeeze(sum(frame.non_rigid_green_dwell,4));
    % get non-rigid rates
    rv = int_corr.*rc./rd;
    rv(~isfinite(rv)) = 0;
    gv = int_corr.*gc./gd;
    gv(~isfinite(gv)) = 0;
    % apply intensity correction
%     fp = frame.intensity_patch;
% 
% 
%     
%     rv = IntensityCorrector.getFittedFramePar(fp,rv,patchGrid,imGrid);
%     gv = IntensityCorrector.getFittedFramePar(fp,gv,patchGrid,imGrid);
%     
    if (dofine)
        % get non-rigid counts and dwell subframes
        rc_fine = frame.non_rigid_red_counts;
        rd_fine = frame.non_rigid_red_dwell;
        gc_fine = frame.non_rigid_green_counts;
        gd_fine = frame.non_rigid_green_dwell;

        %intensity correction can EITHER multiply counts or divide dwell
        rc_fine_ic = rc_fine.*int_corr4;
        rd_fine_ic = rd_fine./int_corr4; 
        gc_fine_ic = gc_fine.*int_corr4;
        gd_fine_ic = gd_fine./int_corr4;

    end
    tx(j) = mean(frame.time_edges);
    ind0 = (j-1)*ntimes;

    for m = 1:ntimes
        txfine(ind0+m) = mean(frame.time_edges([m m+1]));
    end
    
    for n = 1:length(voisets)
        for k = 1:length(voisets(n).vois)
            if (iscell(voisets(n).vois))
                mask = voisets(n).vois{k}.makeMask;
            else
                 mask = voisets(n).vois(k).makeMask;
            end
            voisets(n).red_rate(k,j) = mean(rv(mask),'all', 'omitnan');
            voisets(n).green_rate(k,j) = mean(gv(mask),'all', 'omitnan');
            if (dofine)

                for m = 1:ntimes
                    voisets(n).rcfine(k,ind0+m) = sum(rc_fine(:,:,:,m).*mask,'all');
                    voisets(n).rcfine_ic(k,ind0+m) = sum(rc_fine_ic(:,:,:,m).*mask,'all');
                    voisets(n).rdfine(k,ind0+m) = sum(rd_fine(:,:,:,m).*mask,'all');
                    voisets(n).rdfine_ic(k,ind0+m) = sum(rd_fine_ic(:,:,:,m).*mask,'all');
                    voisets(n).gcfine(k,ind0+m) = sum(gc_fine(:,:,:,m).*mask,'all');
                    voisets(n).gcfine_ic(k,ind0+m) = sum(gc_fine_ic(:,:,:,m).*mask,'all');

                    voisets(n).gdfine(k,ind0+m) = sum(gd_fine(:,:,:,m).*mask,'all');
                    voisets(n).gdfine_ic(k,ind0+m) = sum(gd_fine_ic(:,:,:,m).*mask,'all');

                end
            end
        end
    end
    if (mod(j-1,20) == 0)
        fprintf('%d/%d - %.1f min elapsed %.1f min remain\n', j, length(frameinds), toc(ts1)/60, toc(ts1)/60 * (length(frameinds)-j)/j);
    end
end

activity.tx = tx;
activity.txfine = txfine;
activity.voisets = voisets;

end