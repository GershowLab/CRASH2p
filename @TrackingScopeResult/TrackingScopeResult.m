classdef TrackingScopeResult
    %replaces TrackingMicroscopeResult for new data
    
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tracker;
        stage;
        event_log;
        neuron;
        com;
        behavior;
        video;
        dt = 1e-3;
        d_loglambda = [.01 .1]; %ratio diffusion constant, total fluorescence diffusion constant
        tx;
        rate_fit;
        hasVideo = false;
    end
    properties (Transient = true)
        videocube=[];
    end
    properties (Constant)
        FPGA_CLOCK = VRFileChunk.FPGA_CLOCK;
    end
    
    methods
        function tsr = TrackingScopeResult(varargin)
            if (nargin > 0)
                tic;
                tsr = tsr.loadFile(varargin{:});
                toc; disp ('loaded data');
                try 
                    tsr = tsr.conditionData();
                    tsr = tsr.processPhotonCounts();
                    tsr = tsr.correctRatioForBleaching();
                    toc; disp ('processed data and photon counts');
                catch me
                    disp (me.getReport());
                end
                % do not autoload video - it's slow
%                 try
%                     tsr = tsr.addVideo();
%                 catch me
%                     disp (me.getReport());
%                 end
            end
        end        

        function tsr = loadAnotations(tsr)
            try
                annotationFileName = tsr.behavior.annotationFileName;
                data = load(annotationFileName);
                annotation = data.annotation;
            catch
                [d,f,~] = fileparts(tsr.behavior.videoFileName);
                tsr.behavior.annotationFileName = fullfile(d,[f '_annotation.mat']);
                annotation = repmat(struct('frame',NaN,'gut',[NaN NaN],'spot', [NaN NaN], 'head', [NaN NaN]),0);
            end
            if (~isfield(tsr.behavior, 'annotation') || isempty(tsr.behavior.annotation))
                tsr.behavior.annotation = annotation;
            end
            frameloaded = [annotation.frame];
            currentframe = [tsr.behavior.annotation.frame];
            mergein = ~ismember(frameloaded,currentframe); %do not overwrite annotated frames in tsr
            tsr.behavior.annotation = [tsr.behavior.annotation, annotation(mergein)];
        end

        function tsr = annotateFrame(tsr, et)
          %  function tsr = annotateFrame(tsr, et)
          %  note that function asks for TIME STAMP, not frame number

            if length(et) > 1
                for j = 1:length(et)
                    tsr = tsr.annotateFrame(et);
                end
            end

            [im,framenum] = tsr.getBehaviorImage(et);
            figure(1); clf(1);
            pcolor(im); shading flat; colormap gray; axis equal;
            title ('click (in order) gut, spot, head')
            [x,y] = ginput(3);
            annotation.frame = framenum;
            annotation.gut = [x(1),y(1)];
            annotation.spot = [x(2), y(2)];
            annotation.head = [x(3), y(3)];
            if (~isfield(tsr.behavior, 'annotation') || isempty(tsr.behavior.annotation))
                tsr.behavior.annotation = annotation;
            end
            I = find([tsr.behavior.annotation.frame] == annotation.frame, 1, 'first');
            if isempty(I)
                tsr.behavior.annotation = [tsr.behavior.annotation annotation];
            else
                tsr.behavior.annotation(I) = annotation;
            end
        end

        function saveAnotations(tsr)
            annotation = tsr.behavior.annotation;
            try
                save(tsr.behavior.annotationFileName, 'annotation');
            catch
                [d,f,~] = fileparts(tsr.behavior.videoFileName);
                tsr.behavior.annotationFileName = fullfile(d,[f '_annotation.mat']);
                save(tsr.behavior.annotationFileName, 'annotation');
            end
        end

        function tsr = addSleapResult(tsr, sleapFilePath)
            existsAndDefault('sleapFilePath',[]);
            if isempty(sleapFilePath)
                [tsr,sleapFilePath] = tsr.parseAndSaveSleapH5Result2MatFile;
            else
                if ~isfile(sleapFilePath) || endsWith(sleapFilePath,'.analysis.h5')
                    h5FilePath = strrep(sleapFilePath,'.mat','.analysis.h5');
                    [tsr,sleapFilePath] = tsr.parseAndSaveSleapH5Result2MatFile(h5FilePath);
                end
            end
            load(sleapFilePath,'behav_point','behav_segment','behav_timer');
            tsr.behavior.sleap.filename = sleapFilePath;
            spotind = find(strcmpi({behav_point.name}, 'spot'),1,'first');
            if (isempty(spotind))
                warning ('expecting labeled spot in behav_point - exiting')
                return;
            end
            behavpixelsize = 7.4 / 2.5; %calibrated from sensor and magnifcation -- shouldn't be hardcoded but what the heck -- from aca640-90um pixel @ 2.5 x magnification
            
            % gut jump check: if spot-gut orientation differs too much from
            % local average (lowpassed by 2s), then set corresponding gut
            % info to NaN; doing this removes weird kinks in tsr.com.vfwd later
            k = find(strcmpi({behav_segment.name}, 'spot-gut'),1,'first');
            theta = angle(behav_segment(k).unit_vector(1,:)+1i*behav_segment(k).unit_vector(2,:));
            theta = fillmissing(theta,'previous','EndValues',NaN);
            dtheta = abs(theta-lowpass1D(theta,2*tsr.behavior.framerate_Hz));
            badinds = dtheta>=pi/2;
            gutind = find(strcmpi({behav_point.name}, 'gut'),1,'first');
            behav_point(gutind).pos_px(:,badinds) = NaN;
            behav_point(gutind).pos(:,badinds) = NaN;
            behav_point(gutind).vel(:,badinds) = NaN;
            behav_point(gutind).vel_abs(badinds) = NaN;
            behav_point(gutind).acc(:,badinds) = NaN;
            behav_point(gutind).acc_abs(badinds) = NaN;
            
            sl = behavpixelsize*behav_point(spotind).pos_px; 

            tsr.behavior.sleap.tx = behav_timer.et;
            loc_smooth = mean(tsr.tracker.loc_smooth(1:2,:,:),3);
            tl = interp1(tsr.tx, loc_smooth', tsr.behavior.sleap.tx)'; % use center of mass of all tracked neuron if multiple

            valid = all(isfinite(sl) & isfinite(tl));
            sl = sl(:,valid);
            tl = tl(:,valid);

            s0 = mean(sl,2,'omitnan');
            t0 = mean(tl,2,'omitnan');

            [u,s,v] = svd((sl-s0)*(tl-t0)');
            r = v*u'; %rotation might include a reflection, and that's OK

            at = zeros(3,3);
            at(1:2,1:2) = r;
            at(1:2,3) = t0 - r*s0;
            at(3,3) = 1;

            tsr.behavior.at_px2microns = at*[behavpixelsize 0 0; 0 behavpixelsize 0; 0 0 1];
            tsr.behavior.at_microns2px = inv(at*[behavpixelsize 0 0; 0 behavpixelsize 0; 0 0 1]);

            for j = 1:length(behav_point)
                try
                    x = behavpixelsize*behav_point(j).pos_px; x(3,:) = 1;
                    valid = all(isfinite(x),1);
                    x2 = interp1(tsr.behavior.sleap.tx(valid), (at*x(:,valid))', tsr.behavior.sleap.tx, 'linear', NaN)';
                    valid = all(isfinite(x2),1);
                    
                    x2 = interp1(tsr.behavior.sleap.tx(valid), x2(:,valid)', tsr.behavior.sleap.tx, 'nearest', 'extrap')';
                    %tsr.behavior.sleap.(behav_point(j).name) = at*x;
                    tsr.behavior.sleap.(behav_point(j).name) = x2(1:2,:);
                catch
                    fprintf('failed to add sleap %s labels\n',behav_point(j).name);
                    tsr.behavior.sleap.(behav_point(j).name) = NaN(2,length(tsr.behavior.sleap.tx));
                end
            end


        end
        
        function [tsr,sleapFilePath] = parseAndSaveSleapH5Result2MatFile (tsr, h5FilePath, overwrite)
            existsAndDefault('h5FilePath',[]);
            existsAndDefault('overwrite',false);
            if isempty(tsr.behavior)
                tsr = tsr.addBehaviorVideo;
            end
            % parse file paths
            if isempty(h5FilePath)
                fdir = fileparts(tsr.behavior.videoFileName);
                flist = dir(fullfile(fdir,'*.analysis.h5'));
                if isempty(flist)
                    error('no sleap .h5 file found! export from SLEAP GUI first');
                end
                h5FilePath = fullfile(flist(end).folder,flist(end).name);
            end
            sleapFilePath = strrep(h5FilePath,'.analysis.h5','.mat');
            if ~isfile(h5FilePath)
                error('%s not found! check file path or export from SLEAP GUI first',h5FilePath);
            end
            if isfile(sleapFilePath) && ~overwrite
                disp('sleap .mat file already saved!');
                fprintf('found .h5 file: %s\n',h5FilePath);
                fprintf('found .mat file: %s\n',sleapFilePath);
                return;
            end
            % get metadata
            bvr = VideoReader(tsr.behavior.videoFileName);
            w = bvr.Width;
            h = bvr.Height;
            et = tsr.behavior.et;
            nframes = numel(et);
            framerange = [1 length(et)];
            framestart = min(framerange);
            framestop = max(framerange);
            % get stage positions interped in behavior frame rate
            stagepos = interp1(tsr.tx,tsr.tracker.stageloc(1:2,:)',tsr.behavior.et)';
            behav_pxsize = 7.4/2.5*[1;1]; % camera: aca640-90um pixel @ 2.5x magnification (copied from vr.playMovie)
            % load and parse .h5 point data
            data = h5read(h5FilePath,'/tracks');
            if size(data,1)<nframes
                data = cat(1,data,NaN(nframes-size(data,1),size(data,2),size(data,3)));
            end
            nlabels = size(data,2);
            labelnames = deblank(h5read(h5FilePath,'/node_names'));
            behav_point = struct('name',[],'pos_px',NaN(2,nframes),'pos',NaN(2,nframes));
            behav_point = repmat(behav_point,1,nlabels);
            for j=1:nlabels
                behav_point(j).name = labelnames{j};
                pos_px = squeeze(data(framestart:framestop,j,:))';
                behav_point(j).pos_px(:,framestart:framestop) = pos_px;
                pos_px_flipped = [w-pos_px(1,:);pos_px(2,:)];
                pos_um = pos_px_flipped.*behav_pxsize+stagepos(:,framestart:framestop);
                pos_um = fillmissing(pos_um,'linear',2,'MaxGap',5);
                behav_point(j).pos(:,framestart:framestop) = pos_um;
                v = diff(behav_point(j).pos,1,2);
                v = [v v(:,end)]; % extrapolate last frame using previous frame
                behav_point(j).vel = v;
                behav_point(j).vel_abs = vecnorm(behav_point(j).vel,2,1);
                a = diff(behav_point(j).vel,1,2);
                a = [a a(:,end)];
                behav_point(j).acc = a;
                behav_point(j).acc_abs = vecnorm(behav_point(j).acc,2,1);
            end
            % calculate segment data
            behav_segment = struct('name',[],'unit_vector',NaN(2,nframes),'length',NaN(1,nframes));
            segpairs = [1,2;2,3;1,4;1,5;4,2;5,2];
            nsegs = size(segpairs,1);
            behav_segment = repmat(behav_segment,1,nsegs);
            for j=1:nsegs
                a = segpairs(j,1);
                b = segpairs(j,2);
                behav_segment(j).name = sprintf('%s-%s',behav_point(a).name,behav_point(b).name);
                x = behav_point(b).pos-behav_point(a).pos;
                behav_segment(j).length = vecnorm(x,2,1);
                behav_segment(j).unit_vector = x./vecnorm(x,2,1);
            end
            % get timer data
            behav_timer.behav_vid_time = (1:nframes)/tsr.behavior.framerate_Hz;
            behav_timer.et = et;
            behav_timer.has_label = isfinite(behav_point(1).pos_px(1,:));
            % save parsed results
            save(sleapFilePath,'behav_point','behav_segment','behav_timer');
            disp('parsed and saved SLEAP .h5 data to .mat file');
        end
        
        function tsr = loadFile (tsr, filename, filename2, filename3)
            if (nargin > 1 && ischar(filename))
                if (exist (filename, 'dir'))
                    d1 = dir(fullfile(filename, '*_tracker.txt'));
                    d2 = dir(fullfile(filename, '*_stage_enc.txt'));
                    d3 = dir(fullfile(filename, '*_event_log.txt'));
                    tsr = tsr.loadFile(fullfile(filename, d1.name),...
                        fullfile(filename, d2.name),...
                        fullfile(filename, d3.name));
                    return;
                end
                
                % load _tracker.txt
                raw = TrackingScopeResult.importDataToStruct(filename, true);
                while (any(diff(raw.timer) < 0))
                    ind = find(diff(raw.timer) < 0, 1, 'first');
                    raw.timer((ind+1):end) = raw.timer((ind+1):end) + 2^32;
                end
                
                ks = TrackingScopeResult.getKalmanScalings();
                fn = intersect(fieldnames(ks), fieldnames(raw));
                for j = 1:length(fn)
                    raw.(fn{j}) = raw.(fn{j})*ks.(fn{j});
                end
                tsr.tracker.filename = filename;
            end
            
            % load _stage.txt
            if (nargin > 2)
                if (contains(filename2, '_enc'))
                    if (nargin > 3)
                        enc = parseStageEncoder(filename2, filename3);
                    else
                        enc = parseStageEncoder(filename2);
                    end
                    tsr.stage.enc = enc;
                    if (~isempty(enc))
                        tsr.stage.filename = filename2;
                        tsr.stage.loc = enc.loc;
                        tsr.stage.ms_timer = enc.mstime/1000;
                    else
                        tsr.stage.loc = zeros([3 2]);
                        tsr.stage.ms_timer = [0 2^32];
                    end
                else
                    tsr.stage = TrackingScopeResult.importDataToStruct(filename2);
                    tsr.stage.filename = filename2;
                    tsr.stage.loc = [tsr.stage.x; tsr.stage.y; tsr.stage.z];
                    tsr.stage.ms_timer = tsr.stage.ms_timer/1000; % change to seconds
                end
            end
            
            % load _event_log.txt
            if (nargin > 3)
                try % event_log.txt can be empty
                    tsr.event_log = TrackingScopeResult.importDataToStruct(filename3);
                    tsr.event_log.filename = filename3;
                    tsr.event_log.timer = tsr.event_log.timer/1.6e8; % change to seconds
                    if isfield(tsr.event_log,'ms_timer') % some older data don't have the ms_timer column
                        tsr.event_log.ms_timer = tsr.event_log.ms_timer/1e3; % change to seconds
                    end
                    tsr.event_log.IRflash = tsr.event_log.code==1; % code 1: IR flashes in behavioral camera
                catch me
                    disp(me.getReport());
                end
            end
            
            % deal with invalid values in any raw fields (don't do this in
            % production)
            raw = TrackingScopeResult.removeInvalidValues(raw);
            
            raw = TrackingScopeResult.processStatus(raw);
            
            raw.nspots = ceil(max(raw.axis(raw.status.tracking))/3);
            raw.spotnum = ceil((raw.axis+1)/3);
           
            raw.com = raw.weighted_sum./(raw.photons_ch0 .* raw.status.feedback0 + raw.photons_ch1 .* raw.status.feedback1);
            for j = 1:raw.nspots           
                raw.spot(j).xinds = raw.status.tracking & (raw.axis == 0 + 3*(j-1));
                raw.spot(j).yinds = raw.status.tracking & (raw.axis == 1 + 3*(j-1));
                raw.spot(j).zinds = raw.status.tracking & (raw.axis == 2 + 3*(j-1));
            end
            
            %somehow lots of duplicate entries in data -- check this out
            

            
            tsr.tracker.raw = raw;
        end
            
        function tsr = applySmoother(tsr, gain_times_radius)
            %only works if vel = 0
            try
                if length(gain_times_radius) < 3
                    gain_times_radius = gain_times_radius(1) * [1 1 1];
                end
                g = gain_times_radius(mod(tsr.tracker.raw.axis,3)+1);
                pkk = tsr.tracker.raw.Pxx;
                np = (tsr.tracker.raw.photons_ch0 .* tsr.tracker.raw.status.feedback0 + tsr.tracker.raw.photons_ch1 .* tsr.tracker.raw.status.feedback1);
                r = g.^2 ./ np;
                r(np <= 0) = 1e6;
                if (any ((r < pkk) & tsr.tracker.raw.status.tracking))
                    warning  ([num2str(nnz(r < pkk)) ' places measurement error less than estimated error - this should never happen - maybe bad inputs to gain_times_raidus?']);
                end
                pm1 = pkk.*r ./ (r-pkk);
                cc = pkk./pm1;
                cc(pm1 <= 0) = 1;
                cc(~isfinite(cc)) = 0; %should never happen, but just effectively start smoothing over if it does
                tsr.tracker.raw.loc_smooth = tsr.tracker.raw.loc;
                for j = min(tsr.tracker.raw.axis):max(tsr.tracker.raw.axis)
                    inds = tsr.tracker.raw.axis == j;
                    if (isempty(inds))
                        continue;
                    end
                    x = tsr.tracker.raw.loc(inds);
                    xs = x;
                    c = cc(inds);
                    for k = (length(x)-1):-1:1
                        xs(k) = x(k) + c(k)*(xs(k+1) - x(k));
                    end
                    tsr.tracker.raw.loc_smooth(inds) = xs;
                end
                raw = tsr.tracker.raw;
                et = raw.timer;

                for j = 1:raw.nspots
                    xinds = raw.spot(j).xinds;
                    yinds = raw.spot(j).yinds;
                    zinds = raw.spot(j).zinds;
                    x = resample_to_t(et(xinds), raw.loc_smooth(xinds), tsr.tx);
                    y = resample_to_t(et(yinds), raw.loc_smooth(yinds), tsr.tx);
                    z = resample_to_t(et(zinds), raw.loc_smooth(zinds), tsr.tx);
                    tsr.tracker.loc_smooth(:,:,j) = [x;y;z];
                end
            catch
                disp('not kalman smoothing -- using simple lowpass on location');
                tsr.tracker.loc_smooth = NaN(size(tsr.tracker.loc));
                for j=1:size(tsr.tracker.loc,3) % multiple neurons
                    tsr.tracker.loc_smooth(:,:,j) = lowpass1D(tsr.tracker.loc(:,:,j), 1);
                end
            end

        end
        
        function tsr = conditionData (tsr)
            raw = tsr.tracker.raw;
            et = raw.timer;
            
            tsr.tx = min(et(raw.status.tracking)):tsr.dt:max(et(raw.status.tracking));
            
            for j = 1:raw.nspots
                xinds = raw.spot(j).xinds;
                yinds = raw.spot(j).yinds;
                zinds = raw.spot(j).zinds;
                x = resample_to_t(et(xinds), raw.loc(xinds), tsr.tx);
                y = resample_to_t(et(yinds), raw.loc(yinds), tsr.tx);
                z = resample_to_t(et(zinds), raw.loc(zinds), tsr.tx);
%                 if ~(isempty(x) || isempty(y) || isempty(z))
                    tsr.tracker.loc(:,:,j) = [x;y;z];
%                 else
%                     tsr.tracker.loc(:,:,j) = NaN(3,length(tsr.tx));
%                 end
                dx = resample_to_t(et(xinds), raw.com(xinds), tsr.tx);
                dy = resample_to_t(et(yinds), raw.com(yinds), tsr.tx);
                dz = resample_to_t(et(zinds), raw.com(zinds), tsr.tx);
%                 if ~(isempty(dx) || isempty(dy) || isempty(dz))
                    tsr.tracker.com(:,:,j) = [dx;dy;dz];
%                 else
%                     tsr.tracker.com(:,:,j) = NaN(3,length(tsr.tx));
%                 end
            end
            xinds = any(cat(1,raw.spot.xinds),1);
            yinds = any(cat(1,raw.spot.yinds),1);
            zinds = any(cat(1,raw.spot.zinds),1);
            if (isfield(raw, 'vel'))
                vel = raw.vel_com;
                vx = resample_to_t(et(xinds), vel(xinds), tsr.tx);
                vy = resample_to_t(et(yinds), vel(yinds), tsr.tx);
                vz = resample_to_t(et(zinds), vel(zinds), tsr.tx);
                tsr.tracker.vel = [vx;vy;vz];
            end
            pz = resample_to_t(et(zinds), raw.piezo_val(zinds), tsr.tx);
            tsr.tracker.piezo = pz;
            xinds = any(cat(1,raw.spot.xinds),1);
            si = find(diff([0 tsr.tracker.raw.status.tracking]) > 0);
            ei = find(diff([tsr.tracker.raw.status.tracking 0]) < 0);
            [~,I] = max(ei - si);
            tsr.tracker.startTime = tsr.tracker.raw.timer(si(I));
            tsr.tracker.endTime = tsr.tracker.raw.timer(ei(I));             
            
            ax0 = raw.axis == 0;
            tsr.tracker.tracking = logical(interp1(raw.timer(ax0), double(raw.status.tracking(ax0)), tsr.tx, 'next', false)) & logical(interp1(raw.timer(ax0), double(raw.status.tracking(ax0)), tsr.tx, 'previous', false));
%             [et_adj,inds,~] = unique(et);
%             tsr.tracker.tracking = logical(interp1(et_adj, double(raw.status.tracking(inds)), tsr.tx, 'next', false)) & logical(interp1(et_adj, double(raw.status.tracking(inds)), tsr.tx, 'previous', false));
%             % deal with repetitive entires in raw.timer <- shouldn't
%             % happen since timer is FPGA ticks*256; need to investigate
            
            mstimer = raw.ms_timer(ax0);
            etx0 = raw.timer(ax0);
            inds = (find(diff(mstimer) > 0));
            
            mstimer_adj = interp1(etx0(inds), mstimer(inds), tsr.tx, 'linear', 'extrap');
            valid = all(isfinite(tsr.stage.loc));
            slix = griddedInterpolant(tsr.stage.ms_timer(valid), tsr.stage.loc(1,valid), 'linear', 'nearest');
            sliy = griddedInterpolant(tsr.stage.ms_timer(valid), tsr.stage.loc(2,valid), 'linear', 'nearest');
            sliz = griddedInterpolant(tsr.stage.ms_timer(valid), tsr.stage.loc(3,valid), 'linear', 'nearest');
            
            tsr.tracker.stageloc = -1.0*[slix(mstimer_adj);sliy(mstimer_adj);sliz(mstimer_adj)]; %negative added 8/10/23 by MHG
            %tsr.tracker.stageloc = interp1 (tsr.stage.ms_timer, tsr.stage.loc', mstimer_adj)';
            for j = 1:raw.nspots
                tsr.neuron(j).loc = tsr.tracker.stageloc + squeeze(tsr.tracker.loc(:,:,j));
               % tsr.neuron(j).loc_filt = tsr.lpfilt(tsr.neuron(j).loc);
                tsr.neuron(j).vel = 1/tsr.dt * deriv(tsr.neuron(j).loc, 1);
               % tsr.neuron(j).vel_filt = 1/tsr.dt * deriv(tsr.neuron(j).loc_filt, 1);
                tsr.neuron(j).speed = sqrt(sum(tsr.neuron(j).vel.^2));
               % tsr.neuron(j).speed_filt = sqrt(sum(tsr.neuron(j).vel_filt.^2));
            end
            tsr.com.loc = tsr.neuron(1).loc;
            for j = 2:raw.nspots
                tsr.com.loc = tsr.com.loc + tsr.neuron(j).loc;
            end
            tsr.com.loc = tsr.com.loc/raw.nspots;
            tsr.com.vel = 1/tsr.dt * deriv(tsr.com.loc, 1);
            tsr.com.speed = sqrt(sum(tsr.com.vel.^2));
            
        end
        
        function tsr = extractTimeRange(tsr, trange)
            ts = length(tsr.tx);
            inds = tsr.tx >= min(trange - 2*tsr.dt) & tsr.tx <= max(trange + 2*tsr.dt);
            
            tsr.tracker.raw = [];
            tsr.event_log = [];
            tsr.stage = [];
            tsr.rate_fit = [];
            tsr.video = []; tsr.videocube = []; tsr.hasVideo = false;
            
            fm = {'tracker', 'neuron'};
            for k = 1:length(fm)
                fn = fieldnames(tsr.(fm{k}));
                for j = 1:length(fn)
                    temp = tsr.(fm{k}).(fn{j});
                    if strcmp(fm{k},'neuron')
                        if (size(temp,2) == ts)
                            n = length(tsr.(fm{k}));
                            for l=1:n
                                temp = tsr.(fm{k})(l).(fn{j});
                                tsr.(fm{k})(l).(fn{j}) = temp(:,inds);
                            end
                        end
                    else
                        if (size(temp,2) == ts)
                            tsr.(fm{k}).(fn{j}) = temp(:,inds);
                        end
                    end
                end
            end
            
            tsr.tx = tsr.tx(inds);
            
        end
        
        function tsr = processPhotonCounts (tsr, varargin)
            
            sampletime = 0.01;
            varargin = assignApplicable(varargin);
            
            for j = 1:tsr.tracker.raw.nspots
                
                ratefit = mleratefit(tsr, sampletime, j);
                tsr.neuron(j).redrate = interp1(ratefit.et,ratefit.lambda(:,2), tsr.tx, 'linear', 'extrap');
                tsr.neuron(j).redrate_eb = interp1(ratefit.et,ratefit.sigma_lambda(:,2), tsr.tx, 'linear', 'extrap');
                tsr.neuron(j).greenrate = interp1(ratefit.et,ratefit.lambda(:,1), tsr.tx, 'linear', 'extrap');
                tsr.neuron(j).greenrate_eb = interp1(ratefit.et,ratefit.sigma_lambda(:,1), tsr.tx, 'linear', 'extrap');
                tsr.neuron(j).ratio = interp1(ratefit.et,ratefit.ratio, tsr.tx, 'linear', 'extrap');
                tsr.neuron(j).ratio_eb = interp1(ratefit.et,ratefit.sigma_ratio, tsr.tx, 'linear', 'extrap');
                if (j == 1)
                    tsr.rate_fit = ratefit;
                else
                    tsr.rate_fit(j) = ratefit;
                end
            end
            
            function ratefit = mleratefit(tsr, sampletime, spotnum)
                raw = tsr.tracker.raw;
                txbase = (min(tsr.tx)-sampletime/2):sampletime:(max(tsr.tx) + sampletime/2);
                pinds = (raw.spot(spotnum).xinds & raw.status.tracking);
                
                %bin photons for faster processing
                et = raw.timer(pinds);
                rp = raw.photons_ch1(pinds);
                gp = raw.photons_ch0(pinds);
                [~,~,bin] = histcounts(et, txbase);
                valid = bin > 0 & bin < length(txbase);
                nr = accumarray(bin(valid)', rp(valid)', [length(txbase)-1 1]);
                ng = accumarray(bin(valid)', gp(valid)', [length(txbase)-1 1]);
                tedges = accumarray(bin(valid)', et(valid)', [length(txbase)-1 1], @max);
                valid = tedges > 0;
                nr = nr(valid);
                ng = ng(valid);
                tedges = [txbase(1); tedges(valid)];
                if (tedges(1) >= tedges(2))
                    tedges(1) = tedges(2) - sampletime;
                end
                
                %calculate rate function
                [ratio, sigma_ratio, lambda, sigma_lambda, thetaOut, cmat_thetaOut] = TrackingMicroscopeResult.MLERates(tedges, nr, ng, tsr.d_loglambda(1), tsr.d_loglambda(end));
                ratefit.et = tedges(2:end);
                ratefit.theta = thetaOut;
                ratefit.wtheta = cmat_thetaOut;
                ratefit.ratio = ratio;
                ratefit.sigma_ratio = sigma_ratio;
                ratefit.lambda = lambda;
                ratefit.sigma_lambda = sigma_lambda;
            end
            
        end
        
        function tsr = correctRatioForBleaching(tsr, tbuffer, basepts, varargin)
            %REVISED IN A POTENTIALLY HACKY WAY 4/8/2019 by MHG, CHECK to
            %make sure fits make sense
            % fits the ratio to an exponential in time
            % finds the longest continuous tracking range to do this
            % and clips tbuffer from each end of this range
            % inspired by backcor.m by Vincent Mazet
            % basepts gives an initial list of baseline points, useful for
            % tracks with lots of activity
            
            trange = [];
            lptime = 0.1; %seconds
            ydcutoff = 0.25; %percentile
            
            assignApplicable(varargin);       
            if (isempty(trange)) 
                trange = tsr.tx([1 end]);
                ti = trange(1);
                tf = trange(end);
            
                [val_i, idx_i] = min(abs(tsr.tx - ti));
                [val_f, idx_f] = min(abs(tsr.tx - tf));
                %if multiple bouts of tracking, pick the longest

                nt = ~tsr.tracker.tracking; 
                nt(idx_i) = true; nt(idx_f) = true;
                tnt = tsr.tx(nt);
                [~,I] = max(diff(tnt));
                tstart = tnt(I);
                tend = tnt(I+1);
                existsAndDefault('tbuffer', 5); 
                tbuffer = min(tbuffer, (tend-tstart)/4); %don't clip more than 50% of range
                tstart = tstart + tbuffer;
                tend = tend - tbuffer;
            else
                tstart = min(trange);
                tend = max(trange);
            end
             inds = tsr.tx > tstart & tsr.tx < tend;
            
             xdata = tsr.tx(inds).';
             for k = 1:length(tsr.neuron)
                ydata = tsr.neuron(k).ratio(inds).';
                op = optimoptions('fminunc','Algorithm','trust-region');
                op.GradObj = 'on';
                op.Display = 'off';

                
                 xd = xdata;
                 yd = ydata;
                 
                 for j = 1:3
                     %estimate standard deviation by subtracting off exponential
                     %fit and only using negative values
                     %then discard outlying regions
                      x0 = polyfit(xd, log(yd), 1);
                     yf = exp(polyval(x0, xd));
                     dy = yd - yf;
                     s = std(dy);
                    
                     vrange = find(abs(dy) < s, 1, 'first'):find(abs(dy) < s, 1, 'last');
                     % plot (xdata, ydata, xd(vrange),yd(vrange),xd, yf); title(num2str(j)); pause
                     xd = xd(vrange);
                     yd = yd(vrange);
                     s = sqrt(mean(dy(dy < 0).^2));
                     
                    
                 end
                 expfit1 = exp(polyval(x0,xdata));
                 expfitfull = exp(polyval(x0,tsr.tx));
                 xd = xdata;
                 yd = ydata./expfit1;
                 
                 
                 
                 if (~existsAndDefault('basepts', {}) || length(basepts) < k || (k > 1 && ~iscell(basepts) || isempty(basepts{k})))
                     %go to lptime interval;
                     dx = ceil(lptime/median(diff(xd)));
                     lpy = lowpass1D(ydata./expfit1, dx/2);
                     
                     xdd = xd(1:dx:end);
                     ydd = lpy(1:dx:end);
                     ydd = ordfilt2(ydd, 1, true(10/lptime,1)); %+/-5 second interval
                     maxyd = percentile(ydd,ydcutoff);
                     valid = ydd > 0 & ydd < maxyd;
                     x0 = polyfit(xdd(valid), log(ydd(valid)), 1);  
                     dy = yd - exp(polyval(x0,xd));
                     s = sqrt(mean(dy(dy < 0).^2));
                 else
                     if (~iscell(basepts))
                         basepts = {basepts};
                     end
                     if (size(basepts{k},1) > size(basepts{k},2))
                         basepts{k} = basepts{k}';
                     end
                         
                     x0 = polyfit(basepts{k}(1,:), log(basepts{k}(2,:)), 1);
                     dy = yd - exp(polyval(x0,xd));
                     s = sqrt(mean(dy(dy < 0).^2));
                 end
                 %fit to exponential using assymetric least squares cost
                 %function
                 cfun = @(u) assymTruncLSQ(u, 3*s);
                 [x,yff] = fminArbCostFunction(@expfun, x0, xd, yd, cfun, op);
%                  plot (xd, yd, xd,expfun(x0,xd), xd, yff); pause

                 %estimate standar deviation by subtracting off fit & only
                 %using negative values
                 dy = yd - yff;
                 s = sqrt(mean(dy(dy < 0).^2));
                 
                 

                 %fit to exponential using cost function that falls off after 2 
                 %positive standard deviation to 0 at 4 standard deviations
                 %and has 0 cost for anything 20 standard deviations
                 %outside range in either direction
                 cfun = @(u) assymFallOffLSQ(u, 2*s);
                 %cfun = @(u) assymTruncLSQ(u, 3*s);

                 [x,yf] = fminArbCostFunction(@expfun, x, xd, yd, cfun, op);
             
                 %repeat
                 dy = yd - yf;
                 s = sqrt(mean(dy(dy < 0).^2));
                 cfun = @(u) assymFallOffLSQ(u, s);
%                cfun = @(u) assymTruncLSQ(u, 3*s);

                 [x,yf] = fminArbCostFunction(@expfun, x, xd, yd, cfun, op);
             
                 
                tsr.neuron(k).ratio_baseline = expfun(x, tsr.tx).*expfitfull;
                tsr.neuron(k).ratio_div_baseline = tsr.neuron(k).ratio./tsr.neuron(k).ratio_baseline;
              % plot (xdata, ydata, xd, yf, xd, yff); pause
             end
             
            function [f, dfdx] = expfun(x, xdata)
                f = exp(x(1)*xdata + x(2));
                dfdx = [f, xdata.*f];
            end
            
            function [f, g] = assymFallOffLSQ (x, thresh)
                f = x.^2;
                g = 2*x;
                f(x > thresh & x < 2*thresh) = 2*thresh.^2 - thresh*x(x > thresh & x < 2*thresh);
                g(x > thresh & x < 2*thresh) = -thresh;
                f(x > 2*thresh) = 0;
                g(x > 2*thresh) = 0;
                f(x < 10*thresh) = 0;
                g(x < 10*thresh) = 0;
            end
             
            function [f, g] = assymTruncLSQ (x, thresh)
                f = x.^2;
                g = 2*x;
                f(x > thresh) = thresh.^2;
                g(x > thresh) = 0;
            end
            
            function [f, g] = assymHuber (x, thresh)
                f = x.^2;
                g = 2*x;
                f(x > thresh) = 2*thresh*x(x > thresh) - thresh.^2;
                g(x > thresh) = 2 *thresh;
            end
             
        end
        
        function tsr = addBehaviorVideo(tsr)
            % get all behav_*.txt filenames
            behavdir = fullfile(fileparts(tsr.tracker.filename),'behavior_video');
            flist = dir(fullfile(behavdir,'behav_*.txt')); % default sort order should already be by timestamp
            nvids = numel(flist);
            
            % detect all behavior camera event series in *_event_log.txt
            eventlogfname = strrep(tsr.tracker.filename,'_tracker.txt','_event_log.txt');
            A = importdata(eventlogfname);
            eventlog = struct(A.colheaders{1},[]);
            for i=1:numel(A.colheaders)
                eventlog = setfield(eventlog,A.colheaders{i},A.data(:,i)); %#ok<SFLD>
            end
            CAMERA_CODE = 1;
            inds1 = find(eventlog.code==CAMERA_CODE);
            framenum = eventlog.data(inds1);
            inds2 = intersect(find(diff(framenum)>=1),find(diff(framenum)<=5)); % allows up to 5 dropped frames in event_log.txt
            ei = getSegment(inds2);
            ei = ei+[0 1]; % use inds1(ei(1,*):ei(2,*)) to get actual eventlog inds
            
            % create vr.behavior
            for i=1:nvids
                % 1) get behav_*.avi file name
                behavlogfname = fullfile(flist(i).folder,flist(i).name);
                behav(i).videoFileName = strrep(behavlogfname,'.txt','.avi'); %#ok<AGROW>
                % 2) load corresponding behav_*.txt
                A = importdata(behavlogfname);
                % TEMP column headers are not tabbed, will correct later
                ch = split(A.textdata{1});
                ch = ch(~cellfun('isempty',ch));
                A.colheaders = ch';
                behavlog = struct(A.colheaders{1},[]);
                for j=1:numel(A.colheaders)
                    behavlog.(A.colheaders{j}) = A.data(:,j);
                end
                % removing frames with errors - MHG 7/23
                valid = behavlog.error_code == 0;
                for j=1:numel(A.colheaders)
                    behavlog.(A.colheaders{j}) = behavlog.(A.colheaders{j})(valid);
                end
                behav(i).videoLog = behavlog; %#ok<AGROW>
                % 3) extract cameventlog for this track
                eln = diff(reshape(eventlog.data(inds1(ei)),size(ei)),1,2); % framenum range as in eventlog
                bln = diff(behavlog.frame_number([1 end])); % framenum range as in behavlog
                ii = find(abs(eln-bln)/bln<0.05,1,'first'); % frame disagreement < 5% of behavlog
                % POTENTIAL PROBLEM: multiple segments with similar length,
                % unlikely to happen though; if it does happen, default
                % order should be the correct one
                if isempty(ii)
                    warning('faild to load behavior video: too much disagreement between event_log.txt and behav_log.txt!');
                    continue;
                end
                % matching allows for single dropped frame at beginning and
                % end, as well as multiple dropped frames in the middle
                inds = inds1(ei(ii,1):ei(ii,2));
                cameventlog.frame_number = eventlog.data(inds);
                cameventlog.ms_timer = eventlog.ms_timer(inds);
                elt = eventlog.timer(inds);
                while (any(diff(elt) < 0))
                    ind = find(diff(elt) < 0, 1, 'first');
                    elt(ind+1:end) = elt(ind+1:end) + 2^32;
                end
                cameventlog.ticks = elt;
                cameventlog.et = cameventlog.ticks/tsr.FPGA_CLOCK;
                behav(i).cameraEventLog = cameventlog; %#ok<AGROW>
                % 4) get video et and eti
                % event_log.txt and behav_TIMESTAMP.txt may not cover the same
                % range of frames; assuming .avi always agrees with
                % behav_TIMESTAMP.txt, make et length agree with behavlog
                if length(cameventlog.frame_number)~=length(behavlog.frame_number)
                    et = interp1(cameventlog.frame_number,cameventlog.et,behavlog.frame_number,'linear','extrap');
                else
                    et = cameventlog.et;
                end
                behav(i).et = et(:)'; %#ok<AGROW>
                behav(i).framerate_Hz = 1/median(diff(et)); %#ok<AGROW> %7/14/2019 MHG - changed mode to median, more stable in the event frame rate jumps around
                % make sure .avi and .txt frame number matches
                v = VideoReader(behav(i).videoFileName); %#ok<TNMLP>
                if v.NumberOfFrames~=behavlog.frames_written_to_avi(end) %#ok<VIDREAD>
                    warning('Error: behavior video .txt and .avi disagree on frame number!');
                    behav(i).valid = false; %#ok<AGROW>
                else
                    behav(i).valid = true; %#ok<AGROW>
                end
            end
            bold = tsr.behavior;
            tsr.behavior = behav';
            for j = 1:length(bold)
                fn = setdiff(fieldnames(bold), fieldnames(tsr.behavior(j)));
                for k = 1:length(fn)
                    tsr.behavior(j).(fn{k}) = bold(j).(fn{k});
                end
            end

            disp('added behavior video');
        end
        
        function [im,framenum, xaxis, yaxis] = getBehaviorImage(tsr,v,t,bn)
            %  [im,framenum, xaxis, yaxis] = getBehaviorImage(tsr,v,t,bn)
            %  im = image
            %  framenum = avi frame number
            %  xaxis,yaxis = location of image, offset by stage location;
            %
            if (nargin < 3)
                t = v;
                v = [];
            end
            if (nargin < 4)
                bn = 1;
            end
            if (isempty(v))
                v = VideoReader(tsr.behavior(bn).videoFileName);                
                [im,framenum, xaxis, yaxis] = getBehaviorImage(tsr,v,t,bn);
                return;
            end
            
            framenum = interp1(tsr.behavior(bn).et,1:length(tsr.behavior(bn).et),t,'nearest','extrap');
            if (~isfinite(framenum))
                framenum = 1;%interp1(1:length(vr.behavior.et), vr.behavior.et, t, 'nearest','extrap');
            end
            v.CurrentTime = (framenum-1)/v.FrameRate;
            im = rgb2gray(v.readFrame);
            sl = interp1(tsr.tx, tsr.tracker.stageloc(1:2,:)', t,'nearest','extrap')';
            try
                pxsize = sqrt(abs(det(tsr.behavior.at_px2microns(1:2,1:2))));
                xaxis = ((1:size(im,1)) - tsr.behavior.at_microns2px(1,3))*sign(tsr.behavior.at_px2microns(1,1))*pxsize + sl(1);
                yaxis = ((1:size(im,2)) - tsr.behavior.at_microns2px(2,3))*sign(tsr.behavior.at_px2microns(2,2))*pxsize + sl(2);
            catch
                pxsize = 7.4/2.5;
                xaxis = ((1:size(im,1)) - size(im,2)/2)*pxsize + sl(1);
                yaxis = -((1:size(im,2)) - size(im,1)/2)*pxsize + sl(2);
            end


        end
        
        

        function [ims, frametime, framenum] = getBehaviorSequence(tsr,trange, aligned, tref, v, bn)
            existsAndDefault('tref', mean(trange));
            existsAndDefault('bn',1);
            existsAndDefault('v', VideoReader(tsr.behavior(bn).videoFileName));
            existsAndDefault('aligned', true);

            framenum = find(tsr.behavior(bn).et >= min(trange) & tsr.behavior(bn).et < max(trange));
            frametime = tsr.behavior(bn).et(framenum);
            sl = interp1(tsr.tx, tsr.tracker.stageloc(1:2,:)', frametime)';
            slr = interp1(tsr.tx, tsr.tracker.stageloc(1:2,:)', tref)';
            offset = tsr.behavior.at_microns2px(1:2,1:2)*(sl-slr);
            im1 = tsr.getBehaviorImage(v,frametime(1),bn);
            %offset = offset([2 1],:);

            if (aligned)
                spot = tsr.behavior.at_microns2px(1:2,1:2)*tsr.behavior.sleap.spot(:,framenum) + offset;
                theta = atan2d(spot(2,end)-spot(2,1), spot(1,end)-spot(1,1));
            else
                theta = 0;
            end



            offset = [cosd(theta) sind(theta);-sind(theta) cosd(theta)]*offset;

            mino = min(floor(min(offset,[],2)),0);
            maxo = max(ceil(max(offset,[],2)),0);

            im1 = padarray(padarray(im1, abs(mino([2 1])), 0, 'pre'), maxo([2 1]), 0, 'post');
            imsize = size(im1);
            ims = zeros([imsize, length(frametime)]);
            for j = 1:length(frametime)
                im = tsr.getBehaviorImage(v,frametime(j),bn);    
                ims(:,:,j) = imtranslate(imrotate(padarray(padarray(im, abs(mino([2 1])), 0, 'pre'), maxo([2 1]), 0, 'post'),theta,'crop'),offset(:,j)');
                pcolor(ims(:,:,j)); colormap gray; shading flat; axis equal; drawnow();
            end

            if(aligned)
                bodywidthpx = 100;
                [~,I] = max(mean(ims,[2 3]));
                yrange = I + [-bodywidthpx bodywidthpx];
                yrange = min(yrange, size(ims,1));
                yrange = max(yrange, 1);
                ims = ims(yrange(1):yrange(end),:,:);
            end



        end

        function tsr = addVideo(tsr,varargin) % deprecated
            if (nargin>1 && ischar(varargin{1}) && isfile(varargin{1}))
                % if first argin is .avi filepath, load that video
                [~,~,ext] = fileparts(varargin{1});
                if strcmp(ext,'.avi')
                    fpath = varargin{1};
                end
            elseif (tsr.hasVideo && isfield(tsr.video,'filename') && ~isempty(tsr.video.filename))
                % if tsr is previously saved as .mat, videocube will be
                % purged; reload from video.filename
                % CAUTION: this will reload and realign, updating both
                % tsr.video and tsr.videocube
                fpath = tsr.video.filename;
            else
                % prompt user to select a .avi file
                fdir = fileparts(tsr.tracker.filename);
                [fname,fdir] = uigetfile(fullfile(fdir,'*.avi'),'select behavioral video');
                if fname==0
                    disp('no behavioral video selected');
                    return;
                end
                fpath = fullfile(fdir,fname);
            end
            % load video
            tsr.video.filename = fpath;
            tsr = tsr.loadVideo(fpath);
            % align video to tracker time (tsr.tx) by IR flashes
            tsr = tsr.alignVideo(varargin{:});
            %disp('added behavioral video');
        end
        
        function tsr = loadVideo(tsr,fpath) % deprecated
            try
                v = VideoReader(fpath);
                % want to preallocate videocube, but don't know how long the video is,
                % so load it once to get et only
                ind = 0;
                while (v.hasFrame())
                    ind = ind+1;
                    tsr.video.et(ind) = v.currentTime();
                    fr = v.readFrame();
                    if (size(fr,3)>1) % in case the video is saved as RGB
                        fr = rgb2gray(fr);
                    end
                    tsr.video.intensity(ind) = sum(fr(:));
                end
                % now we know how long it is, can load actual frames
                tsr.videocube = zeros(v.Height,v.Width,length(tsr.video.et));
                v.currentTime = tsr.video.et(1);
                ind = 0;
                while (v.hasFrame())
                    ind = ind+1;
                    fr = v.readFrame();
                    if (size(fr,3)>1)
                        fr = rgb2gray(fr);
                    end
                    tsr.videocube(:,:,ind) = fr;
                end
                tsr.hasVideo = true;
            catch me
                disp(me.getReport());
                tsr.hasVideo = false;
            end
        end
        
        function tsr = alignVideo(tsr,varargin) % deprecated
            smoofac = 50;
            minintensitydrop = 0.2;
            alignwith = 'start';
            verbose = 0;
            assignApplicable(varargin);
            % 1) identify IR flashes
            % get regularized intensity
            I = tsr.video.intensity';
            Ismoo = smooth(I,smoofac);
            I = abs(I-Ismoo);
            I = I/max(I)+eps;
            % get FPGA event log of flashes
            event_flash = tsr.event_log.IRflash;
            % tsr_t = tsr.event_log.timer-tsr.event_log.timer(1);
            event_t = tsr.event_log.timer;
            nmax = sum(event_flash); % <- max number of flashes; actual flashes in video could be less
            [~,flashinds] = findpeaks(I,'NPeaks',nmax,'MinPeakHeight',minintensitydrop);
            % 2) align video time to tsr.tx
            % interpolate video timing from event_log
            event_flashinds = find(event_flash==1);
            dn = length(event_flashinds)-length(flashinds);
            if dn>0 % video is shorter than tsr
                switch alignwith
                    case 'start' % assume video stops early
                        event_flashinds = event_flashinds(1:end-dn);
                    case 'end' % assume video starts late
                        event_flashinds = event_flashinds(dn+1:end);
                end
            end
            tsr.video.missedIRFlashes = dn;
            tsr.video.alignwith = alignwith;
            event_flasht = event_t(event_flashinds);
            p = polyfit(flashinds(:),event_flasht(:),1);
            vidtx = polyval(p,1:length(I));
            tsr.video.tx = vidtx;
            tsr.video.timewarp = range(tsr.video.tx)/range(tsr.video.et);
            % 3) resample neuron data and add to video
            for j=1:length(tsr.neuron)
                vv = isfinite(tsr.video.tx);
                tsr.video.greenrate(j,vv) = resample_to_t(tsr.tx,tsr.neuron(j).greenrate,tsr.video.tx(vv));
                tsr.video.redrate(j,vv) = resample_to_t(tsr.tx,tsr.neuron(j).redrate,tsr.video.tx(vv));
                tsr.video.ratio_div_baseline(j,vv) = resample_to_t(tsr.tx,tsr.neuron(j).ratio_div_baseline,tsr.video.tx(vv));
            end
            if verbose>0 % check visuals
                figure('Units','normalized','Position',[0.1 0.1 0.3 0.6]);
                subplot(3,1,1);
                plot(event_t,event_flash,'k'); hold on;
                plot(tsr.video.et-tsr.video.et(1),I,'r'); hold off;
                title('IR flashes'); legend('tracker time','video time'); xlabel('t (s)')
                subplot(3,1,2);
                plot(1:length(I),I,'k',flashinds,I(flashinds),'ro');
                title('identified IR flashes in video'); xlabel('video frame'); ylabel('normalized intensity');
                subplot(3,1,3);
                plot(event_t,event_flash,'k',vidtx,I,'r');
                title(['IR flashes (aligning with ' alignwith ' of experiment)']); legend('tracker time','video time (aligned)'); xlabel('t (s)');
            end
        end
        
        function getMultiNeuronPlots(tsr, varargin)
            srcdir = fileparts(tsr.tracker.filename);
            t = tsr.tx;
            trange = true(size(t));
            n_neuron = numel(tsr.neuron);
            dstdir = srcdir;
            savefigure = false;
            assignApplicable(varargin);

            ra = RatioAnalysis;

            for i=1:n_neuron
                ra.plot_ratio(i, n_neuron, srcdir, t, trange, tsr, 'savefigure', savefigure, 'dstdir', dstdir);
            end

            ra.plot_ntrajectory(srcdir, n_neuron, trange, tsr, 'savefigure', savefigure, 'dstdir', dstdir);
            ra.plot_nactivity(srcdir, n_neuron, t, trange, tsr, 'savefigure', savefigure, 'dstdir', dstdir);
            ra.rel_locs(srcdir, n_neuron, t, trange, tsr, 'savefigure', savefigure, 'dstdir', dstdir);
            ra.plot_piezo_stage(srcdir, t, trange, tsr, 'savefigure', savefigure, 'dstdir', dstdir);
        end
        
        function playMovie(tsr,varargin)
            % NAME-VALUE inputs:
            %    'neuron_id': default is all tracked neurons
            % 	 'trange': default is tsr.tracker.start/endTime
            % 	 'twidth': time window to display photon rates/ratio (unit:
            %              s); default is 10s before and after 
            %	 'framerate': default is behavior video frame rate (unit: Hz)
            %	 'savemovie',true/[false]: whether to save .mp4 movie to disk
            %	 'fpath': specify full .mp4 path; if savemovie=true and fpath
            %             is not specified, will automatically save to data
            %             folder
            %	 'bn': for data sets with more than one behavior video saved
            
            % TODO add overlaytracker option? (with inset for z)
            
            % default playback and saving settings
            neuron_id = 1:numel(tsr.neuron);
            trange = [tsr.tracker.startTime tsr.tracker.endTime];
            twidth = 10; % time window before and after for neuron activity display (unit: s)
            framerate = [];
            savemovie = false;
            fpath = []; %#ok<NASGU>
            bn = 1; % for when there's multiple behavior videos
            
            assignApplicable(varargin);
            
            % get display params
            if isempty(framerate)
                framerate = tsr.behavior(bn).framerate_Hz;
            end
            jwidth = round(twidth*framerate);
            
            % prepare video reader for behavior images
            v = VideoReader(tsr.behavior(bn).videoFileName);
            
            % prepare video writer if saving movie
            if savemovie
                if isempty(fpath) %#ok<UNRCH>
                    % save to data folder if fpath not specified
                    fpath = strrep(tsr.tracker.filename,'_tracker.txt','_tsr_movie.mp4');
                end
                vw = VideoWriter(fpath,'MPEG-4');
                vw.FrameRate = framerate;
                open(vw);
            end
            
            % prepare interped tracks and neuron activities
            tstart = trange(1);
            tstop = trange(2);
            t = tstart:(1/framerate):tstop;
            for i=neuron_id
                % get interped trajectory and plotting params
                loc = interp1(tsr.tx,lowpass1D(tsr.neuron(i).loc(1:2,:),median(diff(t))/3)',t)';
                valid = interp1(tsr.tx,double(tsr.tracker.tracking),t,'nearest');
                valid(isnan(valid)) = 0;
                loc(:,~valid) = NaN;
                xrange = range(loc(1,:));
                yrange = range(loc(2,:));
                traj_box = round([min(loc(1,:))-xrange*0.05 max(loc(1,:))+xrange*0.05 ...
                    min(loc(2,:))-yrange*0.05 max(loc(2,:))+yrange*0.05]);
                % get interped activity and plotting params
                greenrate = interp1(tsr.tx,tsr.neuron(i).greenrate,t);
                redrate = interp1(tsr.tx,tsr.neuron(i).redrate,t);
                ratio = interp1(tsr.tx,tsr.neuron(i).ratio_div_baseline,t);
                gr_ylim = [0 max(max(greenrate),max(redrate))*1.1];
                ratio_ylim = [0 max(ratio)*1.1];
                % gather into struct
                neuron_interped(i).loc = loc; %#ok<AGROW>
                neuron_interped(i).valid = valid; %#ok<AGROW>
                neuron_interped(i).traj_box = traj_box'; %#ok<AGROW> % transpose for easier concatenation later
                neuron_interped(i).greenrate = greenrate; %#ok<AGROW>
                neuron_interped(i).redrate = redrate; %#ok<AGROW>
                neuron_interped(i).ratio = ratio; %#ok<AGROW>
                neuron_interped(i).gr_ylim = gr_ylim; %#ok<AGROW>
                neuron_interped(i).ratio_ylim = ratio_ylim; %#ok<AGROW>
            end
            traj_box = [neuron_interped(:).traj_box]';
            tmp1 = min(traj_box,[],1);
            tmp2 = max(traj_box,[],1);
            traj_box = [tmp1(1) tmp2(2) tmp1(3) tmp2(4)];
            
            % prepare figure window and axes
            n = numel(neuron_id);
            set(0,'Units','pixels');
            W = get(0,'screensize');
            W = W(3);
            u = min(40,floor(W/(8*n+12)));
            fig = figure('Units','pixels','Position',[1.5 1.5 8*n+9 18]*u);
            ax_behav = axes('Units','pixels','Position',[1 9 7 7]*u);
            ax_traj = axes('Units','pixels','Position',[1 1 7 7]*u);
            for i=1:n
                ax_gr(i) = axes('Units','pixels','Position',[8*i+1 9 7 7]*u); %#ok<AGROW>
                ax_gr_inset(i) = axes('Units','pixels','Position',[8*i+5 14 3 2]*u); %#ok<AGROW>
                ax_ratio(i) = axes('Units','pixels','Position',[8*i+1 1 7 7]*u); %#ok<AGROW>
                ax_ratio_inset(i) = axes('Units','pixels','Position',[8*i+5 6 3 2]*u); %#ok<AGROW>
            end
            
            % prepare first frames and plot params
            axis(ax_behav,'equal','tight','ij');
            hold(ax_traj,'on');
            for i=1:n
                ii = neuron_id(i);
                % whole trajectory
                h(i).traj = plot(ax_traj,neuron_interped(ii).loc(1,:),neuron_interped(ii).loc(2,:),'DisplayName',['neuron ' num2str(ii)]); %#ok<AGROW>
                h(i).tracker = plot(ax_traj,neuron_interped(ii).loc(1,1),neuron_interped(ii).loc(2,1),'ro','LineWidth',1.5); %#ok<AGROW>
                % green and red rates
                plot(ax_gr(i),t,neuron_interped(ii).greenrate,'g'); hold(ax_gr(i),'on');
                plot(ax_gr(i),t,neuron_interped(ii).redrate,'r');
                h(i).gr = plot(ax_gr(i),[t(1) t(1)],neuron_interped(ii).gr_ylim,'k'); %#ok<AGROW>
                hold(ax_gr(i),'off');
                xl = [t(1)-twidth t(1)+twidth];
                yl = [0 max(max(neuron_interped(ii).greenrate(1:jwidth*2)),max(neuron_interped(ii).redrate(1:jwidth*2)))*1.1];
                set(ax_gr(i),'XLim',xl,'YLim',yl,'XTick',[]);
                title(ax_gr(i),['neuron ' num2str(ii)]);
                % green and red rates inset
                plot(ax_gr_inset(i),t,neuron_interped(ii).greenrate,'g'); hold(ax_gr_inset(i),'on');
                plot(ax_gr_inset(i),t,neuron_interped(ii).redrate,'r');
                h(i).gr_inset = plot(ax_gr_inset(i),[t(1) t(1)],neuron_interped(ii).gr_ylim,'k'); %#ok<AGROW>
                hold(ax_gr_inset(i),'off');
                set(ax_gr_inset(i),'XLim',trange,'YLim',neuron_interped(ii).gr_ylim,'XTick',[],'YTick',[]);
                % ratio/baseline
                plot(ax_ratio(i),t,neuron_interped(ii).ratio,'k'); hold(ax_ratio(i),'on');
                h(i).ratio = plot(ax_ratio(i),[t(1) t(1)],neuron_interped(ii).ratio_ylim,'k'); %#ok<AGROW>
                hold(ax_ratio(i),'off');
                yl = [0 max(neuron_interped(ii).ratio(1:jwidth*2))*1.1];
                set(ax_ratio(i),'XLim',xl,'YLim',yl,'XTick',[t(1)-twidth t(1) t(1)+twidth],...
                    'XTickLabel',{['-' sprintf('%0.1f',twidth) 's'],...
                    [sprintf('%0.1f',t(1)) 's'],['+' sprintf('%0.1f',twidth) 's']});
                % ratio/baseline inset
                plot(ax_ratio_inset(i),t,neuron_interped(ii).ratio,'k'); hold(ax_ratio_inset(i),'on');
                h(i).ratio_inset = plot(ax_ratio_inset(i),[t(1) t(1)],neuron_interped(ii).ratio_ylim,'k'); %#ok<AGROW>
                hold(ax_ratio_inset(i),'off');
                set(ax_ratio_inset(i),'XLim',trange,'YLim',neuron_interped(ii).ratio_ylim,'XTick',[],'YTick',[]);
            end
            hold(ax_traj,'off');
            set(ax_traj,'XLim',traj_box(1:2),'YLim',traj_box(3:4),'XTick',[],'YTick',[],'Box','on');
            axis(ax_traj,'equal'); legend(ax_traj,[h(:).traj],'Location','best');
            
            % play movie
            for j=1:length(t)
                if ~ishandle(fig)
                    return;
                end
                % update behavior
                if t(j)>=tsr.behavior(bn).et(1) && t(j)<=tsr.behavior(bn).et(end)
                    behavIm = tsr.getBehaviorImage(v,t(j),bn);
                    behavIm = flip(flip(behavIm,1),2);
                else
                    behavIm = zeros(v.Width,v.Height);
                end
                imagesc(ax_behav,behavIm); colormap(ax_behav,gray(256)); caxis(ax_behav,[0 Inf]);
                set(ax_behav,'XTick',[],'YTick',[]);
                title(ax_behav,['behavior camera (AVI time: ' datestr(v.CurrentTime/(24*60*60),'MM:SS') ')']); % assuming video <1hr
                % update track and activity for each neuron
                for i=1:n
                    ii = neuron_id(i);
                    % update XData for all saved plots
                    h(i).tracker.XData = neuron_interped(ii).loc(1,j);
                    h(i).tracker.YData = neuron_interped(ii).loc(2,j);
                    h(i).gr.XData = [t(j) t(j)];
                    h(i).gr_inset.XData = [t(j) t(j)];
                    h(i).ratio.XData = [t(j) t(j)];
                    h(i).ratio_inset.XData = [t(j) t(j)];
                    % update axes properties
                    xl = [t(j)-twidth t(j)+twidth];
                    xt = [t(j)-twidth t(j) t(j)+twidth];
                    xtl = {['-' sprintf('%0.1f',twidth) 's'],...
                        [sprintf('%0.1f',t(j)) 's'],...
                        ['+' sprintf('%0.1f',twidth) 's']};
                    if j<=jwidth
                        yl_gr = [min(min(neuron_interped(ii).greenrate(1:jwidth*2)),min(neuron_interped(ii).redrate(1:jwidth*2)))...
                            max(max(neuron_interped(ii).greenrate(1:jwidth*2)),max(neuron_interped(ii).redrate(1:jwidth*2)))];
                        yl_ratio = [min(neuron_interped(ii).ratio(1:jwidth*2))...
                            max(neuron_interped(ii).ratio(1:jwidth*2))];
                    elseif j>=length(t)-jwidth
                        yl_gr = [min(min(neuron_interped(ii).greenrate(end-jwidth*2:end)),min(neuron_interped(ii).redrate(end-jwidth*2:end)))...
                            max(max(neuron_interped(ii).greenrate(end-jwidth*2:end)),max(neuron_interped(ii).redrate(end-jwidth*2:end)))];
                        yl_ratio = [min(neuron_interped(ii).ratio(end-jwidth*2:end))...
                            max(neuron_interped(ii).ratio(end-jwidth*2:end))];
                    else
                        yl_gr = [min(min(neuron_interped(ii).greenrate(j-jwidth:j+jwidth)),min(neuron_interped(ii).redrate(j-jwidth:j+jwidth)))...
                            max(max(neuron_interped(ii).greenrate(j-jwidth:j+jwidth)),max(neuron_interped(ii).redrate(j-jwidth:j+jwidth)))];
                        yl_ratio = [min(neuron_interped(ii).ratio(j-jwidth:j+jwidth))...
                            max(neuron_interped(ii).ratio(j-jwidth:j+jwidth))];
                    end
                    set(ax_gr(i),'XLim',xl,'YLim',yl_gr);
                    set(ax_ratio(i),'XLim',xl,'YLim',yl_ratio,'XTick',xt,'XTickLabel',xtl);
                end
                % metadata
                str = {fileparts(tsr.tracker.filename);...
                    ['time: ' sprintf('%0.2f',t(j))]};
                text(ax_behav,'Position',[1.5 8.25]*u,'String',str,'Units','pixels',...
                    'HorizontalAlignment','left','VerticalAlignment','middle',...
                    'FontSize',12,'FontWeight','bold','Interpreter','none');
                drawnow;
                if savemovie
                    writeVideo(vw,getframe(fig)); %#ok<UNRCH>
                end
            end
            
            if savemovie
                close(vw); %#ok<UNRCH>
            end
        end
        
        function playMovie_sameplot(tsr, varargin)
            % NAME-VALUE inputs:
            %    'neuron_id': default is all tracked neurons
            % 	 'trange': default is tsr.tracker.start/endTime
            % 	 'twidth': time window to display photon rates/ratio (unit:
            %              s); default is 10s before and after 
            %	 'framerate': default is behavior video frame rate (unit: Hz)
            %	 'savemovie',true/[false]: whether to save .mp4 movie to disk
            %	 'fpath': specify full .mp4 path; if savemovie=true and fpath
            %             is not specified, will automatically save to data
            %             folder
            %	 'bn': for data sets with more than one behavior video saved
                        
            % default playback and saving settings
            neuron_id = 1:numel(tsr.neuron);
            trange = [tsr.tracker.startTime tsr.tracker.endTime];
            twidth = 10; % time window before and after for neuron activity display (unit: s)
            framerate = [];
            ratio_ylim = [0 5];
            savemovie = false;
            plottitle = fileparts(tsr.tracker.filename);
            fpath = []; %#ok<NASGU>
            bn = 1; % for when there's multiple behavior videos
            
            assignApplicable(varargin);
            
            % get display params
            if isempty(framerate)
                framerate = tsr.behavior(bn).framerate_Hz;
            end
%             jwidth = round(twidth*framerate);
            
            % prepare behavior image loading
            try
                v = VideoReader(tsr.behavior(bn).videoFileName);
                BOX.behav = [0 v.Width 0 v.Height];
            catch
                v = [];
                BOX.behav = [0 480 0 480];
            end
            
            % prepare tracker and activity data
            tstart = trange(1);
            tstop = trange(2);
            t = tstart:(1/framerate):tstop;
            for i=neuron_id
                % Get interped trajectory and plotting params 
                % Get the trajectory in the same frame rate as the video:
                loc = interp1(tsr.tx,lowpass1D(tsr.neuron(i).loc(1:2,:),median(diff(t))/3)',t)'; % 
                valid = interp1(tsr.tx,double(tsr.tracker.tracking),t,'nearest');
                valid(isnan(valid)) = 0;
                loc(:,~valid) = NaN;
                xrange = range(loc(1,:)); % Used to set the plot range of the trajectory in x 
                yrange = range(loc(2,:)); % Used to set the plot range of the trajectory in y
                traj_box = round([min(loc(1,:))-xrange*0.05 max(loc(1,:))+xrange*0.05 ...
                    min(loc(2,:))-yrange*0.05 max(loc(2,:))+yrange*0.05]);
                % get interped activity and plotting params
                greenrate = interp1(tsr.tx,tsr.neuron(i).greenrate,t);
                redrate = interp1(tsr.tx,tsr.neuron(i).redrate,t);
                ratio = interp1(tsr.tx,tsr.neuron(i).ratio_div_baseline,t);
                gr_ylim = [0 max(max(greenrate),max(redrate))*1.1];
%                 ratio_ylim = [0 max(ratio)*1.1];
%                 ratio_ylim = [0 5];
                % gather into struct
                neuron_interped(i).loc = loc; %#ok<AGROW>
                neuron_interped(i).valid = valid; %#ok<AGROW>
                neuron_interped(i).traj_box = traj_box'; %#ok<AGROW> % transpose for easier concatenation later
                neuron_interped(i).greenrate = greenrate; %#ok<AGROW>
                neuron_interped(i).redrate = redrate; %#ok<AGROW>
                neuron_interped(i).ratio = ratio; %#ok<AGROW>
                neuron_interped(i).gr_ylim = gr_ylim; %#ok<AGROW>
                neuron_interped(i).ratio_ylim = ratio_ylim; %#ok<AGROW>
            end
            traj_box = [neuron_interped(:).traj_box]';
            tmp1 = min(traj_box,[],1);
            tmp2 = max(traj_box,[],1);
            BOX.tracker = [tmp1(1) tmp2(2) tmp1(3) tmp2(4)];
            % adjust aspect ratio if too large
            xrange = diff(BOX.tracker(1:2));
            yrange = diff(BOX.tracker(3:4));
            if (xrange/yrange)<1/3
                xm = mean(BOX.tracker(1:2));
                xrange = yrange/3;
                BOX.tracker(1:2) = [xm-xrange/2 xm+xrange/2];
            elseif (xrange/yrange)>3
                ym = mean(BOX.tracker(3:4));
                yrange = xrange/3;
                BOX.tracker(3:4) = [ym-yrange/2 ym+yrange/2];
            end
            
            % prepare VideoWriter
            if savemovie
                if isempty(fpath) %#ok<UNRCH>
                    % save to data folder if fpath not specified
                    fpath = strrep(tsr.tracker.filename,'_tracker.txt','_video.mp4');
                end
                vw = VideoWriter(fpath,'MPEG-4');
                vw.FrameRate = framerate;
                open(vw);
            end
            
            % prepare figure window and axes
            n = numel(neuron_id);
            set(0,'Units','pixels');
            W = get(0,'screensize');
            W = W(3);
            % Condition to get the correct figure size
            if n <= 2 
                u = min(40,floor(W/(8*3+12))); % n in here has to be at least 3 to get the correct figure size
                FIG = figure('Units','pixels','Position',[1.5 1.5 8*(3-1)+9 18]*u);
            else
                u = min(40,floor(W/(8*n+12)));
                FIG = figure('Units','pixels','Position',[1.5 1.5 8*(n-1)+9 18]*u);
            end
            AX.behav = axes('Units','pixels','Position',[1 9 7 7]*u);  % Plot position of the behavior video 
            AX.tracker = axes('Units','pixels','Position',[9 9 7 7]*u);   % Plot position of the trajectory video
            AX.ratio = axes('Units','pixels','Position',[17 9 7 7]*u);
            AX.ratio_inset = axes('Units','pixels','Position',[21 14 3 2]*u);
            for i=1:n
                AX.gr(i) = axes('Units','pixels','Position',[8*i-7 1 7 7]*u); % subplot position for the green and red
                AX.gr_inset(i) = axes('Units','pixels','Position',[8*i+5-8 6 3 2]*u); % inset position (green and red)
            end
            
            % prepare scale bars
            % ======== behavior ========
            behavpixelsize = 7.4/2.5; % aca640-90um pixel @ 2.5x magnification (copied from vr.playMovie)
            bs = 0.25*range(BOX.behav(1:2)*behavpixelsize);
            bs = round(bs/500)*500; % round to nearest .5mm
            BAR.behav.str = [sprintf('%0.1f',bs/1000) ' mm'];
            bs = bs/behavpixelsize;
            x = 0.05*BOX.behav(1)+0.95*BOX.behav(2)-[bs 0];
            y = 0.05*BOX.behav(3)+0.95*BOX.behav(4)+[0 0];
            BAR.behav.x = x;
            BAR.behav.y = y;
            % ======== tracker ========
            bs = 0.25*range(BOX.tracker(1:2));
            % NOTE: depend on trackerbox setting, bs value can range widely; round
            % differently for different order of magnitudes
            if bs>500
                bs = floor(max(1000,bs)/500)*500; % min 1mm; round to nearest .5mm
                BAR.tracker.str = [sprintf('%0.1f',bs/1000) ' mm'];
            elseif bs>50
                bs = floor(max(100,bs)/50)*50; % min 0.1mm; round to nearest .05mm
                BAR.tracker.str = [sprintf('%0.2f',bs/1000) ' mm'];
            else
                bs = floor(max(10,bs)/5)*5; % min 10um; round to nearest 5um
                BAR.tracker.str = [sprintf('%d',bs) ' um'];
            end
            x = 0.05*BOX.tracker(1)+0.95*BOX.tracker(2)-[bs 0];
            y = 0.05*BOX.tracker(3)+0.95*BOX.tracker(4)+[0 0];
            BAR.tracker.x = x;
            BAR.tracker.y = y;

            % prepare first frames and plot cosmetics
            % ======== behavior ========
            axes(AX.behav);
            try
                im_behav = tsr.getBehaviorImage(v,max(t(1),tsr.behavior(bn).et(1)),bn);
                im_behav = flip(flip(im_behav,1),2);
            catch
                im_behav = zeros(v.Width,v.Height);
            end
            h_behav = imagesc(im_behav); colormap(gca,gray(256));
            axis equal ij; axis(BOX.behav); xticks([]); yticks([]); title('behavior camera');
            labeledBar(gca,BAR.behav.x,BAR.behav.y,BAR.behav.str,'below','w',{'lineWidth',4});
            % ======== tracker and activity ========
            axes(AX.tracker);
            hold(AX.tracker,'on');
            colors = viridis(n); % define colors for each neuron
            for i=1:n
                ii = neuron_id(i);
                % ======== tracker ========
                plot(AX.tracker, neuron_interped(ii).loc(1,:),neuron_interped(ii).loc(2,:),...
                    'DisplayName',['neuron ' num2str(ii)], 'Color', colors(i,:));
                h(i).tracker = plot(AX.tracker, neuron_interped(ii).loc(1,1),neuron_interped(ii).loc(2,1),'ro','LineWidth',1.5); %#ok<AGROW>
                % ======== green and red rates ========
                plot(AX.gr(i),t,neuron_interped(ii).greenrate,'g-'); hold(AX.gr(i),'on');
                plot(AX.gr(i),t,neuron_interped(ii).redrate,'r--');
                h(i).gr = xline(AX.gr(i),t(1),'k'); %#ok<AGROW>
                hold(AX.gr(i),'off');
                xl = [t(1)-twidth t(1)+twidth];
%                 yl = [0 max(max(neuron_interped(ii).greenrate(1:jwidth*2)),max(neuron_interped(ii).redrate(1:jwidth*2)))*1.1];
                set(AX.gr(i),'XLim',xl,'YLim',neuron_interped(ii).gr_ylim,'XTick',[],'Box','on');
                title(AX.gr(i),['neuron ' num2str(ii)]);
                % ======== green and red rates inset ========
                plot(AX.gr_inset(i),t,neuron_interped(ii).greenrate,'g'); hold(AX.gr_inset(i),'on');
                plot(AX.gr_inset(i),t,neuron_interped(ii).redrate,'r');
                h(i).gr_inset = xline(AX.gr_inset(i),t(1),'k'); %#ok<AGROW>
                hold(AX.gr_inset(i),'off');
                set(AX.gr_inset(i),'XLim',trange,'YLim',neuron_interped(ii).gr_ylim,'XTick',[],'YTick',[]);
                % ======== ratio/baseline ========
                hold(AX.ratio,'on');
                plot(AX.ratio,t,neuron_interped(ii).ratio,'Color',colors(i,:),'LineWidth',2); 
                h(i).ratio = xline(AX.ratio,t(1)); %#ok<AGROW>
                hold(AX.ratio,'off');
                xl = [t(1)-twidth t(1)+twidth];
                xt = [t(1)-twidth t(1) t(1)+twidth];
                xtl = {['-' sprintf('%0.1f',twidth) 's'],...
                    [sprintf('%0.1f',t(1)) 's'],...
                    ['+' sprintf('%0.1f',twidth) 's']};
%                 yl = [0 max(neuron_interped(ii).ratio(1:jwidth*2))*1.1];
                yl = ratio_ylim;
                set(AX.ratio,'XLim',xl,'YLim',yl,'XTick',xt,'XTickLabel',xtl,'Box','on');
                title(AX.ratio,'Ratio: R/R_0');
                % ======== ratio/baseline inset ========
                hold(AX.ratio_inset(1),'on');
                plot(AX.ratio_inset,t,neuron_interped(ii).ratio,'Color',colors(i,:));
                h(i).ratio_inset = xline(AX.ratio_inset,t(1)); %#ok<AGROW>
                hold(AX.ratio_inset(1),'off');
                set(AX.ratio_inset,'XLim',trange,'YLim',neuron_interped(ii).ratio_ylim,'XTick',[],'YTick',[],'Box','on');
            end
            % scale bar and cosmetics for tracker plot
            axis(AX.tracker,'equal'); axis(AX.tracker,BOX.tracker);
            set(AX.tracker,'XTick',[],'YTick',[],'Box','on');
            labeledBar(AX.tracker,BAR.tracker.x,BAR.tracker.y,BAR.tracker.str,'below','k',{'lineWidth',4});            
            title(AX.tracker,'neuron trajectory');
            hold(AX.tracker,'off');
            % ======== metadata ========
%             plottitle = fileparts(tsr.tracker.filename);
            maxlen = 100;
            if length(plottitle)>maxlen
                plottitle = fullfile(plottitle(1:2),['...' plottitle(end-(maxlen-5):end)]);
            end
            str = {plottitle;
                ['time: ' sprintf('%0.2f',t(1))]};
            h_title = text(AX.gr(1),'Position',[0.5 16.25]*u,'String',str,'Units','pixels',...
                'HorizontalAlignment','left','VerticalAlignment','middle',...
                'FontSize',12,'FontWeight','bold','Interpreter','none');

            % play movie
            for j=1:length(t)
                if ~ishandle(FIG)
                    return;
                end
                % ======== update behavior ========
                if t(j)>=tsr.behavior(bn).et(1) && t(j)<=tsr.behavior(bn).et(end) 
                    im_behav = tsr.getBehaviorImage(v,t(j),bn);
                    im_behav = flip(flip(im_behav,1),2);
                else
                    im_behav = zeros(v.Width,v.Height);
                end
                h_behav.CData = im_behav;
                if ~isempty(v)
                    str = ['behavior camera (AVI time: '...
                        datestr(v.CurrentTime/(24*60*60),'MM:SS') ')'];
                    AX.behav.Title.String = str;
                end
                % ======== update tracker and activity ========
                for i=1:n
                    ii = neuron_id(i);
                    % update XData for all saved plots
                    h(i).tracker.XData = neuron_interped(ii).loc(1,j);
                    h(i).tracker.YData = neuron_interped(ii).loc(2,j);
                    h(i).gr.Value = t(j);
                    h(i).gr_inset.Value = t(j);
                    h(i).ratio.Value = t(j);
                    h(i).ratio_inset.Value = t(j);
                    % update axes properties
                    xl = [t(j)-twidth t(j)+twidth];
                    xt = [t(j)-twidth t(j) t(j)+twidth];
                    xtl = {['-' sprintf('%0.1f',twidth) 's'],...
                        [sprintf('%0.1f',t(j)) 's'],...
                        ['+' sprintf('%0.1f',twidth) 's']};
%                     if j<=jwidth
%                         yl_gr = [min(min(neuron_interped(ii).greenrate(1:jwidth*2)),min(neuron_interped(ii).redrate(1:jwidth*2)))...
%                             max(max(neuron_interped(ii).greenrate(1:jwidth*2)),max(neuron_interped(ii).redrate(1:jwidth*2)))];
% %                         yl_ratio = [min(neuron_interped(ii).ratio(1:jwidth*2))...
% %                             max(neuron_interped(ii).ratio(1:jwidth*2))];
%                     elseif j>=length(t)-jwidth
%                         yl_gr = [min(min(neuron_interped(ii).greenrate(end-jwidth*2:end)),min(neuron_interped(ii).redrate(end-jwidth*2:end)))...
%                             max(max(neuron_interped(ii).greenrate(end-jwidth*2:end)),max(neuron_interped(ii).redrate(end-jwidth*2:end)))];
% %                         yl_ratio = [min(neuron_interped(ii).ratio(end-jwidth*2:end))...
% %                             max(neuron_interped(ii).ratio(end-jwidth*2:end))];
%                     else
%                         yl_gr = [min(min(neuron_interped(ii).greenrate(j-jwidth:j+jwidth)),min(neuron_interped(ii).redrate(j-jwidth:j+jwidth)))...
%                             max(max(neuron_interped(ii).greenrate(j-jwidth:j+jwidth)),max(neuron_interped(ii).redrate(j-jwidth:j+jwidth)))];
% %                         yl_ratio = [min(neuron_interped(ii).ratio(j-jwidth:j+jwidth))...
% %                             max(neuron_interped(ii).ratio(j-jwidth:j+jwidth))];
%                     end
                    set(AX.gr(i),'XLim',xl,'XTick',[]);
                    set(AX.ratio,'XLim',xl,'XTick',xt,'XTickLabel',xtl);
                end
                % ======== update metadata ========
                str = {plottitle;...
                    ['time: ' sprintf('%0.2f',t(j))]};
                h_title.String = str;
                drawnow;
                if savemovie
                    writeVideo(vw,getframe(FIG)); %#ok<UNRCH>
                end
            end
            
            if savemovie
                close(vw); %#ok<UNRCH>
            end
        end
        
        function playMovie_custom_plot(tsr, varargin)
            % NAME-VALUE inputs:
            %    'neuron_id': default is all tracked neurons
            % 	 'trange': default is tsr.tracker.start/endTime
            % 	 'twidth': time window to display photon rates/ratio (unit:
            %              s); default is 10s before and after 
            %	 'framerate': default is behavior video frame rate (unit: Hz)
            %	 'savemovie',true/[false]: whether to save .mp4 movie to disk
            %	 'fpath': specify full .mp4 path; if savemovie=true and fpath
            %             is not specified, will automatically save to data
            %             folder
            %	 'bn': for data sets with more than one behavior video saved
            
            % TODO add overlaytracker option? (with inset for z)
            
            % default playback and saving settings
            neuron_id = 1:numel(tsr.neuron);
            trange = [tsr.tracker.startTime tsr.tracker.endTime];
            twidth = 3; % time window before and after for neuron activity display (unit: s)
            framerate = [];
            savemovie = false;
            fpath = []; %#ok<NASGU>
            bn = 1; % for when there's multiple behavior videos
            
            assignApplicable(varargin);
            
            % get display params
            if isempty(framerate)
                framerate = tsr.behavior(bn).framerate_Hz;
            end
            jwidth = round(twidth*framerate);
            
            % prepare video reader for behavior images
            v = VideoReader(tsr.behavior(bn).videoFileName);
            
            % prepare video writer if saving movie
            if savemovie
                if isempty(fpath) %#ok<UNRCH>
                    % save to data folder if fpath not specified
                    fpath = strrep(tsr.tracker.filename,'_tracker.txt','_custom_video.mp4');
                end
                vw = VideoWriter(fpath,'MPEG-4');
                vw.FrameRate = framerate;
                open(vw);
            end
            
            vt_range = 845:1000;

            % prepare interped tracks and neuron activities
            tstart = trange(1);
            tstop = trange(2);
            t = tstart:(1/framerate):tstop;
            for i=neuron_id
                % Get interped trajectory and plotting params 
                % Get the trajectory in the same frame rate as the video:
                loc = interp1(tsr.tx,lowpass1D(tsr.neuron(i).loc(1:2,:),median(diff(t))/3)',t)'; % 
                valid = interp1(tsr.tx,double(tsr.tracker.tracking),t,'nearest');
                valid(isnan(valid)) = 0;
                loc(:,~valid) = NaN;
                xrange = range(loc(1,:)); % Used to set the plot range of the trajectory in x 
                yrange = range(loc(2,:)); % Used to set the plot range of the trajectory in y
                traj_box = round([min(loc(1,:))-xrange*0.05 max(loc(1,:))+xrange*0.05 ...
                    min(loc(2,:))-yrange*0.05 max(loc(2,:))+yrange*0.05]);
                % get interped activity and plotting params
                ratio = interp1(tsr.tx,tsr.neuron(i).ratio_div_baseline,t);
                ratio_ylim = [0 max(ratio(vt_range))+.2];

                % gather into struct
                neuron_interped(i).loc = loc; %#ok<AGROW>
                neuron_interped(i).valid = valid; %#ok<AGROW>
                neuron_interped(i).traj_box = traj_box'; %#ok<AGROW> % transpose for easier concatenation later
                neuron_interped(i).ratio = ratio; %#ok<AGROW>
                neuron_interped(i).ratio_ylim = ratio_ylim; %#ok<AGROW>
            end
            traj_box = [neuron_interped(:).traj_box]';
            tmp1 = min(traj_box,[],1);
            tmp2 = max(traj_box,[],1);
            traj_box = [tmp1(1) tmp2(2) tmp1(3) tmp2(4)];
            
            % prepare figure window and axes
            n = numel(neuron_id);

            set(0,'Units','pixels');
            W = get(0,'screensize');
            W = W(3);
            
            % Condition to get the correct figure size
            if n <= 2 
                u = min(40,floor(W/(8*3+12))); % n in here has to be at least 3 to get the correct figure size
                fig = figure('Units','pixels','Position',[1.5 1.5 25 10]*u);
            else
                u = min(40,floor(W/(8*n+12)));
                fig = figure('Units','pixels','Position',[1.5 1.5 8*(n-1)+9 18]*u);
            end

            ax_behav = axes('Units','pixels','Position',[1 1 7 7]*u);  % Plot position of the behavior video 
            ax_traj = axes('Units','pixels','Position',[9 1 7 7]*u);   % Plot position of the trajectory video
            ax_ratio = axes('Units','pixels','Position',[17 1 7 7]*u); %#ok<AGROW>
            ax_ratio_inset = axes('Units','pixels','Position',[21 6 3 2]*u); %#ok<AGROW>

            % prepare first frames and plot params
            axis(ax_behav,'equal','tight','ij');
            % Prepare the trajectory panel 
            hold(ax_traj,'on');
            colors = viridis(n); % define colors for each neuron
            
            custom_trange = [t(vt_range(1)) t(vt_range(end))];
            for i=1:n
                ii = neuron_id(i);
                
                % ======== whole trajectory (left bottom plot) ========
                h(i).traj = plot(ax_traj, neuron_interped(ii).loc(1,vt_range),neuron_interped(ii).loc(2,vt_range),...
                    'DisplayName',['neuron ' num2str(ii)], 'Color', colors(i,:)); %#ok<AGROW>
                h(i).tracker = plot(ax_traj, neuron_interped(ii).loc(1,vt_range(1)),neuron_interped(ii).loc(2,vt_range(1)),'ro','LineWidth',1.5); %#ok<AGROW>
                title(ax_traj, 'neuron trajectory');
%                 set(ax_traj, 'XTick',[], 'YTick',[]); 
                
                xl = [t(vt_range(1))-twidth t(vt_range(1))+twidth];
                yl = ratio_ylim;
%                 yl = [0 max(neuron_interped(ii).ratio(1:jwidth*2))*1.1];

                % ======================== ratio/baseline ========================
                hold(ax_ratio,'on');
                plot(ax_ratio, t(vt_range), neuron_interped(ii).ratio(vt_range), 'Color', colors(i,:), 'LineWidth', 2); 
                h(i).ratio = plot(ax_ratio,[t(vt_range(1)) t(vt_range(1))], neuron_interped(ii).ratio_ylim); %#ok<AGROW>
                hold(ax_ratio,'off');
                
                set(ax_ratio,'XLim',xl,'YLim', yl, 'XTick',[t(vt_range(1))-twidth t(vt_range(1)) t(vt_range(1))+twidth],...
                    'XTickLabel',{['-' sprintf('%0.1f',twidth) 's'],...
                    [sprintf('%0.1f',t(vt_range(1))) 's'],['+' sprintf('%0.1f',twidth) 's']});
                title(ax_ratio, 'Ratio: R/R_0');
            
                % ======== ratio/baseline inset ========
                hold(ax_ratio_inset(1),'on');
                plot(ax_ratio_inset,t(vt_range),neuron_interped(ii).ratio(vt_range), 'Color', colors(i,:));
                h(i).ratio_inset = plot(ax_ratio_inset,[t(vt_range(1)) t(vt_range(1))],neuron_interped(ii).ratio_ylim); %#ok<AGROW>
                hold(ax_ratio_inset(1),'off');
                set(ax_ratio_inset,'XLim',custom_trange,'YLim',neuron_interped(ii).ratio_ylim,'XTick',[],'YTick',[]);
            end
            p1 = [80,1250]; p2 = [80,1350]; % scale bar: 65px = 200um
            plot(ax_traj, [p1(2),p2(2)],[p1(1),p2(1)],'Color','k','LineWidth',3);
            text(ax_traj, 1250, 65, '100 \mum','Color', 'k','FontSize',10);

            hold(ax_traj,'off');
            
            % play video
            for j = vt_range % 1:length(t)
                if ~ishandle(fig)
                    return;
                end
                % update behavior
                if t(j)>=tsr.behavior(bn).et(1) && t(j)<=tsr.behavior(bn).et(end) 
                    behavIm = tsr.getBehaviorImage(v,t(j),bn);
                    behavIm = flip(flip(behavIm,1),2);
                else
                    behavIm = zeros(v.Width,v.Height);
                end
                % ======== behavior video ========
                imagesc(ax_behav, behavIm); colormap(ax_behav,gray(256)); caxis(ax_behav,[0 Inf]);
                hold(ax_behav, 'on');
                p1 = [440,550]; p2 = [440,615]; % scale bar: 65px = 200um
                plot(ax_behav, [p1(2),p2(2)],[p1(1),p2(1)],'Color','w','LineWidth',3);
                text(ax_behav, 520, 460, '200 \mum','Color',[1 1 0.99],'FontSize',10);
%                 text(550, 460, '200 \mum','Color',[1 1 0.99],'FontSize',8);
                hold(ax_behav, 'off');
                
                set(ax_behav,'XTick',[],'YTick',[]);
                title(ax_behav,['behavior camera (AVI time: ' datestr(v.CurrentTime/(24*60*60),'MM:SS') ')']); % assuming video <1hr
                % update track and activity for each neuron
                for i=1:n
                    ii = neuron_id(i);
                    % update XData for all saved plots
                    h(i).tracker.XData = neuron_interped(ii).loc(1,j);
                    h(i).tracker.YData = neuron_interped(ii).loc(2,j);
                    h(i).ratio.XData = [t(j) t(j)];
                    h(i).ratio_inset.XData = [t(j) t(j)];
                    % update axes properties
                    xl = [t(j)-twidth t(j)+twidth];
                    xt = [t(j)-twidth t(j) t(j)+twidth];
                    xtl = {['-' sprintf('%0.1f',twidth) 's'],...
                        [sprintf('%0.1f',t(j)) 's'],...
                        ['+' sprintf('%0.1f',twidth) 's']};

                    set(ax_ratio,'XLim',xl,'XTick',xt,'XTickLabel',xtl);
                end
                
                % metadata
                fname = tsr.tracker.filename;
                [filepath,name,ext] = fileparts(fname);
                filepath = strrep(filepath, '\','/');
                strIdx = strfind(filepath, '/');
                title_name = filepath(strIdx(4)+1:end);

                str = {title_name;... % fileparts(fname) is replaced by title_name
                    ['time: ' sprintf('%0.2f',t(j))]};
                text(ax_behav,'Position',[1.5 8.25]*u,'String',str,'Units','pixels',...
                    'HorizontalAlignment','left','VerticalAlignment','middle',...
                    'FontSize',12,'FontWeight','bold','Interpreter','none');
                drawnow;
                if savemovie
                    writeVideo(vw,getframe(fig)); %#ok<UNRCH>
                end
            end
            
            if savemovie
                close(vw); %#ok<UNRCH>
            end
        end
        
        function playMovieOld(tsr,varargin)
            %PLAYMOVIE plays behavioral movie alongside green/red photon rates and ratio/baseline
            % additional options (Name,Value pair arguments):
            %   tstart,tstop  - play fractions of whole experiment (must be within tsr.tx)
            %   twidth        - show neuron activity __seconds before and after current frame
            %   playbackspeed - default is 10x real time; doesn't affect frame rate of saved movie
            %   savemovie     - if set to true, will prompt user to select where to save .mp4 movie
            %   fpath         - can also specify full .mp4 file path
            
            % TODO add overlaytracker option (with inset for z)
            
            % default playback and saving settings
            tstart = tsr.video.tx(1);
            tstop = tsr.video.tx(end);
            twidth = 10; % time window before and after (unit: s)
            playbackspeed = 10; % >1: faster <1: slower
            savemovie = false;
            fpath = []; %#ok<NASGU>
            assignApplicable(varargin);
            % get file path if saving movie
            if savemovie
                if isempty(fpath) %#ok<UNRCH>
                    % prompt user to select where to save movie
                    fdir = fileparts(tsr.tracker.filename);
                    [fname,fdir] = uiputfile(fullfile(fdir,'*.mp4'),'select behavioral video');
                    fpath = fullfile(fdir,fname);
                end
                v = VideoWriter(fpath,'MPEG-4');
                v.FrameRate = round(1/mean(diff(tsr.video.tx)));
                open(v);
            end
            % get helper variables
            jstart = find(tsr.video.tx>=tstart,1,'first');
            jstop = find(tsr.video.tx<=tstop,1,'last');
            jwidth = round(twidth/mean(diff(tsr.video.tx)));
            % set y axis limits
            ylim1 = [0 max(max(tsr.video.greenrate(jstart:jstop)),max(tsr.video.redrate(jstart:jstop)))];
            ylim2 = [0 max(tsr.video.ratio_div_baseline(jstart:jstop))];
            % prepare figure
            fig = figure('Units','normalized','Position',[0.1 0.1 0.4 0.3]);
            ax1 = axes('Position',[0.05 0.1 0.35 0.8]);
            ax2 = axes('Position',[0.45 0.55 0.5 0.35]);
            ax3 = axes('Position',[0.45 0.1 0.5 0.35]);
            % prepare for movie
            j = jstart;
            axes(ax1); % behavioral image
            imagesc(squeeze(tsr.videocube(:,:,j)));
            xticks([]); yticks([]); axis equal tight;
            axes(ax2); % green and red
            plot(1:length(tsr.video.tx),tsr.video.greenrate,'g'); hold on;
            plot(1:length(tsr.video.tx),tsr.video.redrate,'r');
            h1 = plot([j j],ylim1,'k'); hold off;
            xlim([j-jwidth j+jwidth]); ylim(ylim1); xticks([]);
            title('photon counts');
            axes(ax3); % ratio
            plot(tsr.video.ratio_div_baseline,'k'); hold on;
            h2 = plot([j j],ylim2,'k'); hold off;
            xlim([j-jwidth j+jwidth]); ylim(ylim2); xticks([j-jwidth j j+jwidth]);
            xticklabels({['-' sprintf('%0.1f',twidth) 's'],...
                sprintf('%0.1f',tsr.video.tx(j)),...
                ['+' sprintf('%0.1f',twidth) 's']});
            title('ratio/baseline');
            % play movie
            for j=jstart:jstop
                if ~ishandle(fig)
                    return;
                end
                % update behavioral image
                imagesc(ax1,squeeze(tsr.videocube(:,:,j)));
                xticks(ax1,[]); yticks(ax1,[]); axis(ax1,'equal','tight');
                % update green/redrates
                h1.XData = [j j];
                xlim(ax2,[j-jwidth j+jwidth]); ylim(ax2,ylim1); xticks(ax2,[]);
                % update ratio/baseline
                h2.XData = [j j];
                xlim(ax3,[j-jwidth j+jwidth]); ylim(ax3,ylim2); xticks(ax3,[j-jwidth j j+jwidth]);
                xticklabels(ax3,{['-' sprintf('%0.1f',twidth) 's'],...
                    sprintf('%0.1f',tsr.video.tx(j)),...
                    ['+' sprintf('%0.1f',twidth) 's']});
                drawnow;
                if savemovie
                    f = getframe(fig); %#ok<UNRCH>
                    writeVideo(v,f);
                else
                    pause(0.1/playbackspeed);
                end
            end
            if savemovie
                close(v); %#ok<UNRCH>
            end
        end % deprecated



        
        function [y, yset, yset_nan_infill,num_elements] = alignDataToMotionPhase(tsr, t, y, phaseAxis, trange, varargin)
            % [y, yset] = alignDataToMotionPhase(tsr, t, y, phaseAxis, trange)
            % aligns data according to peristaltic phase (must run
            % tsr.findCounterMovements first)
            % t: time associated with y data (could be tsr.tx but doesn't
            % have to be)
            % y: Nxt (or 1D) data to be aligned. matches t in dimension
            % phaseAxis: axis to report aligned data, phase = 0 is at time
            %            brain is farthest back +/-2pi are next/previous
            %            brain furthest backs (set this like: linspace(-2*pi,2*pi,100);)
            % note - phaseAxis must be no more than 1 cycle (2pi) or will
            % overwrite!
            % trange (optional): time range within which to perform this
            %           operation
            % varargin: ignore_valid -> false (if true, use all
            % counterinds)
            % y - average of all readings
            % yset - N x len(t) - each aligned trace

            existsAndDefault('trange',[-Inf Inf]);
            behav_label = [];
            ignore_valid = false;
            assignApplicable(varargin);

            valid = tsr.com.perivalid | ignore_valid;
            if (existsAndDefault('behav_label', []))
                if (ischar(behav_label))
                    if (lower(behav_label(1)) == 'f')
                        behav_label = 1;
                    else
                        behav_label = 2;
                    end
                end
                valid = valid & (tsr.com.behav_label == behav_label);
            end

            ci = tsr.com.counterInd(valid(tsr.com.counterInd));
            ci = ci(any(tsr.tx(ci) >= trange(:,1) & tsr.tx(ci) <= trange(:,2),1));

            if (isempty(ci))
                y = nan([size(y,1:ndims(y)-1) length(phaseAxis)]);
                yset = zeros([size(y,1:ndims(y)-1) 0 length(phaseAxis)]);
                yset_nan_infill = yset;
                num_elements = yset;
                return;
            end

            p = unwrap(tsr.com.periphase); %p is phase vs. tx
            p_y = interp1(tsr.tx, p, t,'linear','extrap'); %p_y is phase vs. t (assoc. w/ y)
            phaseEdges = conv(phaseAxis, [0.5, 0.5], 'valid');
            phaseEdges = [1.5*phaseAxis(1)-0.5*phaseAxis(2) phaseEdges 1.5*phaseAxis(end)-0.5*phaseAxis(end-1)];
            %yys = zeros(length(ci),length(phaseAxis)); %each y per sweep

            %pbins is npts x ndim -- ind1,ind2;ind1,ind2;...
            pbins = NaN([length(p_y), 2]);
            for k = 1:length(ci)
                kbin = discretize(p_y, phaseEdges+p(ci(k)));
                valid = isfinite(kbin);
                pbins(valid,2) = kbin(valid);
                pbins(valid,1) = k;
            end

            in_a_bin = all(isfinite(pbins),2);

            pbins = pbins(in_a_bin, :);

            numphase = length(phaseAxis);


            function m = fastmean(x)
                if (numel(x) <= 1)
                    m = x;
                else
                    m = sum(x)/numel(x);
                end
            end

            if (size(y,2) == 1)
                y = y';
            end
            
            tind = ndims(y);
            S.type = '()';
            S.subs = cell(1,tind);
            [S.subs{1:tind}] = deal(':');
            S.subs{tind} = find(in_a_bin);
            y = y(S.subs{:});
            
            if (true || nnz(size(y) > 1) > 1)
                
                %adapted from
                %https://www.mathworks.com/matlabcentral/answers/332752-n-d-matrix-accumarray
                
                
                tind = ndims(y);
                sz = size(y,1:(tind-1));
                indices = arrayfun(@(s) 1:s, size(y), 'UniformOutput', false); %vector of indices in each dimension
                [indices{:}] = ndgrid(indices{:});
                indices{end} = permute(repmat(pbins(:,1), [1 sz]), circshift(1:tind,-1));
                indices{end+1} = permute(repmat(pbins(:,2), [1 sz]), circshift(1:tind,-1));
              

                inds = cell2mat(cellfun(@(v) v(:), indices, 'UniformOutput', false));
                %apply accumarray with default sum function
               % yset = accumarray(inds, y(:),[sz, length(ci), length(phaseAxis)],@fastmean,NaN);
                yset = accumarray(inds, y(:),[sz, length(ci), length(phaseAxis)],@sum,0);
                num_elements = accumarray(inds, ones(size(y(:))),[sz, length(ci), length(phaseAxis)],@sum,0);
                yset = yset./num_elements;

                if nnz(size(y) > 1) <= 1
                    y = squeeze(mean(yset, tind,"omitnan"))'; %remember tind is the counterind index for yset 
                    yset = squeeze(yset);
                    num_elements = squeeze(num_elements);
                else
                    y = squeeze(mean(yset, tind,"omitnan")); %remember tind is the counterind index for yset 
                end
                if (nargout > 2)
                    yset_nan_infill = nan_infill_nd(yset);
                    valid = isfinite(yset_nan_infill);
                    yset_nan_infill(~valid) = 0;
                end

                return;
            end

            
          
        end

        function [y, yset] = alignDataToMotionTime(tsr, t, y, timeAxis, trange)
            % [y, yset] = alignDataToMotionTime(tsr, t, y, timeAxis, trange)
            % aligns data according to time relative to peristaltic cycle (must run
            % tsr.findCounterMovements first)
            % t: time associated with y data (could be tsr.tx but doesn't
            % have to be)
            % y: 1D data to be aligned. matches t in dimension
            % timeAxis: axis to report aligned data, time = 0 is time
            %            brain is farthest back 
            % trange (optional): time range within which to perform this
            %           operation
            % y - average of all readings
            % yset - N x len(t) - each aligned trace


            existsAndDefault('trange',[-Inf Inf]);
            existsAndDefault('trange',[-Inf Inf]);
            ci = tsr.com.counterInd(tsr.com.perivalid(tsr.com.counterInd));
            ci = ci(any(tsr.tx(ci) >= trange(:,1) & tsr.tx(ci) <= trange(:,2),1));
            
            if (nnz(size(y) > 1) > 1)
                yout = zeros(size(y,1), length(timeAxis));
                yset = zeros(size(y,1), length(timeAxis), length(ci));
                for j = 1:size(y,1)
                    [yy,yys] = alignDataToMotionTime(tsr, t, y(j,:), timeAxis, trange);
                    yout(j,:) = yy;
                    yset(j,:,:) =yys';
                end
                y = yout;
                return;
            end

            
           % ci = tsr.com.counterInd(tsr.com.perivalid(tsr.com.counterInd));
            %ci = ci(tsr.tx(ci) >= trange(1) & tsr.tx(ci) <= trange(2));
            yset = zeros(length(ci),length(timeAxis));
            for j = 1:length(ci)
                yset(j,:) = interp1(t,y,tsr.tx(ci(j)) + timeAxis);
            end
            y = mean(yset,'omitnan');
        end


        function tsr = markForwardStrides(tsr,time_interval,period)
            smoothf = 0.99;
            existsAndDefault('time_interval', [max(min(tsr.behavior.sleap.tx), min(tsr.tx)), min(max(tsr.behavior.sleap.tx), max(tsr.tx))]);
            insleaprange = (tsr.tx > min(tsr.behavior.sleap.tx) & tsr.tx < max(tsr.behavior.sleap.tx));
            if (islogical(time_interval))
                inti = find(time_interval & insleaprange);
            else
                inti = find( insleaprange & any(tsr.tx >= time_interval(:,1) & tsr.tx <=  time_interval(:,2),1) );
            end

            vfwd = tsr.com.vfwd-tsr.com.vfwd_smooth;
            tx = tsr.tx(inti);


            [vmax,locs] = findpeaks(vfwd(inti), "MinPeakHeight",10); 
            vmax = vmax(dl(inti(locs))>0); %only consider points moving forward when spot-gut axis is extended
            cutoff = max(10,mean(vmax) - 2*std(vmax));
            [~,locs] = findpeaks(vfwd(inti), "MinPeakHeight",cutoff);
            locs = unique(locs(dl(inti(locs))>0));

            valid = true(size(locs));
            dt = diff(tx(locs));
            for j = find(dt(:)' < period/2)
                if (vfwd(inti(locs(j))) > vfwd(inti(locs(j+1))))
                    valid(j+1) = false;
                else
                    valid(j) = false;
                end
            end
            locs = locs(valid);

            
%            plot(tx,vfwd(inti), tx(locs), vfwd(inti(locs)),'r.')
            fwdi = 0*locs;
            for j = 1:length(locs)
                fwdi(j) = inti(find(vfwd(inti(1:locs(j))) <= 0, 1, 'last')); %find the last stationary point before current peak
              %  plot(tx,vfwd(inti), tx(locs), vfwd(inti(locs)),'r.', tsr.tx(fwdi(j)), tsr.com.vfwd(fwdi(j)),'go'); pause(0.1)
            end
            %if (~isfield(tsr.com, 'fwdi')|| isempty(tsr.com.fwdi))
                tsr.com.fwdInd = unique(fwdi);
            %else
             %   tsr.com.fwdi = sort([tsr.com.fwdi fwdi]);
            %end
             



             

        end

        function tsr = markBackwardStrides(tsr,time_interval,period)
            smoothf = 0.99;
            existsAndDefault('time_interval', [max(min(tsr.behavior.sleap.tx), min(tsr.tx)), min(max(tsr.behavior.sleap.tx), max(tsr.tx))]);
            insleaprange = (tsr.tx > min(tsr.behavior.sleap.tx) & tsr.tx < max(tsr.behavior.sleap.tx));

            if (islogical(time_interval))
                inti = find(insleaprange & time_interval);
            else
                inti = find(insleaprange & any(tsr.tx >= time_interval(:,1) & tsr.tx <  time_interval(:,2),1) );
            end

            gut = tsr.behavior.sleap.gut;%+ interp1(tsr.tx, tsr.tracker.stageloc(1:2,:)', tsr.behavior.sleap.tx)';
            spot = tsr.behavior.sleap.spot;% + interp1(tsr.tx, tsr.tracker.stageloc(1:2,:)', tsr.behavior.sleap.tx)';
            gut = csaps(tsr.behavior.sleap.tx, gut, smoothf, tsr.tx);
            spot = csaps(tsr.behavior.sleap.tx, spot, smoothf, tsr.tx);

            lgs = sqrt(sum((spot-gut).^2)); %gut spot length
            fwd_dir = (spot-gut)./lgs;

            ltt = gradient(gradient(lgs));
            zcr = [false diff(ltt < 0) ~= 0] & (abs(gradient(lgs,tsr.dt)) > std(gradient(lgs(inti),tsr.dt)));
            midlength = csaps(tsr.tx(zcr), lgs(zcr), .15,tsr.tx);
            dl = lgs-midlength;
            %plot(tsr.tx, dl)


            short_max = findpeaks(-dl(inti),"MinPeakHeight",5); %shortest length valleys in time range
            cutoff = mean(short_max) - 2*std(short_max);
            [~,locs] = findpeaks(-dl(inti),"MinPeakHeight",cutoff);
            tx = tsr.tx(inti);

            existsAndDefault('period',median(diff(tx(locs))));

            %backward crawling - find the local minima of gut body length,
            %then find the fastest forward movement (opposite the direction of backward crawling) between each point,
            %then find the place where the forward velocity is zero before 

            imax = zeros(size(locs));
            locs(end+1) = length(inti);
            for j = 1:length(imax)
                rng = inti(locs(j):locs(j+1));
                rng = rng(tsr.tx(rng)-tsr.tx(rng(1)) < 2*period);
                [~,I] = max(tsr.com.vfwd(rng));
                imax(j) = I + inti(locs(j))-1;
            end

            valid = true(size(imax));
            dt = diff(tsr.tx(imax));
            for j = find(dt(:)' < period/2)
                if (tsr.com.vfwd(imax(j)) > tsr.com.vfwd(imax(j+1)))
                    valid(j+1) = false;
                else
                    valid(j) = false;
                end
            end
            imax = imax(valid);

        
            
            bcki = 0*imax;
            for j = 1:length(imax)
                bcki(j) = find(tsr.com.vfwd(1:imax(j)) <= 0, 1, 'last'); %find the last stationary point before current peak
            end
            tsr.com.bckInd = unique(bcki);
             



             

        end

        function tsr = findCounterMovements (tsr, varargin)
            % tsr = findCounterMovements (tsr, varargin)
            % finds times where neuron motion is counter to prevailing
            % motion; returns zero crossing from low to high (farthest
            % opposite position)
            % varargin = lp_sigma (smoothing for travel direction - default
            % 2)
            % period = expected period of peristalsis (sets lp_sigma)
            % mindist = minimum distance between separate strides (default
            % 20 microns)
            %
            % perivalid flags periods within 50-150% of mean period
            % can go badly when stride duration varies over course of
            % experiment
            period = [];
            varargin = assignApplicable(varargin);
            if (~isempty(period))
                lp_sigma = period;
            else
                lp_sigma = 3; %seconds
            end
            mindist = 20; %microns - minimum travel per stride

                        
            varargin = assignApplicable(varargin);

 
              vl = deriv(tsr.com.loc(1:2,1:10:end), 0.1*lp_sigma/tsr.dt);
              d = cumsum(sqrt(sum(vl.^2)))*tsr.dt*10;
              xs = csaps(d, tsr.com.loc(1:2,1:10:end), .25,d); %spline smoothed location, vs distance traveled
              vsmooth = deriv(xs,0.1*lp_sigma/tsr.dt)/(10*tsr.dt);
              vsmooth =  interp1(tsr.tx(1:10:end), vsmooth', tsr.tx, 'linear','extrap')';
              vhat = vsmooth./sqrt(sum(vsmooth.^2));
 
            if (isfield(tsr.behavior, 'sleap'))
                gut = interp1(tsr.behavior.sleap.tx, tsr.behavior.sleap.gut', tsr.tx)';
                spot = interp1(tsr.behavior.sleap.tx, tsr.behavior.sleap.spot', tsr.tx)';
                fwd_dir = (spot-gut)./sqrt(sum((spot-gut).^2));
                tsr.com.fwd_dp = sum(vhat.*fwd_dir); %long time scale movement alignment to gut-spot axis
                
                
                xsh = csaps(d, tsr.com.loc(1:2,1:10:end), .99,d); %spline smooth tracker location
                v = deriv(xsh, 0.1*lp_sigma/12/tsr.dt)/(10*tsr.dt);
                v = interp1(tsr.tx(1:10:end), v', tsr.tx, 'linear','extrap')';

                tsr.com.vfwd = sum(v.*fwd_dir);
                tsr.com.vfwd_smooth = sum(vsmooth.*fwd_dir);
                tsr.com.fwd_dir = fwd_dir;


                %use behavior time spacing
                gut = tsr.behavior.sleap.gut + interp1(tsr.tx, tsr.tracker.stageloc(1:2,:)', tsr.behavior.sleap.tx,'nearest', 'extrap')';
                dtsleap = median(diff(tsr.behavior.sleap.tx));
                tsr.com.vgut_smooth = sum(interp1(tsr.behavior.sleap.tx, deriv(gut, lp_sigma/dtsleap)', tsr.tx)'/dtsleap.*fwd_dir);
                tsr.com.vgut = sum(interp1(tsr.behavior.sleap.tx, deriv(gut, 0.2*lp_sigma/dtsleap)', tsr.tx)'/dtsleap.*fwd_dir);
                
                % mark and phase align strides
                tsr = tsr.markStridesWithBehaviorLabels(varargin{:});
                tsr = tsr.alignStrides('behav_type','forward',varargin{:});
                tsr = tsr.alignStrides('behav_type','backward',varargin{:});
                tsr.com.counterInd = [];
                for j = 1:length(tsr.com.strides)
                    if (tsr.com.strides(j).behav_label == 1 || tsr.com.strides(j).behav_label == 2)
                        tsr.com.counterInd = [tsr.com.counterInd tsr.com.strides(j).txinds(1)];
                    end
                end

                tsr.com.behav_label = NaN(size(tsr.tx));
                for j = 1:length(tsr.com.strides)
                    tsr.com.behav_label(tsr.com.strides(j).txinds) = tsr.com.strides(j).behav_label;
                end


                tsr.com.perivalid = false(size(tsr.tx));
                tsr.com.perivalid([tsr.com.strides.txinds]) = true;

                %make phase continuously increasing
                t = tsr.tx([tsr.com.strides.txinds]);
                p = unwrap([0 tsr.com.strides.p_warped]);
                p = p(2:end);
                ii = find(p >=0 ,1,'first');
                p(1:(ii-1)) = 0;
                valid = true(size(t));
                while(any(valid & (diff([0 p])<0)))
                    I = find(valid & (diff([0 p])<0), 1, 'first');
                    
                    lowval = p(I);
                    highval = p(I-1);
                    start = find(p(1:I) < lowval, 1, 'last');
                    stop = I + find(p(I:end) > highval, 1, 'first');
                    if (isempty(stop))
                        stop = length(p);
                    end
                    if (isempty(start))
                        start = 1;
                    end
                    valid(start:stop) = false;
                end
%                 sp = spaps(find(valid),p(valid),.1);
%                 p = fnval(sp, 1:length(p));
                p = csaps(find(valid), p(valid), 0.9, 1:length(p));
                ii = find(p >=0 ,1,'first');
                p(1:(ii-1)) = 0;
                %now repeat, but with linear interpolation to fix spline
                %overshoot errors
                valid = true(size(t));
                while(any(valid & (diff([0 p])<0)))
                    I = find(valid & (diff([0 p])<0), 1, 'first');
                    lowval = p(I);
                    highval = p(I-1);
                    start = find(p(1:I) < lowval, 1, 'last');
                    stop = I + find(p(I:end) > highval, 1, 'first');
                    if (isempty(stop))
                        stop = length(p);
                    end
                    if (isempty(start))
                        start = 1;
                    end
                    valid(start:stop) = false;
                end


                tsr.com.periphase = interp1(t(valid), p(valid), tsr.tx, 'linear','extrap');




                


                return;
            else
                fwd_dir = vhat;
            end



% 
%             d =  interp1(tsr.tx(1:10:end), d', tsr.tx, 'linear','extrap')';
             vl = deriv(tsr.com.loc(1:2,1:10:end), 0.1*lp_sigma/tsr.dt);
             vl = interp1(tsr.tx(1:10:end), vl', tsr.tx, 'linear','extrap')';
             d = cumsum(sqrt(sum(vl.^2)))*tsr.dt;
             xsh = csaps(d, tsr.com.loc(1:2,:), .99,d); %much less smoothing
             v = deriv(xsh, .01/tsr.dt);


%            vh = vl./repmat(sqrt(sum(vl.^2,1)),[3 1]);

%             tsr.com.vhat_smooth = vhat;
%             tsr.com.fwd_dir = fwd_dir;
 %           v = lowpass1D(tsr.com.vel(1:2,:),.1/tsr.dt);
            tsr.com.local_speed = sum(v.*fwd_dir,1);

            s = tsr.com.local_speed;% .* sign(tsr.com.fwd_index); %invert when larva is backing up
        

            [lowval,lowloc] = findpeaks(-s,"MinPeakHeight",0); %all negative peaks below 0
            [highval,highloc] = findpeaks(s,"MinPeakHeight",0); %all peaks above 0
            lowvalid = true(size(lowloc));
            highvalid = true(size(highloc));
            
            %consolidate low peaks to pick the lowest between each high
            %peak
            for j = 1:length(lowloc)
                if (~lowvalid(j))
                    continue; %already marked bad
                end
                prevhigh = max(highloc(highloc < lowloc(j)));
                if (isempty(prevhigh))
                    prevhigh = 0;
                end
                nexthigh = min(highloc(highloc > lowloc(j)));
                if (isempty(nexthigh))
                    nexthigh = max(lowloc);
                end
                invalley = lowloc >= prevhigh & lowloc <= nexthigh;
                ind0 = find(invalley,1,'first');
                [~,I] = max(lowval(invalley)); %deepest dip - max is used because of -s
                lowvalid(invalley) = false;
                lowvalid(ind0+I-1) = true;                
            end

            lowval = lowval(lowvalid);
            lowloc = lowloc(lowvalid);

             
            %consolidate high peaks to pick the highest between each low
            %peak
            for j = 1:length(highloc)
                if (~highvalid(j))
                    continue; %already marked bad
                end
                prevlow = max(lowloc(lowloc < highloc(j)));
                if (isempty(prevlow))
                    prevlow = 0;
                end
                nextlow = min(lowloc(lowloc > highloc(j)));
                if (isempty(nextlow))
                    nextlow = max(highloc);
                end
                onmountain = highloc >= prevlow & highloc <= nextlow;
                ind0 = find(onmountain,1,'first');
                [~,I] = max(highval(onmountain)); %deepest dip - max is used because of -s
                highvalid(onmountain) = false;
                highvalid(ind0+I-1) = true;                
            end    

            highval = highval(highvalid);
            highloc = highloc(highvalid);

            phv = interp1(highloc, highval, lowloc,'previous','extrap');
            nhv = interp1(highloc,highval,lowloc,'next','extrap');
            phv(~isfinite(phv)) = nhv(~isfinite(phv));
            nhv(~isfinite(nhv)) = phv(~isfinite(nhv));
            delta = 0.5*(phv+nhv) + lowval; %change in local velocity between previous and next peak and intermediate valley            

            dist = sqrt(sum(diff(tsr.com.loc(:,lowloc),[],2).^2)); %distance between low points
            cdist = [0 cumsum(dist)];

           
            lowvalid = true(size(lowloc));
            for j = 1:(length(lowloc))
                if (~lowvalid(j))
                    continue;
                end
                ingroup = cdist-cdist(j) < mindist & cdist >= cdist(j);
                [~,I] = max(delta(ingroup));
                lowvalid(ingroup) = false;
                lowvalid(j-1+I) = true;
            end
                
            loc = lowloc(lowvalid);
            pkh = delta(lowvalid);


            

            zcr = s>0 & [false s(1:end-1)<0]; %s > 0 and previous element was less than 0

            zcrnext = loc;
            zcrvalid = true(size(zcrnext));
            for j = 1:length(loc)
                ii = find(zcr(loc(j):end), 1, 'first');
                if (~isempty(ii))
                    zcrnext(j) = loc(j)-1+ii;
                else
                    zcrvalid(j) = false;
                end
            end
            zcrnext = zcrnext(zcrvalid & [zcrnext(1:end-1)<loc(2:end) true]); %next zero crossing

            delta = median(diff(zcrnext)); %period in units of index

            rng = round(delta/2); %range to extract dip
            d = cumsum(s);

            dip = 0;
            valid = (zcrnext + rng) < length(d) & zcrnext > (rng + 1);
            for j = find(valid)
                dip = dip + (d((zcrnext(j)-rng):(zcrnext(j)+rng)) - d(zcrnext(j)));
            end
            dip = dip/nnz(valid);


            %i'm not sure what's going on here -- may need to check it out
            nnz(valid)

            [dleft,indmax] = max(dip(1:rng)); %maximum avg value before retrograde movement and time this occurs
            delay = rng-indmax; %time between maximum value and zero crossing 
            dright = dip(delay+rng); %avg value same delay after retrograde movement

            valid = false(size(zcrnext));
            for j = 1:length(zcrnext)
                ind1 = max(1,zcrnext(j)-2*delay);
                ind2 = min(length(d), zcrnext(j) + 2*delay);
                valid(j) = any(d(ind1:zcrnext(j)) > d(zcrnext(j)) + dleft*.5) && any(d(zcrnext(j):ind2) > d(zcrnext(j)) + dright*.5); 
            end
            nnz(valid)

            zcrnext = zcrnext(valid);

            tsr.com.counterInd = zcrnext;
            tsr.com.counterDelta = pkh(valid);
            tcross = tsr.tx(tsr.com.counterInd);

            tau = gradient(tcross);
            tau0 = median(tau);
            
            tauv = (tau > 0.5)*tau0 & tau < 1.5*tau0;

            deltat = tsr.tx - interp1(tcross, tcross, tsr.tx, 'nearest', 'extrap');
            p = interp1(tcross, tau, tsr.tx, 'linear','extrap');
            tsr.com.periphase = 2*pi*deltat./p;
            tsr.com.perivalid = logical(interp1(tcross, double(tauv), tsr.tx, 'nearest', 'extrap'));
            if (numel(tsr.neuron) == 1) %preserve backward compatibility with existing code
                tsr.neuron.counterInd = tsr.com.counterInd;
                tsr.neuron.perivalid = tsr.com.perivalid;
                tsr.neuron.local_speed = tsr.com.local_speed;
                tsr.neuron.counterInd = tsr.com.counterInd;
                tsr.neuron.periphase = tsr.com.periphase;
            end
            return;

        end
        
        tsr = markStridesWithBehaviorLabels(tsr,varargin) % in separate file
        
        function tsr = markStridesWithBehaviorLabels_old(tsr,varargin)
            %DEPRECATED tsr = tsr.markStridesWithBehaviorLabels('valid_trange',trange);
            % trange: dataset-specific time range where larva is uncompressed
            % use Hilbert transform on (vfwd-vfwd_smooth) to detect strides
            
            nphasebins = 36;
            valid_trange = tsr.tx([1 end]);
            vmin = 200; % threshold for behavior labeling
            verbose = 0;
            assignApplicable(varargin);
            
            % get vfwd
            X = tsr.com.vfwd-tsr.com.vfwd_smooth; % subtracting lowpass com movement gives a more sinusoidal signal with a better hilbert transform
            
            valid = isfinite(X);

            X(~isfinite(X)) = 0;

            % do Hilbert transform on -vfwd: phase rollover occurs at vfwd peaks
            sig = hilbert(-X');
%             omega = gradient(unwrap(angle(sig)));
%             inse = abs(sig).^2;
%             insp = angle(sig);
%             fs = 1/tsr.dt; % sampling frequency
%             insf = fs/(2*pi)*omega;
            phi = angle(sig);
            
            % detect cycle by phase cliffs (diff~-2pi)
            dphi = diff(phi);
            dphi = cat(1,NaN,dphi);
            cliffinds = find(valid' & dphi<=-2*pi*0.9);
            stridestartinds = cliffinds;
            
            % segment behavior into forward/backward by thresholding fwd_dp
            behav_labels = zeros(size(tsr.com.fwd_dp));
            behav_labels(tsr.com.fwd_dp>cosd(45) & abs(tsr.com.vfwd_smooth)>=vmin) = 1;
            behav_labels(tsr.com.fwd_dp<-cosd(45) & abs(tsr.com.vfwd_smooth)>=vmin) = 2;
            
            % organize into struct
            nstrides = length(stridestartinds)-1;
            S = repmat(struct('txinds',[],'t',[],'vfwd',[],'tt',[],'vfwd_tt',[],'behav_label',0),nstrides,1);
            for j=1:nstrides
                inds = stridestartinds(j):(stridestartinds(j+1)-1);
                S(j).txinds = inds;
                S(j).vfwd = X(inds);
                t = tsr.tx(inds);
                S(j).t = t; % original tracker time
                S(j).tt = linspace(t(1),t(end),nphasebins); % tracker time resampled to get 1-to-1 map to phase bins
                S(j).vfwd_tt = interp1(S(j).t,S(j).vfwd,S(j).tt);
                k = mode(behav_labels(inds));
                inrange = S(j).t(1)>=valid_trange(1) & S(j).t(end)<=valid_trange(2);
                if all(behav_labels(inds)==k) && inrange
                    % only label as forward/backward if
                    % 1) entire stride is within valid trange (uncompressed), and
                    % 2) entire stride is uncomplicated foward/backward
                    S(j).behav_label = k; 
                end
            end
            
            % save into tsr.com
            tsr.com.tcross = tsr.tx(stridestartinds); % includes both start and end of last stride, so length=nstrides+1
            tsr.com.strides = S;
            tsr.com.phasebins = linspace(0,2*pi,nphasebins);
            
            % plot if verbose=true
            if verbose>0
                figure;
                h(1) = subplot(3,1,1); % vfwd
                plot(tsr.tx,X); hold on;
                xline(tsr.tx(stridestartinds),'r--'); hold off;
                title('data');
                h(2) = subplot(3,1,2); % instantaneous phase
                plot(tsr.tx,phi); hold on;
                xline(tsr.tx(stridestartinds),'r--'); hold off;
                title('\phi'); ylim([-1 1]*1.5*pi); yticks([-1 0 1]*pi); yticklabels({'-\pi','0','\pi'});
                h(3) = subplot(3,1,3); % diff(instantaneous phase)
                plot(tsr.tx,dphi); hold on;
                plot(tsr.tx(cliffinds),dphi(cliffinds),'rx'); hold off;
                title('\Delta\phi'); ylim([-1 1]*2*pi); yticks([-1 0 1]*2*pi); yticklabels({'-2\pi','0','2\pi'});
                linkaxes(h,'x'); xlim(tsr.tx([1 end]));
            end
        end
        
        [S_merged,mergepairs] = detectAndMergeShortStrides(tsr,S,varargin) % in separate file
        
        function tsr = alignStrides(tsr,varargin)
            %tsr = tsr.alignStrides('behav_type','forward')
            % align strides using our own spline time warping code
            % default template: mean of all strides==behav_type
            % can specify template by setting input 'template_trange'
            % verbose = 0: (default) show nothing
            %           1: show final comparison plot
            %           2: show each stride's alignment details with 1s pause
            %           3: save 2 as a movie to behavior_video/stride_alignment
            %              if dstdir is not also specified
            
            max_shift = 1/8;
            nctrpts = 4;
            behav_type = 'forward';
            template_trange = [];
            verbose = 0;
            dstdir = [];
            
            assignApplicable(varargin);
            
            if isempty(dstdir)
                dstdir = fullfile(fileparts(tsr.behavior.videoFileName),'stride_alignment');
            end
            
            S = tsr.com.strides;
            pp = tsr.com.phasebins;
            nphasebins = length(pp);
            nstrides = length(S);
            
            if isfield(tsr.com,'periphase')
                periphase = tsr.com.periphase;
            else
                periphase = NaN(1,length(tsr.tx));
            end
            
            % get valid cycles
            if isempty(template_trange)
                switch behav_type
                    case 'forward'
                        behav_idx = 1;
                    case 'backward'
                        behav_idx = 2;
                    otherwise
                        behav_idx = 0;
                end
                % TODO should probably use BehaviorStateList instead?
                valid_strides = find([S.behav_label]==behav_idx);
            else
                valid_strides = find(tsr.com.tcross>=template_trange(1) & tsr.com.tcross<=template_trange(2));
                valid_strides(end) = []; % last one is only partially within template_trange
            end
            nvalidstrides = length(valid_strides);
            
            if nvalidstrides==0
                fprintf('no %s stride found!\n',behav_type);
                return;
            end
            
            % get mean
            nmargin = round(nphasebins*max_shift);
            nphasebins_ext = nphasebins+nmargin*2;
            XX_ext = NaN(nvalidstrides,nphasebins_ext-2);
            pp_ext = linspace(0-2*pi*max_shift,2*pi*(1+max_shift),nphasebins_ext-2);
            for j=1:nvalidstrides
                jj = valid_strides(j);
                xx = S(jj).vfwd_tt;

                % Rui Wu approach to margins
                % use the previous/next stride whether or not these are contiguous in time
                % 20231110: only use previous/next stride when continuous in time, otherwise leave as NaN
               

                if jj==1 || S(jj-1).behav_label~=behav_idx
                    xx_prev = NaN(1,nmargin);
                else
                    xx_prev = S(jj-1).vfwd_tt(end-nmargin+1:end);
                end
                if jj==nstrides || S(jj+1).behav_label~=behav_idx
                    xx_next = NaN(1,nmargin);
                else
                    xx_next = S(jj+1).vfwd_tt(1:nmargin);
                end
                XX_ext(j,:) = [xx_prev(1:end-1) xx xx_next(2:end)];

                %MHG approach - use phase on either end of actual stride,
                %even if not a 
%                 p0 = tsr.com.periphase(S(jj).txinds(0));
%                 XX_ext(j,:) = interp1(tsr.com.periphase - p0, tsr.com.vfwd, pp_ext);


            end
            xx0_ext = mean(XX_ext,1,'omitnan');
            
            % align strides
            XX_aligned = NaN(nvalidstrides,nphasebins);
            if verbose>0
                titlestr = sprintf('aligning %s strides, maxshift=1/%d stride length, %d control points',behav_type,1/max_shift,nctrpts);
            end
            if verbose>1
                fig = figure('WindowState','maximized');
            end
            if verbose>2
                if ~isfolder(dstdir)
                    mkdir(dstdir);
                end
                fpath = fullfile(dstdir,sprintf('stride-alignement_%s.mp4',behav_type));
                vw = VideoWriter(fpath,'MPEG-4');
                vw.FrameRate = 1;
                open(vw);
            end
            for j=1:nvalidstrides
                jj = valid_strides(j);
                % 1) get data resampled to phase bins
                xx = S(jj).vfwd_tt;
                tt = S(jj).tt;
                % 2) align
                [pp_warped,delta_p,cpts,scale,offset] = spline_interp_time_warp_with_scale_and_offset(pp_ext,xx0_ext,pp,xx,nctrpts,max_shift);
                tt_aligned = interp1(pp_warped,tt,pp,'linear','extrap');
                p = linspace(0,2*pi,length(S(jj).txinds));
                p_warped = interp1(pp,pp_warped,p);
                % 3) gather result
                S(jj).pp_warped = pp_warped;
                S(jj).tt_aligned = tt_aligned;
                S(jj).delta_p = delta_p;
                S(jj).ctr_pts = cpts;
                S(jj).scale = scale;
                S(jj).offset = offset;
                S(jj).p_warped = p_warped;
                periphase(S(jj).txinds) = p_warped;
                xx_aligned = interp1(pp_warped,xx,pp,'linear',NaN); % plot vs pp/tt to get aligned waveform
                XX_aligned(j,:) = xx_aligned;
                if verbose>1
                    subplot(1,3,1);
                    plot(pp_ext,xx0_ext*scale+offset,'k',pp,xx,'b',pp_warped,xx,'r'); axis square; hold on;
                    xline([0 2*pi],'k:'); hold off;
                    xlim(pp_ext([1 end])); xticks(linspace(0,2*pi,5)); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
                    title('compare'); legend({'adjusted template','data','aligned data'},'Location','northeast');
                    subplot(1,3,2);
                    scatter(pp_ext,(xx0_ext-min(xx0_ext))/range(xx0_ext)-1,'Marker','.','CData',(xx0_ext-min(xx0_ext))/range(xx0_ext)); hold on;
                    scatter(pp,(xx-min(xx))/range(xx),'Marker','.','CData',(xx-min(xx))/range(xx)); axis square;
                    kmap = knnsearch(pp_ext',pp_warped');
                    for k=1:nphasebins
                        kk = kmap(k);
                        plot([pp_ext(kk) pp(k)],[(xx0_ext(kk)-min(xx0_ext))/range(xx0_ext)-1 (xx(k)-min(xx))/range(xx)],'Color',[1 1 1]*0.5);
                    end
                    xline([0 2*pi],'k:'); hold off;
                    xlim(pp_ext([1 end])); xticks(linspace(0,2*pi,5)); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
                    ylim([-1.25 1.25]); yticks([]);
                    title('mapping');
                    subplot(1,3,3);
                    plot([0 2*pi],[0 2*pi],'k',pp,pp_warped,'r.'); axis square tight;
                    xlabel('data'); xlim([0 2*pi]); xticks(linspace(0,2*pi,5)); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
                    ylabel('aligned data'); ylim([0 2*pi]); yticks(linspace(0,2*pi,5)); yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
                    title('warp path');
                    sgtitle({titlestr,sprintf('stride %d (%d/%d aligned, t=%.1f-%.1fs)',jj,j,nvalidstrides,tt([1 end]))});
                    drawnow;
                    if verbose>2
                        writeVideo(vw,getframe(fig));
                    else
                        pause(1);
                    end
                end
            end
            if verbose>2
                close(vw);
            end
            
            % handle -ve values in periphase (not sure if the right thing
            % to do)
            inds = periphase<0;
            periphase(inds) = periphase(inds)+2*pi;
            
            tsr.com.strides = S;
            tsr.com.periphase = periphase;
            
            if verbose>0
                figure;
                subplot(2,1,1);
                plot(pp_ext,XX_ext','Color',[1 1 1]*0.6); hold on;
                h1 = plot(pp_ext,xx0_ext,'k','LineWidth',2);
                xline([0 2*pi],'r--'); hold off;
                legend(h1,'template');
                xlim(pp_ext([1 end])); xticks(linspace(0,2*pi,5)); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
                xlabel('\phi'); title('original');
                subplot(2,1,2);
                plot(pp,XX_aligned','Color',[1 1 1]*0.6); hold on;
                h2(1) = plot(pp_ext,xx0_ext,'k','LineWidth',2);
                h2(2) = plot(pp,mean(XX_aligned,1,'omitnan'),'r--','LineWidth',2);
                xline([0 2*pi],'r--'); hold off;
                legend(h2,{'template','mean(aligned)'});
                xlim(pp_ext([1 end])); xticks(linspace(0,2*pi,5)); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
                xlabel('\phi'); title('aligned');
                sgtitle(titlestr);
            end
        end
        
        saveStrideDetectionMovie(tsr,varargin) % in separate file
        
        function saveImageStacksForAlignedStrides(tsr,varargin)
            %tsr.saveImageStacksForAlignedStrides('cropbox',[xmin ymin xmax ymax]);
            % cropbox centered to tracker in pixels
            
            cropbox = []; % [xmin ymin xmax ymax] centered to tracker; unit: pixels; default: no cropping
            dstdir = fullfile(fileparts(tsr.behavior.videoFileName),'stride_alignment');
            verbose = 0;
            assignApplicable(varargin);
            
            % prep work
            load(tsr.behavior.sleap.filename,'behav_point','behav_segment');
            t_behav = tsr.behavior.et;
            S = tsr.com.strides;
            nstrides = length(S);
            nphasebins = length(tsr.com.phasebins);
            if ~isfolder(dstdir)
                mkdir(dstdir);
            end
            
            % 1) get registered behavior images in original time space
            % prep behavior video reader
            bvr = VideoReader(tsr.behavior.videoFileName);
            nframes = bvr.NumFrames;
            sz = [bvr.Width bvr.Height];
            xl = [-1 1]*sz(1)/2; % following matlab convention when handling image
            yl = [-1 1]*sz(2)/2; % i.e. row=y=height, col=x=width
            % prep for rotation
            theta = zeros(1,nframes);
            for i=1:nframes
                u = behav_segment(2).unit_vector(:,i);
                theta(i) = -atan2(u(2),u(1));
            end
            theta_lp = lowpass1D(theta,1);
            theta_lp = fillmissing(theta_lp,'linear');
            % prep for cropping
            if isempty(cropbox)
                cropbox = [xl(1) yl(1) xl(2) yl(2)];
            end
            xl4 = cropbox([1 3]);
            yl4 = cropbox([2 4]);
            cropsz = [diff(xl4) diff(yl4)];
            croprect = [xl4(1) yl4(1) diff(xl4) diff(yl4)];
            % get stride cubes
            strides_XYT = repmat(struct('t',[],'I',[],'translation',[],'rotation',[]),nstrides,1);
            for j=1:nstrides
                % grab behavior images of this cycle
                t_eval = t_behav(t_behav>=S(j).t(1) & t_behav<=S(j).t(end));
                n = length(t_eval);
                I = zeros(cropsz(1),cropsz(2),n);
                shifts = zeros(2,n);
                thetas = zeros(1,n);
                for i=1:n
                    % get original behavior image
                    [im,frameind] = tsr.getBehaviorImage(bvr,t_eval(i));
                    % pad by 1/4 margin to avoid cropping error later
                    im1 = padarray(im,sz/4,0,'both'); % should be enough
                    xl1 = xl*1.5; yl1 = yl*1.5;
                    imref1 = imref2d(size(im1),xl1,yl1);
                    % recenter to scanning square
                    xc = behav_point(2).pos_px(:,frameind)';
                    dx = sz/2-xc;
                    [im2,imref2] = imtranslate(im1,imref1,dx,'OutputView','full');
                    shifts(:,i) = dx;
                    % rotate by spot-gut axis orientation
                    imT2 = [cos(theta_lp(frameind)) -sin(theta_lp(frameind)) 0; sin(theta_lp(frameind)) cos(theta_lp(frameind)) 0; 0 0 1];
                    [im3,imref3] = imwarp(im2,imref2,affine2d(imT2));
                    thetas(i) = theta_lp(frameind);
                    % crop to larva body size
                    im4 = imcrop(imref3.XWorldLimits,imref3.YWorldLimits,im3,croprect);
                    % deal with imcrop() not always exact
                    sz4 = size(im4);
                    padsz = [0 0];
                    if sz4(1)>cropsz(2)
                        im4 = im4(1:cropsz(2),:);
                    elseif sz4(1)<cropsz(2)
                        padsz(1) = cropsz(2)-sz4(1);
                    end
                    if sz4(2)>cropsz(1)
                        im4 = im4(:,1:cropsz(1));
                    elseif sz4(2)<cropsz(1)
                        padsz(2) = cropsz(1)-sz4(2);
                    end
                    if any(padsz)
                        im4 = padarray(im4,padsz,0,'post');
                    end
                    I(:,:,i) = im4'; % return to our convention: first dimension = x
                end
                strides_XYT(j).t = t_eval;
                strides_XYT(j).I = I;
                strides_XYT(j).translation = shifts;
                strides_XYT(j).rotation = thetas;
            end
            % save to disk
            im_ref = imref2d(cropsz',xl4,yl4);
            fprintf('saving original behavior image stacks to %s...\n',dstdir);
            fpath = fullfile(dstdir,'strides_XYT.mat');
            save(fpath,'strides_XYT','im_ref','-v7.3');
            
            % 2) get images in original phase space
            strides_XYP = struct('t',[],'I',[]);
            I = zeros([nstrides im_ref.ImageSize' nphasebins]);
            T = NaN(nstrides,nphasebins);
            for j=1:nstrides
                tt = S(j).tt;
                for i=1:nphasebins
                    k = find(strides_XYT(j).t<=tt(i),1,'last');
                    if ~isempty(k)
                        % check if k is last frame in this stride and too far from tt
                        if k==length(strides_XYT(j).t)
                            if (tt(i)-strides_XYT(j).t(k)) > max(diff(strides_XYT(j).t)) % look in next stride
                                k = find(strides_XYT(j+1).t<=tt(i),1,'last');
                                if ~isempty(k)
                                    if verbose>0
                                        fprintf('stride %d bin %d located in next stride\n',j,i);
                                    end
                                    I(j,:,:,i) = strides_XYT(j+1).I(:,:,k);
                                    T(j,i) = strides_XYT(j+1).t(k);
                                end
                            else % fine to use last frame in this stride
                                I(j,:,:,i) = strides_XYT(j).I(:,:,k);
                                T(j,i) = strides_XYT(j).t(k);
                            end
                        else
                            I(j,:,:,i) = strides_XYT(j).I(:,:,k);
                            T(j,i) = strides_XYT(j).t(k);
                        end
                    else
                        if strides_XYT(j).t(1)>tt(i) && j>1 % look in previous stride
                            k = find(strides_XYT(j-1).t<=tt(i),1,'last');
                            if ~isempty(k)
                                if verbose>0
                                    fprintf('stride %d bin %d located in previous stride\n',j,i);
                                end
                                I(j,:,:,i) = strides_XYT(j-1).I(:,:,k);
                                T(j,i) = strides_XYT(j-1).t(k);
                            end
                        else
                            fprintf('ATTENTION: stride %d bin %d cannot locate behavior frame; this shouldn''t happen unless it''s the first or last stride\n',j,i);
                        end
                    end
                end
            end
            strides_XYP.t = T;
            strides_XYP.I = I;
            % save to disk
            fprintf('saving resampled behavior image stacks to %s...\n',dstdir);
            fpath = fullfile(dstdir,'strides_XYP.mat');
            save(fpath,'strides_XYP','im_ref','-v7.3');
            
            % 3) get images in aligned phase space
            strides_XYP_aligned = repmat(struct('t',[],'I',[],'behav_label',0,'stride_inds',[]),2,1);
            for q=1:2
                strides_XYP_aligned(q).behav_label = q;
                valid_strides = find([S.behav_label]==q);
                nvalidstrides = length(valid_strides);
                strides_XYP_aligned(q).stride_inds = valid_strides;
                I = zeros([nvalidstrides im_ref.ImageSize' nphasebins]);
                T = NaN(nvalidstrides,nphasebins);
                for j=1:nvalidstrides
                    jj = valid_strides(j);
                    tt_aligned = S(jj).tt_aligned;
                    for i=1:nphasebins
                        k = find(strides_XYT(jj).t<=tt_aligned(i),1,'last');
                        if ~isempty(k)
                            % check if k is last frame in this stride and too far from tt
                            if k==length(strides_XYT(jj).t)
                                if (tt_aligned(i)-strides_XYT(jj).t(k)) > max(diff(strides_XYT(jj).t)) % look in next stride
                                    k = find(strides_XYT(jj+1).t<=tt_aligned(i),1,'last');
                                    if ~isempty(k)
                                        if verbose>0
                                            fprintf('stride %d bin %d located in next stride\n',jj,i);
                                        end
                                        I(j,:,:,i) = strides_XYT(jj+1).I(:,:,k);
                                        T(j,i) = strides_XYT(jj+1).t(k);
                                    end
                                else % fine to use last frame in this stride
                                    I(j,:,:,i) = strides_XYT(jj).I(:,:,k);
                                    T(j,i) = strides_XYT(jj).t(k);
                                end
                            else
                                I(j,:,:,i) = strides_XYT(jj).I(:,:,k);
                                T(j,i) = strides_XYT(jj).t(k);
                            end
                        else
                            if strides_XYT(jj).t(1)>tt_aligned(i) && jj>1 % look in previous stride
                                k = find(strides_XYT(jj-1).t<=tt_aligned(i),1,'last');
                                if ~isempty(k)
                                    if verbose>0
                                        fprintf('stride %d bin %d located in previous stride\n',jj,i);
                                    end
                                    I(j,:,:,i) = strides_XYT(jj-1).I(:,:,k);
                                    T(j,i) = strides_XYT(jj-1).t(k);
                                end
                            else
                                fprintf('ATTENTION: stride %d bin %d cannot locate behavior frame; this shouldn''t happen unless it''s the first or last stride\n',jj,i);
                            end
                        end
                    end
                end
                strides_XYP_aligned(q).t = T;
                strides_XYP_aligned(q).I = I;
            end
            % save to disk
            fprintf('saving aligned behavior image stacks to %s...\n',dstdir);
            fpath = fullfile(dstdir,'strides_XYP_aligned.mat');
            save(fpath,'strides_XYP_aligned','im_ref','-v7.3');
            
            disp('...done!');
        end
    end
    
    methods (Static)
        function ks = getKalmanScalings()
            ks.ms_timer = 1e-3; %ms to s
            ks.timer = 256/1.6E8; %ticks/256 to s
            ks.loc = 1/64; %16/16 to 16/10
            ks.vel_com = 8E7/2^22; %16/16 to 16/-6 microns, per tick to per second
            ks.sampling_com = 1/64; %16/16 to 16/10
            ks.piezo_val = 1/64;
            ks.weighted_sum = 0.5;
            ks.integral_term = 2^-16;% (16/16 to 0/16)
            ks.Pxx = 2^-25; %32/32 to 32/7
            ks.Pvv = 2^-58 * (8E7)^2; %32/32 to 32/-26, per tick squared to per second squared
            
        end
        function d = importDataToStruct(fname, stripduplicates, varargin)
            % custom loading function to handle large _tracker.txt files - Rui 2020/8/31
            fsize = dir(fname);
            fsize = fsize.bytes;
            if ispc
                maxarraysize = memory;
                maxarraysize = maxarraysize.MaxPossibleArrayBytes;
            elseif isunix
                [~,w] = unix('free | grep Mem');
                stats = str2double(regexp(w,'[0-9]*','match'));
                memAvailable = (stats(3)+stats(end))*1024;
                maxarraysize = memAvailable;
            end
            if maxarraysize/fsize>10 % not sure if too daring
                data = importdata(fname, '\t');
            else
                data = importLargeTextFile(fname,varargin);
            end
            %for some reason, data files from tracker to have duplicate
            %lines as of 4/8/19 -- this removes them
            %note that we use columns 2:end, because the first column is
            %mstimer, which can change even if data is duplicated
            if (nargin < 2)
                stripduplicates = false;
            end
            if (stripduplicates)
                valid = [true; ~(all(diff(data.data(:,2:end)) == 0,2))];
                data.data = data.data(valid,:);
            end
            if ~isfield(data,'colheaders')
            % ----- Check if colheaders field exists. If not, load it differently. -----
                warning('applying column headers manually to fix a specific bug in October 2021');
                data.colheaders = strsplit('ms_timer timer axis status stim_level opto_level loc velocity piezo_val photons_ch0 photons_ch1 weighted_sum Pii num_phot_detections_from_pmt num_phot_detections_from_pmt_in_spot P00 sampling_com',' ');
            end
            for j = 1:length(data.colheaders)
                d.(data.colheaders{j}) = data.data(:,j)';
            end
        end
        function raw = removeInvalidValues(raw)
            % WARNING: this is only a temporary fix and we shouldn't rely
            % on this for actual data processing
            npts = length(raw.timer);
            % raw.timer should not have any repetition: it's FPGA tick*256
            if (nnz(diff(raw.timer) == 0 & diff(raw.axis) ~= 0) < 10) 
                %6/7/2021 newer tracking method outputs 3 axes silmutaneously 
                %so turn off raw.timer detection if multiple places show
                %same time with different axes
                [~,inds2,~] = unique(raw.timer);
                if length(inds2)~=npts
                    warning('raw.timer has duplicates');
                end
            else
                inds2 = 1:length(raw.timer);
            end
            % max(raw.axis) should be reasonable: say we can track 10
            % neurons at a time -> max(raw.axis)=10*3
            maxnspots = 100;
            inds1 = find(raw.axis<=maxnspots*3);
            if length(inds1)~=npts
                warning('raw.axis has invalid values');
            end
            % combine all valid inds
            inds = intersect(inds1,inds2);
            fn = fieldnames(raw);
            for i=1:length(fn)
                raw.(fn{i}) = raw.(fn{i})(inds);
            end
        end
        function raw = processStatus(raw)
            if (isstruct(raw.status))
                return;
            end
            sb = uint8(raw.status);
            fn = {'tracking', 'feedback0', 'feedback1', 'TTL0', 'TTL1', 'TTL2', 'TTL3'};
            for j = 1:length(fn)
                status.(fn{j}) = logical(bitand(uint8(bitshift(1, j-1)), sb, 'uint8'));
            end
            raw.status = status;
        end
        function [ratio, sigma_ratio, lambda, sigma_lambda, thetaOut, cmat_thetaOut] = MLERates (tedge, nr, ng, D1, D2)
            %function [lamda, wlambda, theta, wtheta] = MLERates (t, nr, ng, D1, D2)
            %tedge is 1 element longer than nr, ng
            %nr(j) is the number of red photons counted between tedge(j)
            %and tedge(j+1)
            %
            t = tedge(:);
            nr = nr(:);
            ng = ng(:);
            deltat = diff(t); %deltat_k = t_k - t_k-1; deltat(2:end)_k = t_k+1 - t_k
            if (any(deltat <= 0))
                error ('time must increase at each time step');
            end
            if (any(nr < 0) || any(ng < 0))
                error ('number of photons must be >= 0');
            end
            problem.objective = @probfun;
            lg = log(max(ng,1)./deltat);
            lr = log(max(nr,1)./deltat);
            problem.x0 = [lg-lr;lr];
            problem.solver = 'fminunc';
            problem.options = optimoptions(problem.solver);
            problem.options.GradObj = 'on';
            problem.options.Hessian = 'on';
            problem.options.Algorithm = 'trust-region';
            
            thetaOut = reshape(fminunc(problem), [], 2);
            
            ratio = exp(thetaOut(:,1));
            lambda = exp(thetaOut*[1 0;1 1]);
            
            %negative hessian = inverse covariance matrix
            [~,~,nlh] = probfun(thetaOut(:));
            nn = length(thetaOut);
            d1 = spdiags(nlh,0);
            d2 = spdiags(nlh,-nn);
            d2 = d2(1:nn);
            dtrm = d1(1:nn).*d1((nn+1):end) - d2.^2;
            cmat_thetaOut(:,1) = d1((nn+1):end)./dtrm;
            cmat_thetaOut(:,2) = d2./dtrm;
            cmat_thetaOut(:,3) = d1(1:nn)./dtrm;
            sigma_ratio = sqrt(ratio.*cmat_thetaOut(:,1));
            sigma_lambda = lambda.*(cmat_thetaOut*[1 0;2 0; 1 1]).^(0.5);
            
            
            function [negLogP, negLogPgrad, negLogPHess] = probfun (theta)
                %lambdaR = exp(theta_2); lambdaG = exp(theta_2 + theta_1)
                %
                %logP(t) = -(lambdaR(t) dT) + num_r(t) * theta_2(t) ...
                % -(lambdaG(t) dT) + num_g(t) * (theta_1(t) + theta_2(t))
                % - delta_theta_1^2 / 4D1t - delta_theta_2^2 / 4D2t
                %
                %first two lines are poisson probability of data given
                %rates
                %last line is continuity prior
                
                theta = reshape(theta, [],2);
                lambdaR = exp(theta(:,2));
                lambdaG = exp(theta(:,1) + theta(:,2));
                %deltat_k = t_k - t_k-1; deltat(2:end)_k = t_k+1 - t_k
                
                dtheta = [[0 0]; diff(theta)];
                
                logP = sum(-(lambdaR+lambdaG).*deltat + ng.*theta(:,1) + (ng + nr).*theta(:,2)) ...
                    - sum(dtheta(:,1).^2./(4*D1*deltat) + dtheta(:,2).^2./(4*D2*deltat));
                
                %deltat_k = t_k - t_k-1; deltat(2:end)_k = t_k+1 - t_k
                logPgrad(:,1) = ng - deltat.*lambdaG + 1./(2*D1)*(dtheta([2:end 1],1)./deltat([2:end 1]) - dtheta(:,1)./deltat);
                logPgrad(:,2) = ng + nr - deltat.*(lambdaR + lambdaG) + 1./(2*D2)*(dtheta([2:end 1],2)./deltat([2:end 1]) - dtheta(:,2)./deltat);
                
                H12 = -lambdaG.*deltat;
                H11 = -lambdaG.*deltat - 1./(2*D1).*([1./deltat(2:end); 0] + 1./deltat);
                H22 = -(lambdaR+lambdaG).*deltat - 1./(2*D2).*([1./deltat(2:end); 0] + 1./deltat);
                n = size(theta, 1);
                
                negLogP = -logP;
                negLogPgrad = -logPgrad(:);
                negLogPHess = spdiags(-[[H12;H12], [H11;H22], [H12;H12]], [-n 0 n], 2*n, 2*n);
            end
           
        end
    end
end