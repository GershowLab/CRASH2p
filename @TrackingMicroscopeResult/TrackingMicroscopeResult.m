classdef TrackingMicroscopeResult
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tracker;
        stage;
        neuron;
        video;
        posture;
        dt = 1e-3;
        d_loglambda = [.01 .1]; %ratio diffusion constant, total fluorescence diffusion constant
        tx;
        bandwidth_3db = 5;
        rate_fit;
        hasVideo = false;
        micronsPerPixel = 8.6;% found 4/26 by tracking flow on bubbly slide 10.65;%found by manual registration of points in video vs. stage movement; 15.2563; %found by measuring size of 5.08 mm objective glass, which is slightly out of focus
        vidToStageTransform;
        stageToVidTransform;
        zflip = +1; %-1 if +z is away from objective (down); used for plotting only
        bsl; %behavioral state list
    end
    properties (Transient = true)
        videocube=[];
        videogradients=[];
    end
    
    methods
        function tmr = TrackingMicroscopeResult(varargin)
            if (nargin > 0)
                tic
                tmr = tmr.loadFile(varargin{:});
                toc; disp ('loaded data');
                try 
                    tmr = tmr.conditionData();
                    tmr = tmr.processPhotonCounts();
                    tmr = tmr.correctRatioForBleaching();
                    toc; disp ('processed data and photon counts');
                   
                catch me
                    disp (me.getReport());
                end
            end
            tmr = tmr.initBSL();
            
        end
    
        function tmr = initBSL(tmr)
            tmr.bsl = BehavioralStateList('forward', 'backward', 'left turn', 'right turn', 'squished', 'hunch', 'other');
            if (isfield(tmr.video, 'et') && ~isempty(tmr.video.et))
                tmr.bsl = tmr.bsl.setLength(length(tmr.video.et));
            end
        end
        
        function tmr = behaviorGUI(tmr)
            try
                addpath (fullfile(pwd, 'TMRgui'));
            catch
            end
            tmr = annotateBehavioralStateGUI(tmr);
        end
        
        function tmr = addVideo(tmr, varargin)
            tic;
            tmr = tmr.findVideo();
            if (~tmr.hasVideo)
                disp ('no video found');
                return;
            end
            tmr = tmr.loadVideo();
            if (isfield(tmr.video, 'et') && ~isempty(tmr.video.et))
                tmr.bsl = tmr.bsl.setLength(length(tmr.video.et));
            end
            if (~tmr.hasVideo)
                disp ('no video loaded');
                return;
            end
            toc; disp ('loaded video');
            tmr = tmr.alignVideo(varargin{:});
            tmr = tmr.trackSpot();
            tmr = tmr.alignVideoToSpot();
            toc; disp ('tracked beam spot');
        end
        function tmr = enforceVideoDeltaT(tmr, enforcedDeltaT)
            tmr = tmr.alignVideo('enforcedDeltaT',enforcedDeltaT);
            tmr = tmr.alignVideoToSpot();
        end
        
        function tmr = loadFile (tmr, filename, filename2)
            if (nargin > 1 && ischar(filename))
                if (exist (filename, 'dir'))
                    d1 = dir(fullfile(filename, '*-tracking.*'));
                    d2 = dir(fullfile(filename, '*-stage.*'));
                    tmr = tmr.loadFile(fullfile(filename, d1.name), fullfile(filename, d2.name));
                    return;
                end
                tmr.tracker.data = importdata(filename);
                tmr.tracker.filename = filename;
            end
            
            
            
            if (size(tmr.tracker.data,1) > size(tmr.tracker.data,2))
                tmr.tracker.data = tmr.tracker.data';
                
            end
            
            if (nargin > 1)
                tmr.stage.filename = filename2;
                tmr.stage.data = importdata(filename2)';
                tmr.stage.mstimer = tmr.stage.data(1,:);
                tmr.stage.loc = tmr.stage.data(2:4,:);
            end
            
            raw.ax = tmr.tracker.data(1,:);
            raw.pos = tmr.tracker.data(2,:);
            raw.vel = tmr.tracker.data(3,:);
            raw.perr = tmr.tracker.data(4,:);
            raw.status = tmr.tracker.data(5,:);
            raw.tracking = bitand(2, raw.status);
            tmr.tracker.hasStatus = (any(raw.status));
            
            
            
%             raw.verr = tmr.tracker.data(5,:);
            raw.gphoton = tmr.tracker.data(6,:);
            raw.rphoton = tmr.tracker.data(7,:);
            raw.et = tmr.tracker.data(8,:);
            raw.mstimer = tmr.tracker.data(9,:);
            
            raw.nspots = ceil(max(raw.ax)/3);
            raw.spotnum = ceil((raw.ax+1)/3);
            for j = 1:raw.nspots           
                %eliminate problem indices where et decreases
                raw.spot(j).xinds = (raw.ax == 0 + 3*(j-1)) & [true diff(raw.et) > 0];
                raw.spot(j).yinds = (raw.ax == 1 + 3*(j-1)) & [true diff(raw.et) > 0];
                raw.spot(j).zinds = (raw.ax == 2 + 3*(j-1)) & [true diff(raw.et) > 0];
                 %assumes feedback on both channels!
                raw.spot(j).xinn = raw.perr(raw.spot(j).xinds)./(raw.gphoton(raw.spot(j).xinds) + raw.rphoton(raw.spot(j).xinds));
                raw.spot(j).yinn = raw.perr(raw.spot(j).yinds)./(raw.gphoton(raw.spot(j).yinds) + raw.rphoton(raw.spot(j).yinds));
                raw.spot(j).xt = raw.et(raw.spot(j).xinds);
                raw.spot(j).yt = raw.et(raw.spot(j).yinds);
            end
            for j = 1:raw.nspots
                anytracking(j) = any(raw.tracking(raw.spot(j).xinds));
            end
            if (any(anytracking))
                raw.spot = raw.spot(anytracking);
                raw.nspots = length(raw.spot);
            else
                warning ('no tracking detected');
            end
            
            
            if (tmr.tracker.hasStatus && ~any(raw.tracking))
                disp ('according to status bit, not tracking -- going with backup plan');
                et = raw.et(raw.ax == 1);
                tracking = diff(et) < 0.001;
                si = find(diff(tracking) > 0);
                ei = find(diff(tracking) < 0);
                if (tracking(1))
                    si = [1 si];
                end
                if (tracking(end))
                    ei = [ei length(et)];
                end
                [~,I] = max(ei - si);
                disp (['detected tracking time is ' num2str(ei(I) - si(I))]);
                raw.tracking = raw.et > et(si(I)) & raw.et < et(ei(I));
            end
            
            
            tmr.tracker.raw = raw;
        end
        
        function tmr = conditionData (tmr)
            raw = tmr.tracker.raw;
            et = raw.et;
            pos = raw.pos;
            
            tmr.tx = min(et):tmr.dt:max(et);
            
            for j = 1:raw.nspots
                xinds = raw.spot(j).xinds;
                yinds = raw.spot(j).yinds;
                zinds = raw.spot(j).zinds;
                x = resample_to_t(et(xinds), pos(xinds), tmr.tx);
                y = resample_to_t(et(yinds), pos(yinds), tmr.tx);
                z = resample_to_t(et(zinds), pos(zinds), tmr.tx);
                tmr.tracker.loc(:,:,j) = [x;y;z];
                tmr.tracker.loc_filt(:,:,j) = tmr.lpfilt(tmr.tracker.loc(:,:,j));
                tmr.tracker.offset(:,:,j) = [resample_to_t(raw.spot(j).xt(isfinite(raw.spot(j).xinn)), raw.spot(j).xinn(isfinite(raw.spot(j).xinn)), tmr.tx); ...
                    resample_to_t(raw.spot(j).yt(isfinite(raw.spot(j).yinn)), raw.spot(j).yinn(isfinite(raw.spot(j).yinn)),tmr.tx)];
                vel = raw.vel;
                vx = resample_to_t(et(xinds), vel(xinds), tmr.tx);
                vy = resample_to_t(et(yinds), vel(yinds), tmr.tx);
                vz = resample_to_t(et(zinds), vel(zinds), tmr.tx);
                tmr.tracker.vel(:,:,j) = [vx;vy;vz];
                tmr.tracker.vel_filt(:,:,j) = tmr.lpfilt(tmr.tracker.vel(:,:,j));
            end
            zinds = any(cat(1,raw.spot.zinds),1);
            pz = resample_to_t(et(zinds), raw.perr(zinds), tmr.tx);
            tmr.tracker.pz = pz;
            
            
            
            if (tmr.tracker.hasStatus)
                xinds = any(cat(1,raw.spot.xinds),1);
                tmr.tracker.status.scanning = logical(interp1(et(xinds), bitand(1, raw.status(xinds)), tmr.tx, 'nearest', 'extrap'));
                tmr.tracker.status.tracking = logical(interp1(et(xinds), bitand(2, raw.status(xinds)), tmr.tx, 'nearest', 'extrap'));
                tmr.tracker.status.outOfRange = logical(interp1(et(xinds), bitand(4, raw.status(xinds)), tmr.tx, 'nearest', 'extrap'));
                tmr.tracker.status.stimOn = logical(interp1(et(xinds), bitand(8, raw.status(xinds)), tmr.tx, 'nearest', 'extrap'));
                si = find(diff([0 tmr.tracker.status.tracking]) > 0);
                ei = find(diff([tmr.tracker.status.tracking 0]) < 0);
                [~,I] = max(ei - si);
                tmr.tracker.startTime = tmr.tx(si(I));
                tmr.tracker.endTime = tmr.tx(ei(I));             
            else
                %guess tracking status by delta_t
                maxdt = 0.5;
                dxt = diff(tmr.tracker.raw.spot(1).xt);
                tmr.tracker.status.tracking = [false (dxt < mean(dxt(dxt < maxdt)) + std(dxt(dxt < maxdt))*6)]; %discard crazy long intervals
                tmr.tracker.status.tracking(1) = tmr.tracker.status.tracking(2);
                
                start = diff(tmr.tracker.status.tracking) > 0;
                stop = diff(tmr.tracker.status.tracking) < 0;
                if (tmr.tracker.status.tracking(1))
                    start(1) = true;
                end
                
                if (tmr.tracker.status.tracking(end))
                    stop(end) = true;
                end
                
                [~,II] = max(tmr.tracker.raw.spot(1).xt(stop)-tmr.tracker.raw.spot(1).xt(start));
                start = find(start); stop = find(stop);
                
                tmr.tracker.startTime = tmr.tracker.raw.spot(1).xt(start(II));
                tmr.tracker.endTime = tmr.tracker.raw.spot(1).xt(stop(II));
                tmr.tracker.status.tracking = tmr.tx >= tmr.tracker.startTime & tmr.tx < tmr.tracker.endTime;
            end
            
            
            
            
            mstimer = raw.mstimer;
            inds = (find(diff(mstimer) > 0));
            mstimer_adj = interp1(et(inds), mstimer(inds), tmr.tx, 'linear', 'extrap');
            %             stagex = interp1(tmr.mstimer_stage, tmr.stagex, mstimer_adj);
            %             stagey = interp1(tmr.mstimer_stage, tmr.stagey, mstimer_adj);
            %             stagez = interp1(tmr.mstimer_stage, tmr.stagez, mstimer_adj);
            tmr.tracker.stageloc = interp1 (tmr.stage.mstimer, tmr.stage.loc', mstimer_adj)';
            for j = 1:raw.nspots
                tmr.neuron(j).loc = tmr.tracker.stageloc + squeeze(tmr.tracker.loc(:,:,j));
                tmr.neuron(j).loc_filt = tmr.lpfilt(tmr.neuron(j).loc);
                tmr.neuron(j).vel = 1/tmr.dt * deriv(tmr.neuron(j).loc, 1);
                tmr.neuron(j).vel_filt = 1/tmr.dt * deriv(tmr.neuron(j).loc_filt, 1);
                tmr.neuron(j).speed = sqrt(sum(tmr.neuron(j).vel.^2));
                tmr.neuron(j).speed_filt = sqrt(sum(tmr.neuron(j).vel_filt.^2));
            end
        end
        
        function tmr = processPhotonCounts (tmr, varargin)
            raw = tmr.tracker.raw;
            sampletime = 0.005;
            varargin = assignApplicable(varargin);
            txbase = min(tmr.tx):sampletime:max(tmr.tx);
            for j = 1:raw.nspots
                
                pinds = (raw.spot(j).xinds & raw.tracking) | (raw.spot(j).yinds &~raw.tracking);
                et = raw.et(pinds);
                rp = raw.rphoton(pinds);
                gp = raw.gphoton(pinds);
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
                [ratio, sigma_ratio, lambda, sigma_lambda, thetaOut, cmat_thetaOut] = TrackingMicroscopeResult.MLERates(tedges, nr, ng, tmr.d_loglambda(1), tmr.d_loglambda(end));
                ratefit.et = tedges(2:end);
                ratefit.theta = thetaOut;
                ratefit.wtheta = cmat_thetaOut;
                ratefit.ratio = ratio;
                ratefit.sigma_ratio = sigma_ratio;
                ratefit.lambda = lambda;
                ratefit.sigma_lambda = sigma_lambda;
                tmr.neuron(j).redrate = interp1(ratefit.et,ratefit.lambda(:,2), tmr.tx, 'linear');
                tmr.neuron(j).redrate_eb = interp1(ratefit.et,ratefit.sigma_lambda(:,2), tmr.tx, 'linear');
                tmr.neuron(j).greenrate = interp1(ratefit.et,ratefit.lambda(:,1), tmr.tx, 'linear');
                tmr.neuron(j).greenrate_eb = interp1(ratefit.et,ratefit.sigma_lambda(:,1), tmr.tx, 'linear');
                tmr.neuron(j).ratio = interp1(ratefit.et,ratefit.ratio, tmr.tx, 'linear');
                tmr.neuron(j).ratio_eb = interp1(ratefit.et,ratefit.sigma_ratio, tmr.tx, 'linear');
                if (j == 1)
                    tmr.rate_fit = ratefit;
                else
                    tmr.rate_fit(j) = ratefit;
                end
            end
        end
%
        function tmr = processPhotonCountsDeprecated (tmr, presmooth)
            existsAndDefault('presmooth', 0);
            raw = tmr.tracker.raw;
            trackingbuffer = 2*ceil(20e-3/tmr.dt);
            tracking = logical(interp1(tmr.tx, double(imerode(tmr.tracker.status.tracking,ones([1 trackingbuffer]))), raw.et, 'nearest', 'extrap'));
            anyy = any(cat(1,raw.spot.yinds), 1);
            dto = diff(raw.et(anyy)); 
            dto = [dto(1) dto]; 
            dt_obs_all = ones(size(raw.et))*median(diff(raw.et));
            dt_obs_all(anyy & ~tracking) = dto(~tracking(anyy));
            anyx = any(cat(1,raw.spot.xinds), 1);
            dto = diff(raw.et(anyx)); 
            dto = [dto dto(end)]; 
            dt_obs_all(anyx & tracking) = dto(tracking(anyx));
            
            for j = 1:raw.nspots
                
                pinds = (raw.spot(j).yinds &~tracking) | (raw.spot(j).xinds & tracking);
                et = raw.et(pinds);
                dt_obs = dt_obs_all(pinds);
                track = tracking(pinds);

                dtnotrack = 0.2;
                rp = raw.rphoton(pinds);
                gp = raw.gphoton(pinds);
                [tauobs, tf] = resample_to_t(et, dt_obs, min(et):dtnotrack:max(et));
                rpf = resample_to_t(et, rp, min(et):dtnotrack:max(et))./tauobs * dtnotrack;
                
                gpf = resample_to_t(et, raw.gphoton(pinds), min(et):dtnotrack:max(et))./tauobs * dtnotrack;
                not_track_f = ~logical(interp1(tmr.tx, double(tmr.tracker.status.tracking), tf, 'nearest', 'extrap'));
               

                et = et(track);
                gp = gp(track);
                rp = rp(track);
                dt_obs = dt_obs(track);
                rate_fit_notrack = photonSPPFtworates([tf(not_track_f);dtnotrack*ones(size(tf(not_track_f)))], rpf(not_track_f), gpf(not_track_f), tmr.d_loglambda, min(2,presmooth));
                if (j == 1)
                    tmr.rate_fit = photonSPPFtworates([et;dt_obs], rp, gp, tmr.d_loglambda, presmooth);
                else
                    tmr.rate_fit(j) = photonSPPFtworates([et;dt_obs], rp, gp, tmr.d_loglambda, presmooth);
                end
                et = [rate_fit_notrack.et tmr.rate_fit(j).et];
                lr = [rate_fit_notrack.lambda_red_s tmr.rate_fit(j).lambda_red_s];
                slr = [rate_fit_notrack.slambda_red_s tmr.rate_fit(j).slambda_red_s];
                lg = [rate_fit_notrack.lambda_green_s tmr.rate_fit(j).lambda_green_s];
                slg = [rate_fit_notrack.slambda_green_s tmr.rate_fit(j).slambda_green_s];
                rat = [rate_fit_notrack.ratio_s tmr.rate_fit(j).ratio_s];
                srat= [rate_fit_notrack.sratio_s tmr.rate_fit(j).sratio_s];
                [~,I] = sort(et);
                lr = lr(I);
                slr = slr(I);
                lg = lg(I);
                slg = slg(I);
                rat = rat(I);
                srat = srat(I);
                tmr.neuron(j).redrate = resample_to_t(et,lr, tmr.tx);
                tmr.neuron(j).redrate_eb = resample_to_t(et,slr, tmr.tx);
                tmr.neuron(j).greenrate = resample_to_t(et,lg, tmr.tx);
                tmr.neuron(j).greenrate_eb = resample_to_t(et,slg, tmr.tx);
                tmr.neuron(j).ratio = resample_to_t(et,rat, tmr.tx);
                tmr.neuron(j).ratio_eb = resample_to_t(et,srat, tmr.tx);
            end
        end
        function fdata = lpfilt (tmr, data)
            %function fdata = lpfilt (tmr, data)
            fdata = lowpass1D(data, 1/(7.3*tmr.bandwidth_3db*tmr.dt));
        end
        
        function tmr = findVideo(tmr, filename)
            if (nargin < 2 || isempty(filename))
                dd = fileparts(tmr.tracker.filename);
                d = dir(fullfile(dd, 'acA640*.avi'));
                if (isempty(d))
                    d = dir(fullfile(dd, '*.avi'));
                    if (length(d) ~= 1)
                        disp('could not identify movie file automatically.');
                        disp('here are some candidates');
                        disp({d.name});
                        tmr.hasVideo = false;
                        return;
                    end
                end
                filename = fullfile(dd, d.name);
                if (~isempty(strfind(lower(d.name), 'aca2040')))
                    disp (['only movie file found: ' d.name ' looks like an image stack instead of a track movie']);
                    disp (['if this is correct, type: tmr.video.filename = ' filename]);
                    disp ('where tmr is the name of your object, and continue with tmr.loadVideo()');
                    tmr.hasVideo = false;
                    return;
                end
            end
            tmr.hasVideo = true;
            tmr.video.filename = filename;
        end
        
        function tmr = loadVideo(tmr, filename)
         % function tmr = loadVideo(tmr, filename)
            existsAndDefault('filename', '');
            if (~isfield(tmr.video, 'filename') || isempty(tmr.video.filename))
                tmr = tmr.findVideo(filename);
                if (~tmr.hasVideo)
                    return;
                end
            end
            try
                v = VideoReader(tmr.video.filename);
                if (~isfield(tmr.video, 'et') || isempty(tmr.video.et) || ...
                        ~isfield(tmr.video, 'intensity') || isempty(tmr.video.intensity) || ...
                        ~isfield(tmr.video, 'startFrame') || isempty(tmr.video.startFrame) || tmr.video.startFrame <= 0 ||...
                        ~isfield(tmr.video, 'endFrame') || isempty(tmr.video.endFrame) || tmr.video.endFrame <= 0)
                    ind = 0;
                    while (v.hasFrame())
                        ind = ind + 1;
                        tmr.video.et(ind) = v.currentTime();
                        fr = v.readFrame();
                        if (size(fr, 3) > 1)
                            fr = rgb2gray(fr);
                        end
                        tmr.video.intensity(ind) = sum(sum(fr));
                    end
                end
                if (v.hasFrame())
                    fr = v.readFrame();
                end
                tmr.videocube = zeros(v.Height, v.Width, length(tmr.video.et));
                v.currentTime = tmr.video.et(1);
                ind = 0;
                while (v.hasFrame())
                    ind = ind+1;
                    fr = v.readFrame();
                    if (size(fr, 3) > 1)
                        fr = rgb2gray(fr);
                    end
                    tmr.videocube(:,:,ind) = fr;
                end
            catch me
                disp (me.getReport());
                tmr.hasVideo = false;
            end
            tmr.hasVideo = true;         
        end
        
        %aligns video time scale to tracker time scale using intensity
        function tmr = alignVideo(tmr, varargin)
            enforcedDeltaT = [];
            varargin = assignApplicable(varargin);
            if (~isfield(tmr.tracker, 'startTime'))
                dxt = diff(tmr.tracker.raw.spot(1).xt);
                si = find(diff(dxt > mean(dxt) + std(dxt)*6) < 0, 1, 'first');
                if (isempty(si))
                    si = 1;
                end
                ei = si-1 + find(diff(dxt(si:end) > mean(dxt) + std(dxt)*6) > 0, 1, 'first');
                if (isempty(ei))
                    ei = length(dxt);
                end
                tmr.tracker.startTime = tmr.tracker.raw.spot(1).xt(si);

                tmr.tracker.endTime = tmr.tracker.raw.spot(1).xt(ei);              
            end
            
            ddi = conv2(tmr.video.intensity, [-1 16 -30 16 -1]/12, 'valid');
            [~,I] = sort(ddi, 'descend');
            deltaI = floor((length(tmr.video.intensity) - length(ddi))/2);
            I = I + deltaI;
            
            %find first point separate by at least "5" seconds, scare
            %quotes because video time may not be reliable
            
            I2 = I(find(abs(tmr.video.et(I) - tmr.video.et(I(1))) > 5, 1, 'first'));
            I1 = min(I(1), I2);
            I2 = max(I(1), I2);
            
            tmr = tmr.setVideoTrackingAlignment(I1, I2, enforcedDeltaT);            
        end
        
        function tmr = setVideoTrackingAlignment (tmr, startFrame, endFrame, enforcedDeltaT)
            existsAndDefault('enforcedDeltaT', []);
            if (tmr.tracker.hasStatus && ~any(tmr.tracker.status.tracking))
                warning ('all status = not tracking');
                return;
            end
            tmr.video.startFrame =  startFrame;
            if (~isempty(enforcedDeltaT) && enforcedDeltaT > 0)
                if (nargin < 2) || isempty(endFrame) || endFrame < startFrame
                    error ('enforcedDeltaT requires start and end frames');
                end
                tmr.video.tx = ((1:length(tmr.video.et))-startFrame)*enforcedDeltaT + tmr.tracker.startTime;
                tmr.video.endFrame = find(tmr.video.tx <= tmr.tracker.endTime,1,'last');
                tmr.video.tx(tmr.video.tx > max(tmr.tx) | tmr.video.tx < min(tmr.tx)) = NaN;
                tmr.video.timewarp = (tmr.tracker.endTime-tmr.tracker.startTime)/((endFrame-startFrame)*enforcedDeltaT);

            else

                if (nargin > 2)
                    tmr.video.endFrame = endFrame;
                else
                    [~, tmr.video.endFrame] = min(abs(tmr.video.et - tmr.video.et(startFrame) - tmr.tracker.endTime + tmr.tracker.startTime));
                end

                p = polyfit([tmr.video.et(tmr.video.startFrame) tmr.video.et(tmr.video.endFrame)], [tmr.tracker.startTime tmr.tracker.endTime], 1);
                tmr.video.tx = polyval(p, tmr.video.et);
                tmr.video.timewarp = p(1);
            end
            for j = 1:length(tmr.neuron)
                vv = isfinite(tmr.video.tx);
                tmr.video.ratio(j,vv) = resample_to_t(tmr.tx, tmr.neuron(j).ratio, tmr.video.tx(vv));
                tmr.video.ratio_div_baseline(j,vv) = resample_to_t(tmr.tx, tmr.neuron(j).ratio_div_baseline, tmr.video.tx(vv));
            end
            tmr.video.tind = interp1(tmr.tx, 1:length(tmr.tx), tmr.video.tx, 'nearest', 'extrap');
        end
        
        function tmr = trackSpot(tmr)
            imfrac = .25;
            x1 = floor(size(tmr.videocube,2)*(1-imfrac)/2);
            x2 = ceil(size(tmr.videocube,2)*(1+imfrac)/2);
            y1 = floor(size(tmr.videocube,1)*(1-imfrac)/2);
            y2 = ceil(size(tmr.videocube,1)*(1+imfrac)/2);
            vc = tmr.videocube(y1:y2,x1:x2,:);
            thresh = percentile(vc(:), .99);
            vcbw = vc >= thresh;
            for j = 1:size(vcbw,3)
                nopen = 1;
                im = bwmorph(bwmorph(vcbw(:,:,j), 'erode',nopen), 'dilate',nopen);
                stats = regionprops(im, 'area', 'centroid');
                if (isempty(stats))
                    tmr.video.spotLoc(:,j) = [NaN NaN];
                else
            %    imagesc(im);
                   [~,I] = max([stats.Area]);
                    tmr.video.spotLoc(:,j) = [x1 y1] - 1 + stats(I).Centroid;
                end
            end
            
            valid = all(isfinite(tmr.video.spotLoc));
            if (any(~valid))    
                tmr.video.spotLoc(:,~valid) = interp1(find(valid), tmr.video.spotLoc(:,valid)', find(~valid), 'nearest')';
            end        
            
            return;
            
            mi = mean(tmr.videocube, 3)/255;
            myfun = @(x, xdata) min(1,x(1) * exp(-((xdata(:,1)- x(2)).^2 + (xdata(:,2)-x(3)).^2)/(2*x(4))) + x(5));

            op = optimoptions('lsqcurvefit', 'display', 'off');
        
            xi = floor(size(tmr.videocube, 2)*.425):ceil(size(tmr.videocube, 2)*.575);
            yi = floor(size(tmr.videocube, 1)*.425):ceil(size(tmr.videocube, 1)*.575);
            cube = tmr.videocube(yi,xi,:)/255;
            [xx,yy] = meshgrid(xi,yi);
            x0 = [1 mean(xi) mean(yi) 100 0.1];
            lb = [0 min(xi) min(yi) 1 0];
            ub = [100 max(xi) max(yi) size(cube,1).^2 1];
            xdata = [reshape(xx,[],1) reshape(yy,[],1)];
            tmr.video.xfit = repmat(x0, [size(cube,3) 1]);
            
            xlast = x0;
            tic
            for j = 1:size(cube,3)
                if sqrt((xlast(2) - x0(2)).^2 + (xlast(3) - x0(3)).^2) < 10
                    xs = xlast;
                else
                    xs = x0;
                end
                tmr.video.xfit(j,:) = lsqcurvefit(myfun, xs, xdata, reshape(cube(:,:,j), [], 1), lb, ub, op);
                xlast = tmr.video.xfit(j,:);
            end
            tmr.video.spotLoc = tmr.video.xfit(:,2:3)';
        end
        
        %inds = {yinds, xinds, tinds}
        function tmr = trackVideoMotion(tmr, inds)
            if (isempty(tmr.videocube))
                error ('need to load videocube first');
            end
            tmr.video.vflow = zeros([2 size(tmr.videocube,3)]);
            tinds = (tmr.video.startFrame+2):(tmr.video.endFrame-2);
            vc2 = uint8(tmr.videocube(:,:,1:ceil(size(tmr.videocube,3)/1000):size(tmr.videocube,3)));
            for j = 1:size(vc2,3), vc2(:,:,j) = adapthisteq(vc2(:,:,j)); end
            medim = median(vc2,3,'omitnan');
            valid = medim < percentile(medim, 0.8);
            valid = imerode(valid, ones(3));
           
            im = adapthisteq(uint8(tmr.getVidFrame(tinds(1))));
            points = detectMinEigenFeatures(im);
            points = double(points.Location);
            pv = logical(interp2(double(valid), points(:,1), points(:,2), '*nearest'));
            points = points(pv,:);
            pointTracker = vision.PointTracker('MaxBidirectionalError', 1);
            initialize(pointTracker, points, im);         
            npts = length(points);
            for j = tinds(2:end)
                oldim = im;
                % Initialize the tracker with the initial point locations and the initial
                % video frame.
                im = adapthisteq(uint8(tmr.getVidFrame(j)));
                oldPoints = points;
                if (~isempty(points))
                    [points,isFound] = pointTracker.step(im);
                    validPoints = points(isFound,:);
                    validOldPoints = oldPoints(isFound,:);
                else
                    validOldPoints = [];
                    validPoints = [];
                end
                if (length(validPoints) < npts/4)
                    points = detectMinEigenFeatures(oldim);
                    points = double(points.Location);
                    pv = logical(interp2(double(valid), points(:,1), points(:,2), '*nearest'));
                    points = points(pv,:);
                    npts = length(points);
                    oldPoints = points;
                    pointTracker = vision.PointTracker('MaxBidirectionalError', 1);
                    initialize(pointTracker, points, oldim);
                    [points,isFound] = pointTracker.step(im);
                    validPoints = points(isFound,:);
                    validOldPoints = oldPoints(isFound,:);
                end
                
                [xform, ~,points] = estimateGeometricTransform(validOldPoints, validPoints, 'similarity','MaxDistance', 1); 
                pointTracker.setPoints(points);
                 tmr.video.vflow(:,j-1) = xform.T(3,1:2)';
                 if (mod(j,100) == 0)
                     plot (tmr.video.tx(tinds(1):j), cumsum(tmr.video.vflow(1,tinds(1):j)')*tmr.micronsPerPixel, tmr.tx, tmr.tracker.stageloc(1,:)); pause(0.01);
                 end
                
            end
            
            
%             if (~existsAndDefault('inds',{}))
%                 medim = median (tmr.videocube(:,:,tinds), 3);
%                 valid = medim < percentile(medim, 0.98);
%                 goodrow = all(valid,2); 
%                 goodcolumn = all(valid, 1);
%                 if (nnz(goodrow) * size(valid,2) > nnz(goodcolumn) * size(valid,1)) %find direction that yields most good pixels
%                     inds{1} = find(goodrow);
%                     inds{2} = 1:size(valid,2);
%                 else
%                     inds{1} = 1:size(valid,1);
%                     inds{2} = find(goodcolumn);
%                 end
%             end
%             
%             vc = tmr.videocube(inds{1}, inds{2}, tinds);
%             for j = 2:(size(vc,3) - 1)
%                  output = dftregistration(fft2(vc(:,:,j-1)), fft2(vc(:,:,j+1)), 20); %[output, greg] for im registration
%                  tmr.video.vflow(:,tinds(1)+j-1) = -0.5*output([4 3]);               
%             end
            tmr.video.pixelDisplacement = cumsum(tmr.video.vflow, 2);
        end
        
        function tmr = getVideoToStageMag(tmr)
        %    frameTime = median(diff(tmr.video.tx));
            f = 0.5;
            vinds = floor((0.5 + f/2)*tmr.video.startFrame + (0.5-f/2)*tmr.video.endFrame):ceil((0.5-f/2)*tmr.video.startFrame + (0.5+f/2)*tmr.video.endFrame);
    
            sl = interp1(tmr.tx, tmr.tracker.stageloc(1:2,:)', tmr.video.tx(vinds))';
            vd = cumsum(tmr.video.vflow(:,vinds),2);
            vd(3,:) = 1;
            at = sl/vd;
            vdf = at*vd;
            plot (sl(1,:), sl(2,:), vdf(1,:), vdf(2,:));
            
            tmr.vidToStageTransform = at(1:2,1:2);
            tmr.micronsPerPixel = sqrt(det(tmr.vidToStageTransform));
            sl(3,:) = 1;
            vd = vd(1:2,:);
            at = vd/sl;
            tmr.stageToVidTransform = at(1:2,1:2);
            tmr.video.micronDisplacement = tmr.vidToStageTransform * tmr.video.pixelDisplacement;
            %dsl = 0.5*(sl(3:end,:) - sl(1:(end-2),:));
                    %dsl = dsl((vinds-1),:);
            
        end
        
        function tmr = findVideoToStageTimeWarp(tmr)
            existsAndDefault('scale', 1);
            %inds = (ceil(tmr.video.startFrame+scale/2):scale:floor(tmr.video.endFrame-scale/2))';
            inds = (tmr.video.startFrame+1):(tmr.video.endFrame - 1);
            dtvid = (tmr.tracker.endTime - tmr.tracker.startTime)/(tmr.video.endFrame - tmr.video.startFrame);   
            vvid = tmr.micronsPerPixel*tmr.video.vflow(:,inds) / dtvid;
            sl = tmr.tracker.stageloc(1:2,:);
            sl(:,any(~isfinite(sl))) = interp1(find(all(isfinite(sl))), sl(:,all(isfinite(sl)))', find(any(~isfinite(sl))), 'nearest', 'extrap')';
            vstage = deriv(sl, 2/3 * dtvid/tmr.dt)/tmr.dt;
            rfac = floor(dtvid/(tmr.dt));
            vsr = lowpass1D(vstage', floor(rfac/2));
            vsr = downsample(vsr, rfac, floor(rfac/2));
            txr = downsample(tmr.tx, rfac, floor(rfac/2));
%             inds = find(txr < tmr.tracker.startTime | txr <= txr(1), 1, 'last'):find(txr > tmr.tracker.endTime | txr >= txr(end), 1, 'first');
%             txr = txr(inds);
%             vsr = vsr(inds, :);
            sspeed = sqrt(sum(vsr.^2,2));
            percentile(sspeed, 0.9)
            [tfit,mag] = energyMinTAlignMultiScale(txr, vsr, vvid, [tmr.tracker.startTime tmr.tracker.endTime], dtvid*0.8, 20*percentile(sspeed, 0.9)^2);
            tmr.video.tx(inds) = tfit;
            tmr.video.micronsPerPixel = mag*tmr.micronsPerPixel;
        end
            
        
        %{
        function [tfit,inds] = findVideoToStageTimeWarp(tmr, scale)
            existsAndDefault('scale', 1);
            %inds = (ceil(tmr.video.startFrame+scale/2):scale:floor(tmr.video.endFrame-scale/2))';
            inds = (tmr.video.startFrame+1):(tmr.video.endFrame - 1);
            dtvid = (tmr.tracker.endTime - tmr.tracker.startTime)/(tmr.video.endFrame - tmr.video.startFrame);   
            vvid = tmr.micronsPerPixel*tmr.video.vflow / dtvid;
            vstage = deriv(tmr.tracker.stageloc(1:2,:), 2/3 * dtvid/tmr.dt)/tmr.dt;
            
            if (length(inds) > 100)
                [tf,ii] = tmr.findVideoToStageTimeWarp(scale*2);
                ti = interp1([tmr.video.startFrame; ii; tmr.video.endFrame], [tmr.tracker.startTime; tf; tmr.tracker.endTime], inds, 'linear');
            else
                ti = [];
            end
            if (~isempty(ti))
                 tends = [tmr.tracker.startTime, tmr.tracker.endTime];
                dtvid =  diff(tends)/(length(inds) + 1);
            
                dtstage = dtvid/3;
                tx1 = max(tends(1) - dtvid * 6, tmr.tx(1));
                tx2 = min(tends(2) + dtvid * 6, tmr.tx(end));
                txs = tx1:dtstage:tx2;
                sl = interp1(tmr.tx, tmr.tracker.stageloc(1:2,:)', txs)';
                vs = deriv(sl, 2)/dtstage;
            
                
               % vu = vu(:,inds);
                figure();
                plot (txs, vs(1,:), 'b.-', ti, vu(1,inds),'k.-', tf, vu(1,ii), 'r.-'); title(['scale = ' num2str(scale)]);
                 xlim([2440 2460]);
            
           
        %    plot (txs, vs(1,:), tinitial, vu(1,:)); pause
     %           figure(3);
    %            plot (ii,tf(2:end-1),inds,ti); pause
            end
            tfit = findVideoToStageTimeWarpSingleScale(tmr, inds, ti);
            title(['after fit at scale = ' num2str(scale)])
        end
        
        function tfit = findVideoToStageTimeWarpSingleScale (tmr, inds, tintitial, vidvel)
            tends = [tmr.tracker.startTime, tmr.tracker.endTime];
            dtvid =  diff(tends)/(length(inds) + 1);         
            
             dtstage = dtvid/3;
             tx1 = max(tends(1) - dtvid * 6, tmr.tx(1));
             tx2 = min(tends(2) + dtvid * 6, tmr.tx(end));
             txs = tx1:dtstage:tx2;
             sl = interp1(tmr.tx, tmr.tracker.stageloc(1:2,:)', txs)';
             vs = deriv(sl, 2)/dtstage;
            
            vidv = lowpass1D(vidvel, median(diff(inds))/2);
            vu = vidv(:,inds);
            existsAndDefault('tinitial', linspace(tends(1) + dtvid/2, tends(2)-dtvid/2, length(inds)));
        %    plot (txs, vs(1,:), tinitial, vu(1,:)); pause
            tfit = energyMinTAlign(txs, vs, vu, tends, dtvid/3, tinitial);
            
            figure();
            plot (txs, vs(1,:), 'b.-', tinitial, vu(1,:), 'k.-', tfit, vu(1,:), 'r.-'); 
             xlim([2440 2460]);
            
            
        end
        %}
        function tmr = findVideoToStageMagViterbi(tmr, vinds, fps)
            existsAndDefault('fps', round((tmr.video.endFrame - tmr.video.startFrame)/(tmr.tracker.endTime - tmr.tracker.startTime)));
            expectedRatio = fps/((tmr.video.endFrame - tmr.video.startFrame)/(tmr.tracker.endTime - tmr.tracker.startTime));
            alpha = min(expectedRatio, 1) - 0.1;
            beta = 1-alpha;
            existsAndDefault('vinds', floor(alpha*tmr.video.startFrame + beta*tmr.video.endFrame):ceil(beta*tmr.video.startFrame + alpha*tmr.video.endFrame));
             frametime = 1/fps;
            maxforward = 3; %allow it to skip 1 or two frames
            minback = 1; %make it go forward at least 1 frame every time
            sl = interp1(tmr.tx, tmr.tracker.stageloc(1:2,:)', (tmr.tracker.startTime-0.5*frametime):frametime:(tmr.tracker.endTime+0.5*frametime));
            dsl = diff(sl(2:end,1:2));
            vf = tmr.video.vflow(1:2,vinds)'; 
            
            %first alignment is to direction of movement only -- no need to
            %determine magnification
            svn = dsl./repmat(sqrt(sum(dsl.^2,2)) + eps, [1 2]);
            vfn = vf./repmat(sqrt(sum(vf.^2,2)) + eps, [1 2]);
            distfun = @(x,y) sum((x-y).^2,2);
            
            path = viterbi(svn, vfn, distfun, minback, maxforward);
            
            vidtx = interp1(tmr.tracker.startTime:frametime:tmr.tracker.endTime, path);
            
            
            %find magnification by matching displacement
            %assume that the middle 50% is most stable

            vfd = cumsum(vf,1);
            slfit = interp1(tmr.tx, tmr.tracker.stageloc', vidtx);
            
            
            inds = floor(length(vidtx)/4):ceil(3*length(vidtx)/4);
            vfd=vfd(inds,:); vfd = vfd - repmat(mean(vfd,1), [size(vfd,1) 1]);
            slfit=slfit(inds,1:2); slfit = slfit - repmat(mean(slfit,1), [size(slfit,1) 1]);
            
            tmr.vidToStageTransform = vfd\slfit;
            tmr.stageToVidTransform = slfit\vfd;
            vfr = vfd*tmr.vidToStageTransform;
            figure(1), plot (slfit(:,1), slfit(:,2), vfr(:,1), vfr(:,2)); 
            
            %now redo alignment with new estimate of magnification and
            %actual displacements
            oldpath = path;
            path = viterbi(dsl, vf*tmr.vidToStageTransform, distfun, minback, maxforward);
            
            vidtx = interp1(tmr.tracker.startTime:frametime:tmr.tracker.endTime, path);
           % tmr.video.tx(vinds) = vidtx;
            
            vfd = cumsum(vf,1);
            slfit = interp1(tmr.tx, tmr.tracker.stageloc', vidtx);
            
            
            inds = floor(length(vidtx)/4):ceil(3*length(vidtx)/4);
            
            vfd=vfd(inds,:);vfd = vfd - repmat(mean(vfd,1), [size(vfd,1) 1]);
            slfit=slfit(inds,1:2); slfit = slfit - repmat(mean(slfit,1), [size(slfit,1) 1]);
            tmr.vidToStageTransform = vfd\slfit;
            tmr.stageToVidTransform = slfit\vfd;
            figure(2), plot (vinds, tmr.video.tx(vinds), vinds, vidtx); 
            vidtx = interp1(tmr.tracker.startTime:frametime:tmr.tracker.endTime, path);
            tmr.video.tx(vinds) = vidtx;
            vfr = vfd*tmr.vidToStageTransform;
            figure(3),plot (slfit(:,1), slfit(:,2), vfr(:,1), vfr(:,2)); 
            
        end
        
        
        function tmr = trackVideoMotionOld(tmr)
            tmr.video.vflow = zeros(2, size(tmr.videocube,3));
            minradsq = min(size(tmr.videocube(:,:,1))).^2/4; %            size(tmr.videocube,1)*size(tmr.videocube,2) / (2 * pi);
            xc = size(tmr.videocube,2)/2;
            yc = size(tmr.videocube,1)/2;
            if (isempty(tmr.videogradients))
                [xd,yd,it] = imgradient(tmr.videocube,1);
                tmr.videogradients.xd = xd;
                tmr.videogradients.yd = yd;
                tmr.videogradients.it = it;
            else
                xd = tmr.videogradients.xd;
                yd = tmr.videogradients.yd;
                it = tmr.videogradients.it;
            end
            [xx,yy] = meshgrid(1:size(tmr.videocube,1), 1:size(tmr.videocube,2));
            ll = [xx(:) yy(:)];
            valid = (ll(:,1)-xc).^2 + (ll(:,2) - yc).^2 > minradsq & ll(:,1) > 1 & ll(:,2) > 1 & ll(:,1) <= size(tmr.videocube,2) & ll(:,2) <= size(tmr.videocube,1); 
            ll = ll(valid,:);
            for j = 1:size(tmr.videocube,3)
%                 pts = detectMinEigenFeatures(tmr.videocube(:,:,j));
%                 ll = double([pts.Location]);
%                 dx = ones(size(ll)); dx(:,2) = 0;
%                 dy = ones(size(ll)); dy(:,1) = 0;
%                 ll = [ll;ll+dx;ll+dy;ll-dx;ll-dy;ll-dx-dy;ll-dx+dy;ll+dx-dy;ll+dx+dy]; %#ok<AGROW>
                
                tmr.video.vflow(:,j) = TrackingMicroscopeResult.dirtyLK(xd(:,:,j), yd(:,:,j), it(:,:,j), ll);
%                 if (mod(j,100) == 0)
%                     x = cumsum(tmr.video.vflow(1,1:j)); y = cumsum(tmr.video.vflow(2,1:j));
%                     plot(x,y); axis equal; title (num2str(j)); pause(0.01);
%                 end
            end
        end
        function tmr = alignVideoToSpot(tmr)
            if (~isfield(tmr.video, 'spotLoc'))
                tmr = tmr.trackSpot();
            end
            nloc = interp1(tmr.tx, mean(tmr.tracker.loc(1:2,:,:),3)', tmr.video.tx')';
            nl = nloc; nl(3,:) = 1;
            sloc = tmr.video.spotLoc;
            inds = all(isfinite(nl)) & all(isfinite(sloc));
            
            at_las_v = tmr.video.spotLoc(:,inds)*pinv(nl(:,inds));
            sloc2 = at_las_v*nl;
            dl = sqrt(sum((sloc - sloc2).^2));
            valid = dl <= 3*mean(dl(inds)) & inds;
            
            at_las_v = tmr.video.spotLoc(:,valid)*pinv(nl(:,valid));
            
            tmr.video.at_las_v = at_las_v;
            tmr.video.spotLocLaser = interp1(tmr.video.tx(inds), (at_las_v * nl(:,inds))', tmr.video.tx')';
            for k = 1:size(tmr.tracker.loc,3)
                nloc = interp1(tmr.tx, tmr.tracker.loc(1:2,:,k)', tmr.video.tx')';
                nl = nloc; nl(3,:) = 1;
                tmr.video.nLocLaser(:,:,k) = interp1(tmr.video.tx(inds), (at_las_v * nl(:,inds))', tmr.video.tx')';
            end
            tmr.video.micronsPerPixel = 1/sqrt(abs(det(at_las_v(1:2,1:2))));
        end
        
        function tmr = alignVideoToStage(tmr)
            if (~isfield(tmr.video, 'vflow'))
                tmr = tmr.trackVideoMotion();
            end
            sloc = interp1(tmr.tx, tmr.tracker.stageloc', tmr.video.tx)';
            dstage100 = sloc(1:2,101:end)- sloc(1:2,1:(end-100)); inds = all(isfinite(dstage100));
            xs = cumsum(tmr.video.vflow,2);
            dxs100 = xs(:,101:end)-xs(:,1:(end-100));
            
            s = dstage100(:,inds);
            v = dxs100(:,inds);
            
            ds = interp1((tmr.tx(1:end-1) + tmr.tx(2:end))*.5, diff(tmr.tracker.stageloc',1), tmr.video.tx')';
            inds = all(isfinite(ds));
            s = lowpass1D(ds(1:2,inds), 10);
            v = lowpass1D(tmr.video.vflow(:,inds), 10);
            
            %find points that are roughly a rotation and scaling apart
            sc = s - repmat(mean(s,2),[1 size(s,2)]);
            vc = v - repmat(mean(v,2),[1 size(v,2)]);
            sn = mean(sqrt(sum(sc.^2)));
            vn = mean(sqrt(sum(vc.^2)));
            m = sn/vn; %mag: pixels per micron
            [U,~,V] = svd( (sc/sn)*(vc/vn).');
            R = V*U.'; R = [1 0; 0 det(R)]*R;
            vnfit = 1/m * R * sc;
            dx = vc - vnfit; valid = sqrt(sum(dx.^2)) < 1.5*mean(sqrt(sum(dx.^2)));

            
             tmr.video.at_v_s = s(1:2,valid)*pinv(v(:,valid)); %vid to stage
             tmr.video.at_s_v = v(:,valid)*pinv(s(1:2,valid)); %stage to vid
            tmr.video.micronsPerPixel = sqrt(det(tmr.video.at_v_s));
            tmr.video.stageLoc = sloc(1:2,:);
            tmr.video.stageLocPix = tmr.video.at_s_v * sloc(1:2,:);
            tmr.video.vidLocPix = cumsum(tmr.video.vflow, 2);
        end
        
        %gets options that result in movie playback 
        %when larva is moving, zoomed in, and 
        function ops = moviePlaybackOptions_default (tmr, varargin)
%             lptime = 2;
%             spthresh = 10; %um/sec
            rectsize = 1.5; %mm
            scalebarlength = 0.5; %mm
            varargin = assignApplicable(varargin);
%             lpspeed = lowpass1D(tmr.neuron(1).speed_filt, lptime/tmr.dt);
%             if (percentile(lpspeed, 0.5) < spthresh)
%                 warning ('looks like this does not move very fast -- check for problems');
%             end
%             spthresh = min(spthresh, percentile(lpspeed, 0.5)); %avoid something weird if it never moves
%             
%             lpi = floor((3*lptime)/tmr.dt):ceil(length(lpspeed) - (3*lptime)/tmr.dt);
%             startind = find(lpspeed(lpi) > spthresh, 1, 'first') + lpi(1);
%             endind = find(lpspeed(lpi) > spthresh, 1, 'last') + lpi(1);
%             timerange = tmr.tx([startind endind]);
            sl = tmr.tracker.stageloc; sl1 = sl(:,find(all(isfinite(sl),1), 1, 'first'));
            slend = sl(:,find(all(isfinite(sl),1), 1, 'last'));

            sd = sqrt(sum((tmr.tracker.stageloc - repmat(sl1, [1 size(tmr.tracker.stageloc,2)])).^2));
            ind1 = find(sd > 50, 1, 'first'); %when stage moves 50 microns from start
            sd = sqrt(sum((tmr.tracker.stageloc - repmat(slend, [1 size(tmr.tracker.stageloc,2)])).^2));
            ind2 = find(sd > 50, 1, 'last'); %when stage moves 50 microns from start
            timerange = tmr.tx([ind1 ind2]);

            vcenter = median(tmr.video.spotLoc, 2, 'omitnan');
            xl = vcenter(1) + [-500 500]*rectsize/tmr.micronsPerPixel;
            yl = vcenter(2) + [-500 500]*rectsize/tmr.micronsPerPixel;
            lpfilt = true;
            ratrange = percentile(tmr.lpfilt(tmr.neuron(1).ratio_div_baseline(ind1:ind2)), [0.01 0.99]);
           
            ops = {'timerange', timerange, 'xl', xl, 'yl', yl, 'rectsize', [], 'overlaytime', true, 'scalebarlength', scalebarlength, 'ratrange', ratrange, 'lpfilt', lpfilt};
            
        end
        
        function playMovie (tmr, inds, vw, varargin)
            scalebarlength = 1;
            timerange = [];
            stagelpfreq = 10;
            ratrange = [0.5 3.5];
            flipx = false;
            lpfilt = true;
            spotcolors = {[1 1 1], [0 1 1], [1 0 1], [1 1 0]};
            spotnames = {};
            overlayFluorescence = true;
            varargin = assignApplicable(varargin);
            
            
            if (~tmr.hasVideo)
                disp ('no movie');
                return;
            end
            if (~isfield(tmr.video, 'tind'))
                tmr.video.tind = interp1(tmr.tx, 1:length(tmr.tx), tmr.video.tx, 'nearest', 'extrap');
            end
            if (~isfield(tmr.video, 'micronDisplacement'))
                tmr.video.micronDisplacement = interp1(tmr.tx, lowpass1D(tmr.tracker.stageloc(1:2,:), 1000/(stagelpfreq*7.3))', tmr.video.tx)';
            end
            
           
            
            if (~existsAndDefault('inds', tmr.video.startFrame:tmr.video.endFrame))
                if (~isempty(timerange))
                    inds = tmr.video.tx > min(timerange) & tmr.video.tx < max(timerange);
                end
            end
                 
            if (islogical(inds))
                inds = find(inds);
            end
            closevw = false;
            if (existsAndDefault('vw',[]))
                if (ischar(vw))
                    if (strcmpi(vw, 'record'))
                        ddd = fileparts(tmr.video.filename);
                        vw = fullfile(ddd, 'tracking-movie.mp4');
                        jjj = 0;
                        while (exist(vw, 'file'))
                            jjj = jjj+1;
                            vw = fullfile(ddd, ['tracking-movie-', num2str(jjj), '.mp4']);
                        end
                    end
                    vw = VideoWriter(vw, 'mpeg-4');
                    vw.Quality = 95;
                    vw.open();
                    closevw = true;
                end
            end
            fluorOnTop = false;%diff(tmr.video.tx(inds([1 end]))) < 15 && length(tmr.neuron) > 1;
            
            try        
            
            clf;
            xpix = size(tmr.videocube,2);
            ypix = size(tmr.videocube,1);
            
            ax_movie = axes();
            ax_movie.Units = 'pixels';
            pp = ax_movie.Position;
            dw = xpix - pp(3);
            dh = ypix - pp(4);
            ax_movie.Position = pp + [-dw/2 -dh/2 dw dh];
            ax_movie.Units = 'normalized';
           
            ax_movie = tmr.drawFrame(inds(1), ax_movie, varargin{:}, 'ratrange', ratrange);
            ax_movie.Units = 'pixels';
            avirect = ax_movie.Position;
            movrect = avirect;
            movrect(4) = movrect(4)*1.25;
            ax_movie.Units = 'normalized';
           
            pp = ax_movie.Position;
            toploc = [pp(1) pp(2)+pp(4) pp(3) pp(4)/4];
            
            xl = get(ax_movie, 'xlim');
            yl = get(ax_movie, 'ylim');
             if (length(tmr.neuron) == 1)
                npl = 3;
            else
                npl = length(tmr.neuron);
             end
             if (overlayFluorescence)
                if (fluorOnTop)
                    ll = toploc;
                    avirect = movrect;
                    ax_pl = axes('Position', ll, 'Color', 'none');                %#ok<AGROW>
                else
                    ll = dsxy2figxy_marc(ax_movie, [xl(1)*.95+xl(2)*.05 yl(1)*.25+yl(2)*.75 diff(xl)/4 diff(yl)/4]);
                    h = ll(4)/npl;
                    for j = 1:npl
                        ax_pl(j) = axes('Position', [ll(1) ll(2)+(npl-j)*h ll(3) h], 'Color', 'none');                %#ok<AGROW>
                    end
                end
             
           
                tspread = [-10 0];

                if (fluorOnTop)
                    trange = tmr.video.tx([min(inds) max(inds)]);
                    tmr.plotFluorescence(ax_pl, trange, ratrange, lpfilt, spotcolors);
                    set(ax_pl, 'XLim', trange, 'Color', 'w');
                else
                    trange = tmr.video.tx([min(inds) max(inds)]) + tspread;
                    tmr.plotFluorescence(ax_pl, trange, ratrange, lpfilt,spotcolors);
                    set(ax_pl, 'XLim', tmr.video.tx(inds(1)) +  tspread, 'Color', 'k');
                end
          
            
            ctrloc = dsxy2figxy_marc(ax_pl(1), [tmr.video.tx(inds(1)) mean(get(ax_pl(1), 'YLim'))],[tmr.video.tx(inds(1)) mean(get(ax_pl(1), 'YLim'))]);
            x = ctrloc([1 1]);
            la = annotation('line', x, [ll(2) ll(2)+ll(4)], 'Color', 'w', 'LineStyle', ':', 'LineWidth', 2);
%             
           
            
            if (fluorOnTop)
                numsec = 0.5;
                yll = [0.4 max(max(tmr.video.ratio_div_baseline(:,inds)))*1.1];
                ylim(ax_pl(1),yll);
                y0 = 1.75;%max(get(ax_pl(1), 'ylim'))- numsec*1.25;
                numx = ceil((yll(2)-y0)/2)/2;
                
                x0 = tmr.video.tx(inds(1)) + .25;
                
                [tb,tt] = labeledBar(ax_pl(1), x0 + [0 numsec], [y0 y0], [num2str(numsec) 's'],  'below', 'k', {'LineWidth', 2}, {'FontSize', 8});
                [yb,yt] = labeledBar(ax_pl(1), [x0 x0], y0 + [0 numx], [num2str(numx) 'x'],  'right', 'k', {'LineWidth', 2}, {'FontSize', 8});
                %yt.HorizontalAlignment = 'left'; yt.Margin = 0;
            else
                numsec = 2;
                numx = 1;
                 x0 = tmr.video.tx(inds(1)) + tspread(end)+0.05*diff(tspread);
                 y0 = 1;
                 if (isempty(spotnames))
                    [tb,tt] = labeledBar(ax_pl(1), x0 + [0 numsec], [y0 y0], [num2str(numsec) 's'],  'below', 'w', {'LineWidth', 2}, {'FontSize', 8});
                    [yb,yt] = labeledBar(ax_pl(1), [x0 x0], y0 + [0 numx], [num2str(numx) 'x'],  'above', 'w', {'LineWidth', 2}, {'FontSize', 8});
                    yt.HorizontalAlignment = 'left'; yt.Margin = 0;
                 else
                     for j = 1:length(spotnames)
                         labeledBar(ax_pl(j), [x0 x0], [y0 y0], spotnames{j},  'right', spotcolors{j},'noline', {'FontSize', 10, 'FontWeight', 'bold'});
                     end
                     x0 = x0-numsec;
                     y0 = ax_pl(end).YLim(1) - numx;
                     [tb,tt] = labeledBar(ax_pl(end), x0 + [0 numsec], [y0 y0], [num2str(numsec) 's'],  'below', 'w', {'LineWidth', 2}, {'FontSize', 8});
                    [yb,yt] = labeledBar(ax_pl(end), [x0 x0], y0 + [0 numx], [num2str(numx) 'x'],  'right', 'w', {'LineWidth', 2}, {'FontSize', 8});
                    yt.HorizontalAlignment = 'right'; %yt.Margin = 0;
                 end
            end
             end
            xl = get(ax_movie, 'XLim');
            yl = get(ax_movie, 'YLim');
            
            [sb,st] = labeledBar(ax_movie, xl(2)*.95 + [-1000*scalebarlength/tmr.micronsPerPixel 0], yl(2)*.95 + 0.05*yl(1) + [0 0], [num2str(scalebarlength) ' mm'], 'below', 'w', {'LineWidth', 5}, {'Margin', 5});
            if (tmr.tracker.hasStatus && any(tmr.tracker.status.stimOn))
                [stimbar, stimtext] = labeledBar(ax_movie, mean(xl) + diff(xl) * [-.1 .1], yl(1)*.95 + 0.05*yl(2) + [0 0], 'stim on', 'above', 'c', {'LineWidth', 5}, {'Margin', 5});            
                sstate = logical(interp1(tmr.tx, double(tmr.tracker.status.stimOn), tmr.video.tx, 'nearest', 'extrap'));
                sinfo = true;
            else
                sinfo = false;
            end
            
            
            for j = 1:length(inds)
                tmr.drawFrame(inds(j), ax_movie, varargin{:},'ratrange', ratrange, 'spotcolors', spotcolors);
                if (flipx)
                    set(ax_movie, 'XDir', 'reverse');
                end
                if (overlayFluorescence)
                    if (fluorOnTop)
                        c = {[1 1 1], [0 1 1], [1 0 1], [1 0 0]};

                        for k = 1:length(tmr.neuron);    
                            hold(ax_pl, 'on');
                            hht(k) = plot (ax_pl, tmr.video.tx(inds(j)), tmr.video.ratio_div_baseline(k,inds(j)), 'ko', 'Color', c{mod(k,4)+1}, 'LineWidth', 2); 
                        end
                        hht(end+1) = plot (ax_pl, tmr.video.tx(inds(j)) + [0 0], get(ax_pl,'YLim'), 'k:', 'Color', [.6 .6 .6], 'LineWidth', 2); 


                    else
                        set(ax_pl, 'XLim', tmr.video.tx(inds(j)) +  tspread, 'Color', 'k');
                    end
                end
                if (sinfo)
                    if (sstate(inds(j)))
                        stimbar.Visible = 'on';
                        stimtext.Visible = 'on';
                    else
                        stimbar.Visible = 'off';
                        stimtext.Visible = 'off';
                    end
                end
                
                writeFrameOrPause(vw, avirect, 20);
                if (fluorOnTop && overlayFluorescence)
                    delete(hht);
                    clear hht;
                end
            end
            catch me
                disp(me.getReport());
                if (closevw)
                    vw.close();
                end
                return;
            end
            if (closevw)
                vw.close();
            end
            
        end
        
        function [tmr,path,sl,xform] = manualRegistration (tmr, inds, rad,tr)
           
            if (~existsAndDefault('inds', []))
                existsAndDefault('tr', [tmr.tracker.startTime, tmr.tracker.endTime]);
                segmentLength = 100;
                nsegments = 5;
                stageloc = interp1(tmr.tx, tmr.tracker.stageloc(1:2,:)', tmr.video.tx)';
                dstageloc = sum((stageloc(:,(1+segmentLength):end) - stageloc(:,1:(end-segmentLength))).^2);
                segtime = tmr.video.tx((1+segmentLength):end);
                [~,I] = sort(dstageloc, 'descend');
                I = I(isfinite(dstageloc(I)) & segtime > min(tr) & segtime < max(tr));
                si(1) = I(1);
                for j = 2:(nsegments-2)
                    valid = true(size(I));
                    for k = 1:(j-1)
                        valid = valid & abs(I - si(k)) > segmentLength;
                    end
                    si(j) = I(find(valid, 1, 'first'));
                end
                dstageloc = abs((stageloc(1,(1+segmentLength):end) - stageloc(1,1:(end-segmentLength)))); %find biggest x change not looked at yet
                [~,I] = sort(dstageloc, 'descend');
                I = I(isfinite(dstageloc(I)));
                j = nsegments-1;
                valid = true(size(I));
                for k = 1:(j-1)
                    valid = valid & abs(I - si(k)) > segmentLength;
                end
                si(j) = I(find(valid, 1, 'first'));
                    
                dstageloc = abs(stageloc(2,(1+segmentLength):end) - stageloc(2,1:(end-segmentLength))); %find biggest y change not looked at yet
                [~,I] = sort(dstageloc, 'descend');
                I = I(isfinite(dstageloc(I)));
                j = nsegments;
                valid = true(size(I));
                for k = 1:(j-1)
                    valid = valid & abs(I - si(k)) > segmentLength;
                end
                si(j) = I(find(valid, 1, 'first'));
                
                
                for j = 1:nsegments
                    inds{j} = si(j):(si(j) + segmentLength);
                end
            end
                    
            
            existsAndDefault('rad', 20);
            if ~(iscell(inds))
                inds = {inds};
            end
            for k = 1:length(inds)
                path{k} = tmr.trackPoints(inds{k}, inf, rad, 'high contrast fixed point on slide');
                valid = all(isfinite(path{k}));
                path{k} = path{k}(:,valid);
                inds{k} = inds{k}(valid);
                sl{k} = interp1(tmr.tx, tmr.tracker.stageloc(1:2,:)', tmr.video.tx(inds{k}))';
                useme(k) = ~isempty(path(k));
            end
            path = path(useme);
            sl = sl(useme);
            for j = 1:length(path)
                pp{j} = path{j}-repmat(mean(path{j},2), [1 size(path{j},2)]);
                sll{j} = sl{j}-repmat(mean(sl{j},2), [1 size(sl{j},2)]);
            end
            ppp = cat(2, pp{:})';
            slll = cat(2, sll{:})';
            xform = estimateGeometricTransform(slll, ppp, 'similarity');
            tmr.stageToVidTransform = xform.T(1:2,1:2);
            
            tmr.vidToStageTransform = inv(tmr.stageToVidTransform);
            tmr.micronsPerPixel = sqrt(det(tmr.vidToStageTransform));
%             stageToPixel = estimateRotScalingTransformation(dsl', dp', 458, 1);
%             tmr.stageToVidTransform = stageToPixel;
%             tmr.vidToStageTransform = inv(stageToPixel);
            
        
        end
        
        function tmr = addTailSpotVec(tmr)
            if (~isfield(tmr.video, 'spotLocLaser'))
                disp ('need spotLocLaser')
                return;
            end
            if (~isfield(tmr.video, 'tailloc'))
                disp ('need tailloc');
                return;
            end
            if (isempty(tmr.vidToStageTransform))
                disp ('need vidToStageTransform');
                return;
            end
            valid = all(isfinite(tmr.video.tailloc)) & all(isfinite(tmr.video.spotLocLaser));
            tsvec = -tmr.vidToStageTransform'*(tmr.video.spotLocLaser(:,valid) - tmr.video.tailloc(:,valid)); %not sure why the - sign is here but will think about it
            tmr.posture.tailSpotVec = interp1(tmr.video.tx(valid), tsvec', tmr.tx, 'linear', NaN)';
            tmr.posture.tailSpotDist = sqrt(sum(tmr.posture.tailSpotVec.^2));
            tmr.posture.tailSpotDir = tmr.posture.tailSpotVec./(repmat(tmr.posture.tailSpotDist, [2 1]));
             for j = 1:length(tmr.neuron)
                 tmr.neuron(j).vel_rel_ts = dot(tmr.posture.tailSpotDir,tmr.neuron(j).vel(1:2,:)); 
             end
            
        end
        
        function [tmr,a, ptlist] = manualRegistrationOld (tmr, inds)
            
            if (~tmr.hasVideo)
                disp ('no movie');
                return;
            end
            if (~isfield(tmr.video, 'tind'))
                tmr.video.tind = interp1(tmr.tx, 1:length(tmr.tx), tmr.video.tx, 'nearest', 'extrap');
            end
            
            
            try        
            
                clf;
                xpix = size(tmr.videocube,2);
                ypix = size(tmr.videocube,1);
                
                ax_movie = axes();
                ax_movie.Units = 'pixels';
                pp = ax_movie.Position;
                dw = xpix - pp(3);
                dh = ypix - pp(4);
                ax_movie.Position = pp + [-dw/2 -dh/2 dw dh];
                ax_movie.Units = 'normalized';
                
                ax_movie = tmr.drawFrame(inds(1), ax_movie, 'breadcrumbs', false);
                [x,y] = getpts();
                ptlist(:,1) = [x(end);y(end)];
                for j = 2:length(inds)
                    tmr.drawFrame(inds(j), ax_movie, 'breadcrumbs', false);
                    [x,y] = getpts();
                    ptlist(:,j) = [x(end);y(end)];
                end
                stageloc = tmr.tracker.stageloc(1:2,tmr.video.tind(inds));
                dst = stageloc - repmat(stageloc(:,1),[1 size(stageloc,2)]);
                dpt = ptlist - repmat(ptlist(:,1),[1 size(ptlist,2)]);
                a = dst*pinv(dpt)
                micronsPerPixel = sqrt(det(a))
                [v,d] = eig(a)
                tmr.micronsPerPixel = micronsPerPixel;
            catch me
                disp(me.getReport());
            end
            
            
        end
        
        function ax = drawFrame (tmr, j, ax, varargin)
            if (nargin < 3)
                ax = gca();
            end
            
            im = tmr.getVidFrame(j);
            
            rectsize = min(3000/tmr.micronsPerPixel, min(size(im,1),size(im,2)));
            xl = [];
            yl = [];
            flipy = false;
            breadcrumbs = true;
            crumb0 = [];
            edgesigma = -1;
            ratrange = [0.2 4.2];
            overlaytime = false;
            overlayt0 = 0;
            overlayspot = true;
            overlaybehavior = true;
            spotcolors = {[1 1 1], [0 1 1], [1 0 1], [1 1 0]};
            marktail = false;
            varargin = assignApplicable(varargin);
            
            if ~isempty(rectsize)
                if (length(rectsize) == 1)
                    rectsize = rectsize([1 1 ]);
                end
                sc = median(tmr.video.spotLocLaser(:,all(isfinite(tmr.video.spotLocLaser))),2);
                xl = sc(1) + rectsize(1)*[-.5 .5];
                yl = sc(2) + rectsize(2)*[-.5 .5];
            end

            if (edgesigma > 0)
                [xd,yd] = imgradient(im, edgesigma);
                pcolor(ax, sqrt(xd.^2 + yd.^2));
            else
                pcolor(ax, adapthisteq(uint8(im)));
            end
            shading (ax, 'flat'); 
            if (flipy)
                axis(ax, 'ij');
            end
            if (~isempty(xl))
                xlim(ax, xl);
            end
            if (~isempty(yl))
                ylim(ax, yl);
            end
            axis (ax, 'equal'); 
            if (~isempty(xl))
                xlim(ax, xl);
            end
            if (~isempty(yl))
                ylim(ax, yl);
            end
            if (isempty(xl) && isempty(yl))
                axis (ax, 'tight'); 
            end
            
            axis (ax, 'off');
            colormap (ax, 'gray'); 
            hold (ax, 'on');
            
            cmap = hot(256);
            
            nspots = length(tmr.neuron);
            
            basespotr = 5;
            if (nspots > 1)
                offset = squeeze(tmr.video.nLocLaser(1:2,j,:)); 
                offset = offset - repmat(mean(offset,2), [1 size(offset,2)]);
                offset = 1.25*basespotr * offset./repmat(sqrt(sum(offset.^2,1)), [size(offset,1) 1]);
            else
                offset = [0; 0];
            end
            
            %rr = percentile(tmr.video.ratio(1:100:end), [0.05 .5 .95]);
             clist = spotcolors;
             if (overlayspot)
                for k = 1:nspots
                   % mr = median(tmr.video.ratio(k,:));
                    rr = ratrange;
                    spotr = basespotr*sqrt(min(max(rr(1),tmr.video.ratio_div_baseline(k,j)),rr(2)));

                    c = round(256*(tmr.video.ratio_div_baseline(k,j)-rr(1))/diff(rr)); c = max(c,1); c = min(c, 256);
                    fc = cmap(c,:);
                    if (nspots == 1)
                        ec = fc;
                    else
                        ec = clist{mod(k-1,length(clist))+1};
                    end
                    npp = 24;
                    patch(tmr.video.spotLocLaser(1,j) + offset(1,k) + spotr*cos(linspace(0,2*pi,npp)), tmr.video.spotLocLaser(2,j) + offset(2,k) + spotr*sin(linspace(0,2*pi,npp)), fc, 'Parent', ax, 'EdgeColor', ec);
                    %patch(tmr.video.nLocLaser(1,j,k) + spotr*cos(linspace(0,2*pi,npp)), tmr.video.nLocLaser(2,j,k) + spotr*sin(linspace(0,2*pi,npp)), cmap(c,:), 'Parent', ax, 'EdgeColor', cmap(c,:), 'FaceColor', 'none');

                end
             end
             if (breadcrumbs && isfield(tmr.video, 'micronDisplacement'))
                gridSpacing = 1000;
                offset = tmr.video.micronDisplacement(:,j);
%                 di = round(0.5*(tmr.video.tx(j+1)-tmr.video.tx(j))/tmr.dt); %offset to stage location midframe
%                 if (isfield(tmr.tracker, 'stageloc_lp'))
%                     offset = tmr.tracker.stageloc_lp(1:2, min(tmr.video.tind(end), tmr.video.tind(j)+di)); %lowpass calculated in playmovie
%                 else
%                     offset = tmr.tracker.stageloc(1:2, min(tmr.video.tind(end), tmr.video.tind(j)+di));
%                 end
                if (~isempty(crumb0))
                    crumb0 = repmat(offset,[1, size(crumb0,2)]) - crumb0;
                end
                offset = offset - floor(offset/gridSpacing)*gridSpacing;
                xL = ax.XLim;
                yL = ax.YLim;
%                 xmax = floor(tmr.micronsPerPixel*diff(ax.XLim)/(gridSpacing))*gridSpacing;
%                 ymax = floor(tmr.micronsPerPixel*diff(ax.YLim)/(gridSpacing))*gridSpacing;
                xmax = diff(xL)*tmr.micronsPerPixel;
                ymax = diff(yL)*tmr.micronsPerPixel;
                [xx,yy] = meshgrid(0:gridSpacing:xmax, 0:gridSpacing:ymax);
                xx = (xx + offset(1))/tmr.micronsPerPixel;
                yy = (yy + offset(2))/tmr.micronsPerPixel;
                valid = xx > xL(1) & xx < xL(2) & yy > yL(1) & yy < yL(2);
                xx = xx(valid); yy = yy(valid);
                plot (ax, xx, yy, 'c.');
                if (~isempty(crumb0))
                    plot(ax, size(im,2)/2 + crumb0(1,:)/tmr.micronsPerPixel, size(im,1)/2 +crumb0(2,:)/tmr.micronsPerPixel, 'co');
                end
            end
            if (overlaytime)
                if (flipy)
                    text (min(ax.XLim), max(ax.YLim), num2str(tmr.video.tx(j)-overlayt0,'%.1f'), 'Parent', ax, 'Color', 'c', 'VerticalAlignment', 'bottom');
                else
                    text (min(ax.XLim), min(ax.YLim), num2str(tmr.video.tx(j)-overlayt0,'%.1f'), 'Parent', ax, 'Color', 'c', 'VerticalAlignment', 'bottom');
                end
            end
            if (overlaybehavior && ~isempty(tmr.bsl) && isa(tmr.bsl, 'BehavioralStateList') && ~isempty(tmr.bsl.getStateFrame(j)) )
                txt = text (max(ax.XLim), min(ax.YLim), tmr.bsl.getStateFrame(j), 'Parent', ax, 'Color', 'c', 'VerticalAlignment', 'bottom','HorizontalAlignment', 'right');
            end
            if (marktail)
                text (tmr.video.tailloc(1,j), tmr.video.tailloc(2,j), 'T', 'Color', 'M', 'Parent', ax);
            end
            hold (ax, 'off');
            
        end
            
        function ax = plotFluorescence (tmr, ax, trange, ratrange, lpfilt, spotcolors)
            existsAndDefault('ratrange', [0.2 4.2]);
            existsAndDefault('lpfilt', false);
            existsAndDefault('spotcolors', {[0 1 1], [1 0 1], [1 1 0],[1 1 1]});
            nspots = length(tmr.neuron);
            if (nargin < 2 || isempty(ax))
                if (nspots == 1)
                    for j = 1:3
                        ax(j) = subplot(3, 1, j);
                    end
                else
                    for j = 1:nspots
                        ax(j) = subplot(nspots, 1, j);
                    end
                end
                
            end
            if (nargin >= 3)
                inds = tmr.tx >= trange(1) & tmr.tx <= trange(end);
            else
                inds = true(size(tmr.tx));
            end
           
            if (nspots == 1)            
                dd = {tmr.neuron.ratio_div_baseline(inds), tmr.neuron.greenrate(inds), tmr.neuron.redrate(inds)};
    %            dd_eb = {tmr.neuron.ratio_eb(inds), tmr.neuron.greenrate_eb(inds), tmr.neuron.redrate_eb(inds)};
                c = {[1 1 1], [0 1 0], [1 0 0]};
                for j = 1:3
                    axnum = min(j, length(ax));
                    if (axnum < j)
                        hold (ax(axnum), 'on');
                    end
                    if (lpfilt)
                        dd{j} = tmr.lpfilt(dd{j});
                    end
                    if (j == 1)
                        yy = dd{j};
                    else
                        yy = dd{j}/median(dd{j});
                    end
                    plot (ax(axnum), tmr.tx(inds), yy, 'k-', 'Color', c{j}, 'LineWidth', 2); 
    %                 ylim(ax(j), [0.2 10]*median(dd{j}));
    %                 axis(ax(j), 'tight');
                     ylim(ax(axnum), ratrange);
                    %shadedErrorPlot(tmr.tx(inds), dd{j}, dd_eb{j}, dd_eb{j}, c{j}, 'ax', ax(j));
                    axis(ax(axnum), 'off');
                end
                rax = 1;
            else
                c = spotcolors;
                for j = 1:nspots
                    axnum = min(j, length(ax));
                    if (axnum < j)
                        hold (ax(axnum), 'on');
                    end
                    dd = tmr.neuron(j).ratio_div_baseline(inds);
                    if (lpfilt)
                        dd = tmr.lpfilt(dd);
                    end
                    plot (ax(axnum), tmr.tx(inds), dd, 'k-', 'Color', c{mod(j-1,length(c))+1}, 'LineWidth', 2); %/median(dd)
                    ylim(ax(axnum),ratrange);
                    axis(ax(axnum), 'off');
                end
                rax = 1:nspots;
            end
             if (tmr.tracker.hasStatus && any(tmr.tracker.status.stimOn))
                 for k = 1:rax
                    hold (ax(k), 'on');

                    yl = get(ax(k), 'YLim');
                    si = find(diff(tmr.tracker.status.stimOn) > 0);
                    ei = find(diff(tmr.tracker.status.stimOn) < 0);
                    if (length(si) > length(ei) && si(end) > ei(end))
                        ei = [ei length(tmr.tracker.status.stimOn)];
                    end
                    if (length(ei) > length(si) && si(1) > ei(1))
                        si = [1 si];
                    end
                    for j = 1:min(length(si),length(ei))
                        patch(tmr.tx([si(j) ei(j) ei(j) si(j)]), [yl(1) yl(1) yl(2) yl(2)],'c', 'Parent', ax(k), 'EdgeColor', 'none');
                    end
                    ch = get(ax(k), 'Children'); set(ax(1), 'Children', ch(end:-1:1));
                    hold (ax(k), 'off');
                 end
            end
            
            set(ax, 'Color', 'k');
            
        end
        
        function tmr = correctRatioForBleaching(tmr, tbuffer, basepts)
            % fits the ratio to an exponential in time
            % finds the longest continuous tracking range to do this
            % and clips tbuffer from each end of this range
            % inspired by backcor.m by Vincent Mazet
            % basepts gives an initial list of baseline points, useful for
            % tracks with lots of activity
            nt = ~tmr.tracker.status.tracking; nt(1) = true; nt(end) = true;
            tnt = tmr.tx(nt);
            [~,I] = max(diff(tnt));
            tstart = tnt(I);
            tend = tnt(I+1);
              existsAndDefault('tbuffer', 5); 
             tbuffer = min(tbuffer, (tend-tstart)/4); %don't clip more than 50% of range
                 tstart = tstart + tbuffer;
                 tend = tend - tbuffer;
             
             inds = tmr.tx > tstart & tmr.tx < tend;
            
             xdata = tmr.tx(inds).';
             for k = 1:length(tmr.neuron)
                ydata = tmr.neuron(k).ratio(inds).';
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
              %        plot (xdata, ydata, xd(vrange),yd(vrange),xd, yf); title(num2str(j)); pause
                     xd = xd(vrange);
                     yd = yd(vrange);
                     s = sqrt(mean(dy(dy < 0).^2));
                    
                 end
                 
                 xd = xdata;
                 yd = ydata;
                 
                 if (~existsAndDefault('basepts', {}) || length(basepts) < k || (k > 1 && ~iscell(basepts) || isempty(basepts{k})));
                     x0 = polyfit(xd, log(yd), 1);        
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
             %    plot (xd, yd, xd,expfun(x0,xd), xd, yff); pause

                 %estimate standar deviation by subtracting off fit & only
                 %using negative values
                 dy = yd - yff;
                 s = sqrt(mean(dy(dy < 0).^2));
                 
                 

                 %fit to exponential using cost function that falls off after 2 
                 %positive standard deviation to 0 at 4 standard deviations
                 %and has 0 cost for anything 20 standard deviations
                 %outside range in either direction
                 cfun = @(u) assymFallOffLSQ(u, 2*s);
                 [x,yf] = fminArbCostFunction(@expfun, x, xd, yd, cfun, op);
             
                 %repeat
                 dy = yd - yf;
                 s = sqrt(mean(dy(dy < 0).^2));
                 cfun = @(u) assymFallOffLSQ(u, s);
                 [x,yf] = fminArbCostFunction(@expfun, x, xd, yd, cfun, op);
             
                 
                tmr.neuron(k).ratio_baseline = expfun(x, tmr.tx);
                tmr.neuron(k).ratio_div_baseline = tmr.neuron(k).ratio./tmr.neuron(k).ratio_baseline;
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
        
        
        function ah = overviewFigure(tmr, f, varargin)
            
            existsAndDefault('f', []);
            if (length(tmr.neuron) > 1)
                ah = tmr.overviewFigureMultiN(f, varargin{:});
                return;
            end
            trange = [];
            figTitle = 'overview figure';
            t0 = 0;
            varargin = assignApplicable(varargin);
             if (isempty(trange))
                trange = [tmr.tracker.startTime tmr.tracker.endTime];
            end
            
            trinds = tmr.tx > min(trange) & tmr.tx < max(trange);
            
           
            [ad, po] = blank8x10Figure(f, 'hspace', .5/11, 'topmargin', 0.35/11, 'bottommargin', 0.25/11, 'leftmargin', 1/11, 'rightmargin', 1/11, 'nrows', 5);
            po.fontsize = 10;
            po.axesopts = {'FontName', po.font, 'FontSize', po.fontsize, 'LineWidth', po.lineWidth/2, 'box', 'off'};

            po.redcolor = [1 0 0];
            po.greencolor = [0 .8 .5];
            annotation('textbox', [ad.lx2 ad.h0 ad.rx2+ad.w2-ad.lx2 1-ad.h0], 'String', figTitle, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'EdgeColor', 'none');
            an = 1;
            ah(an).pos = [ad.lx2, ad.h0 - ad.h, ad.rx2+ad.w2-ad.lx2, ad.h];
            ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
            ah(an).ax(2) = axes('Position', ah(an).pos, po.axesopts{:});
            ah(an).handles(1) = plot (ah(an).ax(2), tmr.tx(trinds) - t0, tmr.neuron.redrate(trinds), 'r', 'LineWidth', po.lineWidth, 'Color', po.redcolor);
            ah(an).handles(2) = plot (ah(an).ax(1), tmr.tx(trinds) - t0, tmr.neuron.greenrate(trinds), 'g', 'LineWidth', po.lineWidth, 'Color', po.greencolor);
            set(ah(an).ax, po.axesopts{:});
            set(ah(an).ax, 'XLim', trange - t0);
            
            
            
            mg = median(tmr.neuron.greenrate(trinds),'omitnan');
            mr = median(tmr.neuron.redrate(trinds),'omitnan');
            ylr = get(ah(an).ax(2), 'YLim'); 
            ylg = get(ah(an).ax(1), 'YLim'); 
            yls = [min([ylr/mr ylg/mg]) min(10,max([ylr/mr ylg/mg]))];
            set(ah(an).ax(1), 'YLim', yls*mg, po.axesopts{:}, 'YColor', po.greencolor);
            set(ah(an).ax(2), 'YLim', yls*mr, po.axesopts{:}, 'Color', 'none', 'YAxisLocation', 'right', 'YColor', po.redcolor, 'XTick', []);
%             set(ah(an).handles, 'LineWidth', po.lineWidth);
%             axis(ah(an).ax, 'tight');
            xlabel(ah(an).ax(1), 'time (s)');
            ylabel(ah(an).ax(2), 'red rate (/s)');
            ylabel(ah(an).ax(1), 'green rate (/s)');
            
            an = an+1;
            ah(an).pos = ah(an-1).pos - [0 ad.dh 0 0];
            ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
            ah(an).handles = plot (ah(an).ax, tmr.tx(trinds) - t0, tmr.neuron.ratio_div_baseline(trinds));
            set(ah(an).ax, po.axesopts{:});
            set(ah(an).handles, 'LineWidth', po.lineWidth, 'Color', po.color);
            axis(ah(an).ax, 'tight');
            ah = addStimShading (tmr, ah, an, t0);
            
            xlabel(ah(an).ax, 'time (s)');
            ylabel(ah(an).ax, 'ratio / baseline ratio');
            
            axl = {'x', 'y', 'z'};
            
            for k = 1:3
                an = an+1;
                hmod = 0.1*ad.h;
                if (k == 1)
                   ah(an).pos = ah(an-1).pos - [0 ad.dh-hmod 0 hmod];
                else
                    ah(an).pos = ah(an-1).pos - [0 ad.dh-hmod 0 0];
                end
                

                ti = max(min(trange),min(tmr.tx)):0.1:min(max(tmr.tx),max(trange));
                nl = interp1(tmr.tx, tmr.neuron.loc_filt(k,:), ti, 'linear');
                dnl = nl - lowpass1D(medfilt1(nl, 100),10);
                
                if (k == 3)
                    nl = nl*tmr.zflip;
                    dnl = dnl*tmr.zflip;
                end

                ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
                ah(an).ax(2) = axes('Position', ah(an).pos, po.axesopts{:});
                ah(an).handles(1) = plot (ah(an).ax(2), ti - t0, nl, 'b', 'LineWidth', po.lineWidth, 'Color', po.bluecolor);
                ah(an).handles(2) = plot (ah(an).ax(1), ti - t0, dnl, 'k', 'LineWidth', po.lineWidth, 'Color', po.color);
                set(ah(an).ax, po.axesopts{:});
                set(ah(an).ax, 'XLim', trange - t0);
                set(ah(an).ax(1), 'YColor', po.color);
                set(ah(an).ax(2), 'Color', 'none', 'YAxisLocation', 'right', 'YColor', po.bluecolor, 'XTick', []);
                xlabel(ah(an).ax(1), 'time (s)');
                ylabel(ah(an).ax(2), [axl{k} ' ($\mu$m)'], 'Interpreter', 'Latex');
                ylabel(ah(an).ax(1), ['local $\Delta$ ' axl{k} ' ($\mu$ m)'], 'Interpreter', 'Latex');
                ah = addStimShading (tmr, ah, an, t0);
            end

        end
        
        function ah = pathAndActivity(tmr, f, varargin)
            
            existsAndDefault('f', []);
            if (length(tmr.neuron) > 1)
                ah = tmr.overviewFigureMultiN(f, varargin{:});
                return;
            end
            trange = [];
            figTitle = 'path and activity';
            t0 = 0;
            varargin = assignApplicable(varargin);
            mseq = {'o', 's','p', 'v', 'd', '^', '<','h','>'};
            cseq = {[0 .8 0], [.8 0 0], [0 0 1],  [0.8 0.8 0], [1 0 1]};
            nsym = length(mseq)*length(cseq);
            
            
             if (isempty(trange))
                 sl = tmr.tracker.stageloc; sl1 = sl(:,find(all(isfinite(sl),1), 1, 'first'));
                 slend = sl(:,find(all(isfinite(sl),1), 1, 'last'));
                 
                 sd = sqrt(sum((tmr.tracker.stageloc - repmat(sl1, [1 size(tmr.tracker.stageloc,2)])).^2));
                 ind1 = find(sd > 50, 1, 'first'); %when stage moves 50 microns from start
                 sd = sqrt(sum((tmr.tracker.stageloc - repmat(slend, [1 size(tmr.tracker.stageloc,2)])).^2));
                 ind2 = find(sd > 50, 1, 'last'); %when stage moves 50 microns from start  
                 trange = tmr.tx([ind1 ind2]);
                 
%                trange = [tmr.tracker.startTime tmr.tracker.endTime];
            end
            
            trinds = tmr.tx >= min(trange) & tmr.tx <= max(trange);
            
            deltaTSym = 15;
            imi = ceil((max(trange)-min(trange))/(nsym*deltaTSym))*deltaTSym;
            minds = interp1(tmr.tx(trinds), (1:nnz(trinds))', min(trange):imi:max(trange), 'nearest');
            
           
            [ad, po] = blank8x10Figure(f, 'hspace', .5/11, 'topmargin', 0.35/11, 'bottommargin', 0.25/11, 'leftmargin', 1/11, 'rightmargin', 1/11, 'nrows', 2);
            po.fontsize = 10;
            po.axesopts = {'FontName', po.font, 'FontSize', po.fontsize, 'LineWidth', po.lineWidth/2, 'box', 'off'};

            
            annotation('textbox', [ad.lx2 ad.h0 ad.rx2+ad.w2-ad.lx2 1-ad.h0], 'String', figTitle, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'EdgeColor', 'none');
            
            
            an = 1;
            ah(an).pos = [ad.lx2, ad.h0 - ad.h, ad.rx2+ad.w2 - ad.lx2, ad.h];
            ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
            xd = tmr.neuron(1).loc_filt(1,trinds)/1000;
            yd = tmr.neuron(1).loc_filt(2,trinds)/1000;
            %ah(an).handles = plot (ah(an).ax, xd,yd, 'k-', 'LineWidth', po.lineWidth, 'Color', po.color); hold (ah(an).ax, 'on');
            ah(an).handles = plot(ah(an).ax, xd(1),yd(1), 'g.', 'MarkerSize', 20); hold(ah(an).ax, 'on');
            ah(an).handles = [ah(an).handles plot(ah(an).ax, xd(end),yd(end), 'rs', 'MarkerSize', 10)]; 
            for j = 1:length(minds)
                ci = mod(j-1, length(cseq)) + 1;
                mi = ceil(j / length(cseq));
                ah(an).handles = [ah(an).handles plot(ah(an).ax, xd(minds(j)),yd(minds(j)), 'k.', 'MarkerSize', 7, 'Marker', mseq{mi}, 'MarkerFaceColor', cseq{ci}, 'MarkerEdgeColor', 'none')];
                 if (j < length(minds))
                    ni = minds(j+1);
                else
                    ni = length(xd);
                end
                ah(an).handles = [ah(an).handles plot(ah(an).ax, xd(minds(j):ni),yd(minds(j):ni), 'k-', 'LineWidth', po.lineWidth, 'Color', cseq{ci})];
             
            end
            hold (ah(an).ax, 'off');
            
            axis(ah(an).ax, 'equal');
            axis(ah(an).ax, 'off');
            xl = get(ah(an).ax, 'XLim');           
            yl = get(ah(an).ax, 'YLim');
            
            aw = diff(xl);
            if (aw < 3)
                sbwidth = ceil(aw*2)/10;
            else
                sbwidth = ceil(aw/10);
            end
            
            x = xl + [sbwidth, -sbwidth];
            y = yl;
            mind = zeros([4 1]);
            for j = 1:2
                for k = 1:2
                    mind(j + 2*(k-1)) = min((x(j)-xd).^2 + (y(j)-yd).^2);
                end
            end
            [~,I] = max(mind);
            switch(I)
                case 1
                    sblocx = xl(1) + [0 sbwidth];
                    sblocy = [yl(1) yl(1)];
                    offset = 'above';
                case 2
                    sblocx = xl(2) + [-sbwidth 0];
                    sblocy = [yl(1) yl(1)];
                    offset = 'above';
                case 3
                    sblocx = xl(1) + [0 sbwidth];
                    sblocy = [yl(2) yl(2)];
                    offset = 'below';
                case 4
                     sblocx = xl(2) + [-sbwidth 0];
                     sblocy = [yl(2) yl(2)];
                    offset = 'below';
            end
            labeledBar(ah(an).ax, sblocx, sblocy, [num2str(sbwidth) ' mm'], offset, po.color, {'lineWidth', 5}, {'FontName', po.font, 'FontSize', po.fontsize});
           
            
            an = an+1;
            ah(an).pos = ah(an-1).pos;
            ah(an).pos(2) = ah(an).pos(2) - ad.dh;
            ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
            xd = tmr.tx(trinds);
            yd = tmr.neuron(1).ratio_div_baseline(trinds);
            %ah(an).handles = plot (ah(an).ax, xd,yd, 'k-', 'LineWidth', po.lineWidth, 'Color', po.color); hold (ah(an).ax, 'on');
            ah(an).handles = plot(ah(an).ax, xd(1),yd(1), 'g.', 'MarkerSize', 20); hold(ah(an).ax, 'on');
            ah(an).handles = [ah(an).handles plot(ah(an).ax, xd(end),yd(end), 'rs', 'MarkerSize', 10)]; 
            for j = 1:length(minds)
                ci = mod(j-1, length(cseq)) + 1;
                mi = ceil(j / length(cseq));
                ah(an).handles = [ah(an).handles plot(ah(an).ax, xd(minds(j)),max(yd), 'k.', 'MarkerSize', 7, 'Marker', mseq{mi}, 'MarkerFaceColor', cseq{ci}, 'MarkerEdgeColor', 'none')];
                if (j < length(minds))
                    ni = minds(j+1);
                else
                    ni = length(xd);
                end
                ah(an).handles = [ah(an).handles plot(ah(an).ax, xd(minds(j):ni),yd(minds(j):ni), 'k-', 'LineWidth', po.lineWidth, 'Color', cseq{ci})];
                
            end
            hold (ah(an).ax, 'off');
        end
        
        function ah = activityVsBehaviorFigure(tmr, f, varargin)
                        
            existsAndDefault('f', []);
            figTitle = 'activity vs. behavior';
            colorlist = {[0 1 1], [1 .5 0], [1 0 1], [0 1 0], [1 0 0], [0 0 1], [0 0 0]};
            varargin = assignApplicable(varargin);
            sn = tmr.bsl.stateNames;
            for j = 1:length(sn)
                valid(j) = any([tmr.bsl.stateVector.(sn{j})]);
            end
            sn = sn(valid);
            
            [ad, po] = blank8x10Figure(f, 'hspace', .5/11, 'topmargin', 0.35/11, 'bottommargin', 0.25/11, 'leftmargin', 1/11, 'rightmargin', 1/11, 'nrows', nnz(valid) + 2);
            po.fontsize = 10;
            po.axesopts = {'FontName', po.font, 'FontSize', po.fontsize, 'LineWidth', po.lineWidth/2, 'box', 'off'};
            annotation('textbox', [ad.lx2 ad.h0 ad.rx2+ad.w2-ad.lx2 1-ad.h0], 'String', figTitle, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'EdgeColor', 'none');
            
            for j = 1:length(sn)
                ah(j).ax = axes('Position', [ad.lx1, ad.h0 - ad.dh*j, ad.w1, ad.h], po.axesopts{:});
                ah(j).handles = tmr.plotValueForBehavioralState(tmr.tx, tmr.neuron(1).ratio_div_baseline, sn{j}, [], ah(j).ax, 'Color', colorlist{j});
                title(ah(j).ax, regexprep(sn{j}, '_',' '));
            end
            j = j+1;            
            ah(j).ax = axes('Position', [ad.lx1, ad.h0 - ad.dh*j, ad.w1, ad.h], po.axesopts{:});
            ah(j).handles = tmr.plotValueForBehavioralState(tmr.tx, tmr.neuron(1).ratio_div_baseline, {'forward'}, {}, ah(j).ax, 'Color', colorlist{1}); hold(ah(j).ax, 'on')
            ah(j).handles = tmr.plotValueForBehavioralState(tmr.tx, tmr.neuron(1).ratio_div_baseline, {'backward'}, {}, ah(j).ax, 'Color', colorlist{2});
            ah(j).handles = tmr.plotValueForBehavioralState(tmr.tx, tmr.neuron(1).ratio_div_baseline, {},{'backward', 'forward'},  ah(j).ax, 'Color', 'k');
            set([ah.ax], 'XLim', ah(j).ax.XLim, 'YLim', [ah(j).ax.YLim(1) max(tmr.neuron(1).ratio_div_baseline)], po.axesopts{:});
%            title(ah(j).ax, sn{j});
        end
            
        
        function h = plotValueForBehavioralState(tmr, tx, val, trueStateNames, falseStateNames, ax, varargin)
            existsAndDefault('falseStateNames', {});
            valid = logical(interp1(tmr.video.tx, double(tmr.bsl.getValidStates(trueStateNames, falseStateNames)), tx, 'nearest', false));
            si = find(diff(valid) > 0);
            if (valid(1))
                si = [1 si];
            end
            ei = find(diff(valid) < 0);
            if (valid(end))
                ei = [ei length(valid)];
            end
            existsAndDefault('ax', gca);
            np = ax.NextPlot;
            
            for j = 1:length(si)
                h(j) = plot(ax, tx(si(j):ei(j)), val(si(j):ei(j)), varargin{:}); %#ok<AGROW>
                ax.NextPlot = 'add';
            end
            ax.NextPlot = np;
        end
        
        function ah = overviewFigureMultiN(tmr, f, varargin)
            
            figTitle = 'overview figure';
            trange = [min(tmr.tx(tmr.tracker.status.tracking)) max(tmr.tx(tmr.tracker.status.tracking))];
            varargin = assignApplicable(varargin);
            if (isempty(trange))
                trange = [min(tmr.tx) max(tmr.tx)];
            end
            
            trinds = tmr.tx > min(trange) & tmr.tx < max(trange);
            
            existsAndDefault('f', []);
            
           
            [ad, po] = blank8x10Figure(f, 'hspace', .5/11, 'topmargin', 0.35/11, 'bottommargin', 0.25/11, 'leftmargin', 1/11, 'rightmargin', 1/11, 'nrows', 5);
            po.redcolor = [1 0 0];
            po.greencolor = [0 1 .5];
            annotation('textbox', [ad.lx2 ad.h0 ad.rx2+ad.w2-ad.lx2 1-ad.h0], 'String', figTitle, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'EdgeColor', 'none');
            an = 1;
            ah(an).pos = [ad.lx2, ad.h0 - ad.h, ad.rx2+ad.w2-ad.lx2, ad.h];
            ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
            ah(an).ax(2) = axes('Position', ah(an).pos, po.axesopts{:});

            rc = {[1 0 0], [.8 .2 .2], [.6 0 .4], [.7 .4 .1]};
            gc = {[0 1 0], [.2 .8 .2], [0 .6 .4], [.1 .7 .4]};
            ls = {'-', '-.', '--',':'};
            for k = 1:length(tmr.neuron)
                ah(an).handles(2*k-1) = plot (ah(an).ax(2), tmr.tx(trinds), tmr.neuron(k).redrate(trinds), 'r', 'LineWidth', po.lineWidth, 'Color', rc{mod(k-1,4)+1}, 'LineStyle', ls{mod(k-1,4)+1});
                ah(an).handles(2*k) = plot (ah(an).ax(1), tmr.tx(trinds), tmr.neuron(k).greenrate(trinds), 'g', 'LineWidth', po.lineWidth, 'Color', gc{mod(k-1,4)+1}, 'LineStyle', ls{mod(k-1,4)+1});
                hold (ah(an).ax(1), 'on');
                hold (ah(an).ax(2), 'on');
            end
            set(ah(an).ax, po.axesopts{:});
            set(ah(an).ax, 'XLim', trange);
            set(ah(an).ax(1), po.axesopts{:}, 'YColor', po.greencolor);
            set(ah(an).ax(2), po.axesopts{:}, 'Color', 'none', 'YAxisLocation', 'right', 'YColor', po.redcolor, 'XTick', []);
            xlabel(ah(an).ax(1), 'time (s)');
            ylabel(ah(an).ax(2), 'red rate (/s)');
            ylabel(ah(an).ax(1), 'green rate (/s)');
            
            an = an+1;
            ah(an).pos = ah(an-1).pos - [0 ad.dh 0 0];
            ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
            yd =  cat(1,tmr.neuron.ratio_div_baseline);
            ah(an).handles = plot (ah(an).ax, tmr.tx(trinds),yd(:,trinds));
            set(ah(an).ax, po.axesopts{:});
            set(ah(an).handles, 'LineWidth', po.lineWidth, 'Color', po.color);
            ratc = {[0 0 0], [0 0 1], [0 .5 .5], [.5 0 .5]};
            for k = 1:length(ah(an).handles)
                set (ah(an).handles(k), 'Color', ratc{mod(k-1,4)+1}, 'LineStyle', ls{mod(k-1,4)+1});
            end
            axis(ah(an).ax, 'tight');
            ah = addStimShading (tmr, ah, an);
            
            xlabel(ah(an).ax, 'time (s)');
            ylabel(ah(an).ax, 'ratio / baseline ratio');
            
            axl = {'x', 'y', 'z'};
            
            for k = 1:3
                an = an+1;
                hmod = 0.1*ad.h;
                if (k == 1)
                   ah(an).pos = ah(an-1).pos - [0 ad.dh-hmod 0 hmod];
                else
                    ah(an).pos = ah(an-1).pos - [0 ad.dh-hmod 0 0];
                end
                
                ah(an).ax = axes('Position', ah(an).pos, po.axesopts{:});
                ah(an).ax(2) = axes('Position', ah(an).pos, po.axesopts{:});
                for j = 1:length(tmr.neuron)
                    ti = min(trange):0.1:max(trange);
                    nl = interp1(tmr.tx, tmr.neuron(j).loc_filt(k,:), ti, 'linear');
                    dnl = nl - lowpass1D(medfilt1(nl, 100),10);
                    
                    if (k == 3)
                        nl = nl*tmr.zflip;
                        dnl = dnl*tmr.zflip;
                    end
                    ah(an).handles(2*j-1) = plot (ah(an).ax(2), ti, nl, 'b', 'LineWidth', po.lineWidth, 'Color', po.bluecolor, 'LineStyle', ls{mod(j-1,4)+1});
                    ah(an).handles(2*j) = plot (ah(an).ax(1), ti, dnl, 'k', 'LineWidth', po.lineWidth, 'Color', po.color,'LineStyle', ls{mod(j-1,4)+1});
                    hold (ah(an).ax(1), 'on'); hold (ah(an).ax(2), 'on');
                end
                 set(ah(an).ax, po.axesopts{:});
                set(ah(an).ax, 'XLim', trange);
                set(ah(an).ax(1), 'YColor', po.color);
                set(ah(an).ax(2), 'Color', 'none', 'YAxisLocation', 'right', 'YColor', po.bluecolor, 'XTick', []);

                xlabel(ah(an).ax(1), 'time (s)');
                ylabel(ah(an).ax(2), [axl{k} ' ($\mu$m)'], 'Interpreter', 'Latex');
                ylabel(ah(an).ax(1), ['local $\Delta$ ' axl{k} ' ($\mu$ m)'], 'Interpreter', 'Latex');
                ah = addStimShading (tmr, ah, an);
            end
        end
        function ah = addStimShading (tmr, ah, an, t0)
            existsAndDefault('t0', 0);
            if (all(ishandle(ah)) && all(isa(ah, 'matlab.graphics.axis.Axes')))
                ax = ah;
                if (nargin == 3)
                    t0 = an;
                end
                updatechildren = false;
            else
                ax = ah(an).ax(1);
                updatechildren = true;
            end
            if (tmr.tracker.hasStatus && any(tmr.tracker.status.stimOn))
                hold (ax, 'on');

                yl = get(ax, 'YLim');
                si = find(diff(tmr.tracker.status.stimOn) > 0);
                ei = find(diff(tmr.tracker.status.stimOn) < 0);
                if (length(si) > length(ei) && si(end) > ei(end))
                    ei = [ei length(tmr.tracker.status.stimOn)];
                end
                if (length(ei) > length(si) && si(1) > ei(1))
                    si = [1 si];
                end
%                 for j = 1:min(length(si),length(ei))
%                     ah(an).handles = [ah(an).handles plot(tmr.tx([si(j) ei(j)]), (0.02*yl(1) + 0.98*yl(2))*[1 1], 'c-', 'LineWidth', 3)];
%                 end
                xl = get(ax, 'Xlim');
                valid = tmr.tx(ei) - t0 > xl(1) & tmr.tx(si) - t0 < xl(2);
                si = si(valid);
                ei = ei(valid);
                
                if (any(valid))
                    oldhandles = get(ax, 'Children');
                    for j = 1:min(length(si),length(ei))
                        newhandles(j) =  patch(tmr.tx([si(j) ei(j) ei(j) si(j)]) - t0, [yl(1) yl(1) yl(2) yl(2)],'c', 'Parent', ax, 'EdgeColor', 'none'); %#ok<AGROW>
                    end
                    set(ax, 'Children', [oldhandles(:);newhandles(:)]);
                    if (updatechildren)
                        ah(an).handles = [ah(an).handles(:);newhandles(:)];
                    end
                end
                xlim(ax, xl);
                hold (ax, 'off');
            end
        end
        function [t, val] = getRaw(tmr, type, spotnum, ax)
            existsAndDefault('spotnum', 1);
            existsAndDefault('ax', 'x');
            switch (lower(ax))
                case 'x'
                    inds = tmr.tracker.raw.spot(spotnum).xinds;
                case 'y'
                    inds = tmr.tracker.raw.spot(spotnum).yinds;
                case 'z'
                    inds = tmr.tracker.raw.spot(spotnum).yinds;
                otherwise
                    warning ('tmr:badax', 'ax not recognized, using x inds');
                    inds = tmr.tracker.raw.spot(spotnum).xinds;
            end
            t = tmr.tracker.raw.et(inds);
            try
                val = tmr.tracker.raw.(type); val = val(inds);
            catch me %#ok<NASGU>
                fn = fieldnames(tmr.tracker.raw); for j = 1:length(fn), fn{j} = [' ' fn{j}]; end
                disp(['valid type names are' fn{:}]);            
            end
        end
        function [f, psd_data, psd_shot, psd_shot_fit] = calculatePSD (tmr, type, trange, spotnum, f, normtomean)
         %   function [f, psd_data, psd_shot, psd_shot_fit] = calculatePSD (tmr, type, trange, spotnum, f, normtomean)
            existsAndDefault('normtomean', true);
            df = 0.1;
            existsAndDefault('f', 0:df:200);
            
            df = median(diff(f));
            t = tmr.tx;
            existsAndDefault('spotnum', 1);
            
            rr = tmr.neuron(spotnum).redrate;
            rg = tmr.neuron(spotnum).greenrate;
            switch lower(type)
                case 'rawgreen'
                    [t,v] = tmr.getRaw('gphoton', spotnum, 'x');
                case 'rawred'
                    [t,v] = tmr.getRaw('rphoton', spotnum, 'x');
               case 'green'
                   v = rg;
                case 'red'
                   v = rr;
               case 'ratio'
                   v = tmr.neuron(spotnum).ratio;
                case 'ratio_div_baseline'
                   v = tmr.neuron(spotnum).ratio_div_baseline;
               otherwise
                   error ('unknown type: choices are rawgreen rawred green red ratio');
            end
            
            
           existsAndDefault('trange', [tmr.tracker.startTime tmr.tracker.endTime]);
           inds = t >= min(trange) & t < max(trange);
           t = t(inds);
           v = v(inds);
           if (length(inds) == length(rr))
               rr = rr(inds);
               rg = rg(inds);
           end
           deltat = median(diff(t));
           fs = 1./deltat;
           window = min(length(v), ceil(fs/df));
           noverlap = ceil(window/2);
           [psd_data, f] = pwelch(v - mean(v), window, noverlap, f, fs);
           if (normtomean)
               psd_data = psd_data / mean(v)^2;
           end
           if (nargout < 3)
               return;
           end
           
           if (any(strcmpi(type, {'rawgreen', 'rawred'})))
               lambda = sum(v) / (t(end) - t(1));
               psd_shot = lambda*deltat.^2 * ones(size(f));
               psd_shot_fit = psd_shot;              
           else
               rp = poissrnd(mean(rr)*deltat, size(t));
               gp = poissrnd(mean(rg)*deltat, size(t));
               fitstruct = photonSPPFtworates(t, rp, gp, tmr.d_loglambda);
               switch(lower(type))
                   case 'green'
                        vsim = fitstruct.lambda_green_s;
                   case 'red'
                       vsim = fitstruct.lambda_red_s;
                   case 'ratio'
                       vsim = fitstruct.ratio_s;
                   case 'ratio_div_baseline'
                       vsim = fitstruct.ratio_s/mean(tmr.neuron.ratio_baseline(inds));
                    otherwise
                        error ('unknown type: choices are rawgreen rawred green red ratio');
               end
               psd_shot = pwelch(vsim-mean(vsim), window, noverlap, f, fs);
               myfunl = @(xxx,xdata) log(xxx(1)) - log(xxx(2) + xdata.^2);
               xfitl = lsqcurvefit(myfunl, [1, 500], f, log(psd_shot));
               psd_shot_fit = xfitl(1)./(xfitl(2) + f.^2);
           end
           if (normtomean)
               psd_shot = psd_shot / mean(v)^2;
               psd_shot_fit = psd_shot_fit / mean(v)^2;
           end

        end
        
        function tmr = changeDirectory(tmr, newdirectory)
            %function tmr = changeDirectory(tmr, newdirectory)
        
            [~,f,ext] = fileparts(tmr.stage.filename);
            tmr.stage.filename = fullfile(newdirectory, [f ext]);
            [~,f,ext] = fileparts(tmr.tracker.filename);
            tmr.tracker.filename = fullfile(newdirectory, [f ext]);
            if (isfield(tmr.video, 'filename'))
                [~,f,ext] = fileparts(tmr.video.filename);
                tmr.video.filename = fullfile(newdirectory, [f ext]);
            end   
            
        end
        
        function saveBSL(tmr, filename, condensed)
            existsAndDefault('filename', fullfile(fileparts(tmr.video.filename), ['bsl-autosave-' datestr(now,30) '.csv']));
            existsAndDefault('condensed', true);           
            tmr.bsl.toCSV(filename, condensed);
        end
        
        function tmr = loadBSL(tmr, filename)
            if (nargin < 2)
                [filename,pathname] = uigetfile(fullfile(fileparts(tmr.video.filename), '*.csv'));
                if (~isequal(filename, 0))
                    filename = fullfile(pathname, filename);
                end
            end
            if (exist(filename, 'file'))
                try 
                    bsl = BehavioralStateList.fromCSV(filename);
                catch me
                    disp (me.getReport());
                    return;
                end
                tmr.bsl = bsl;
            end
        end
        
        function tmr = openVideo(tmr)
            try 
                if (tmr.video.vidreader.hasFrame())
                    return;
                end
                tmr.video.vidreader.CurrentTime = tmr.video.et(1);
                if (tmr.video.vidreader.hasFrame())
                    return;
                end
            catch 
                tmr.video.vidreader = VideoReader(tmr.video.filename);
            end
        end
                     
        
        function im = getVidFrame(tmr, vidframe)
        %function im = getVidFrame(tmr, vidframe)

            if (size(tmr.videocube,3) >= vidframe)
                im = tmr.videocube(:,:,vidframe);
                return;
            end
            tmr = tmr.openVideo();
            tmr.video.vidreader.CurrentTime = tmr.video.et(vidframe);
            im = tmr.video.vidreader.readFrame();    
        end
        
        function ah = overlayFluorescence(tmr, vidframe, fluorim)
            %function ah = overlayFluorescence(tmr, vidframe, fluorim)
            %vidframe default is 2 frames after start of tracking
            %fluorim default is to open a search directory
            existsAndDefault('vidframe', tmr.video.startFrame + 2); %make sure you don't get the flash
            if (nargout > 1)
                ah = [];
            end
            if (~existsAndDefault('fluorim',[]))
                [filename,pathname] = uigetfile({'*.jpg;*.bmp;*.tif;*.tiff;*.png'}, 'select an image', fullfile(fileparts(tmr.video.filename)));
                if (~isequal(filename, 0))
                    filename = fullfile(pathname, filename);
                end
                if (exist(filename, 'file'))
                    try
                        fluorim = imread(filename);
                    catch me
                        disp('could not load image');
                        disp(me.getReport());
                        return;
                    end
                else
                    disp('file not found');
                    return;
                end
            end
            fmicronsPerPixel = 5.5/40;
            ydir = -1;
            xdir = 1;
            offset = [-76 27];
            vidim = tmr.getVidFrame(vidframe);
            xx = xdir*(1:size(vidim, 2))*tmr.micronsPerPixel;
            yy = ydir*(1:size(vidim, 1))*tmr.micronsPerPixel;
            
            h = pcolor (xx - mean(xx) - offset(1),yy - mean(yy) - offset(2), adapthisteq(uint8(vidim))); shading flat; axis equal;
            if size(fluorim, 3) > 1
                fluorim = rgb2gray(fluorim);
            end
            ax = gca;
            ax2 = axes('Position', ax.Position, 'Color', 'none'); 
            xxf = (1:size(fluorim,2))*fmicronsPerPixel;
            yyf = (1:size(fluorim,1))*fmicronsPerPixel;
            hh = pcolor(ax2, xxf-mean(xxf), yyf-mean(yyf), adapthisteq(uint8(fluorim))); shading flat; 
            xl = [-1 1]*max(xxf);
            yl = [-1 1]*max(yyf);
            ax2.XLim = xl;
            ax2.YLim = yl;
            axis(ax2, 'equal');
            
            
            set(hh, 'FaceAlpha', 0.7);
            colormap(ax, 'gray');
            colormap(ax2, 'hot');
            
            ax.XLim = ax2.XLim;
            ax.YLim = ax2.YLim;
            axis(ax, 'off', 'ij');
            axis(ax2, 'off', 'ij');
            
            if (nargout > 0)
                ah.ax = [ax ax2];
                ah.h = [h hh];
            end
            
        end
        
        function tmr = locateTail(tmr, inds, varargin)
            updateInterval = 1000;
            rad = 12;
            varargin = assignApplicable(varargin);
             if (~existsAndDefault('inds', true(size(tmr.video.et))))
                 inds(tmr.video.startFrame + (-1:1)) = false;
                 inds(tmr.video.endFrame + (-1:1)) = false;
                 %inds(tmr.video.startFrame+1) = false;
                 
                if (~isfield(tmr.video, 'tailloc'))
                    tmr.video.tailloc = NaN([2 length(tmr.video.et)]);
                    tmr.video.tailclicked = false(1, length(tmr.video.et));
                end
            end
            if (islogical(inds));
                inds = find(inds);
            end
            [tmr.video.tailloc(:,inds), ~, tmr.video.tailclicked(inds)] = tmr.trackPoints(inds, updateInterval, rad, 'Tail', true);
%             while (any(any(~isfinite(tmr.video.tailloc(:,inds)))))
%                 notfound = ~any(isfinite(tmr.video.tailloc(:,inds)));
%                 inds = inds(notfound | inds  > find(notfound,1,'first') - updateInterval);
%                 updateInterval = max(100,updateInterval/2);
%                 [tmr.video.tailloc(:,inds), ~, tmr.video.tailclicked(inds)] = tmr.trackPoints(inds, updateInterval, rad, 'Tail', true);
%             end    
            
        end
        
        function tmr = locateHead(tmr, inds)
            %function tmr = locateHead(tmr, inds)
            
            if (~existsAndDefault('inds', 1:length(tmr.video.et)))
                if (~isfield(tmr.video, 'headloc'))
                    tmr.video.headloc = NaN([2 length(tmr.video.et)]);
                else
                    inds = find(~all(isfinite(tmr.video.headloc),1) & tmr.video.tx > tmr.tracker.startTime & tmr.video.tx < tmr.tracker.endTime);
                end
            end
            if (islogical(inds));
                inds = find(inds);
            end
            inds = inds(:)';
            figure(1); clf;
            for j = inds
                %pcolor(adapthisteq(uint8(tmr.getVidFrame(j)))); shading flat; axis equal; colormap gray;
                ax = tmr.drawFrame(j,gca,'overlayspot', false, 'overlaytime', true);  title('click on head - right click to finish');
                [x,y,button] = ginput(1);
                tmr.video.headloc(:,j) = [x;y];
                if (button > 1)
                    break;
                end
            end
            
        end
        function writeTailLoc(tmr, fname)
            %function tmr = locateHead(tmr, inds)
            existsAndDefault('fname', fullfile(fileparts(tmr.tracker.filename), 'tailloc.csv'));
            fid = fopen(fname, 'w');
            for j = 1:size(tmr.video.tailloc, 2)
                if (all(isfinite(tmr.video.tailloc(:,j))))
                    fprintf (fid, '%d,%.1f,%.1f, %d\n', j, tmr.video.tailloc(1,j), tmr.video.tailloc(2,j), double(logical(tmr.video.tailclicked(j))));
                end
            end
            fclose(fid);
        end
        function tmr = readTailLoc(tmr, fname)
            existsAndDefault('fname', fullfile(fileparts(tmr.tracker.filename), 'tailloc.csv'));
            data = importdata(fname);
            tmr.video.tailloc(:,data(:,1)) = data(:,2:3)';
            tmr.video.tailclicked(data(:,1)) = data(:,4)';
        end
        
        function writeHeadLoc(tmr, fname)
            %function tmr = locateHead(tmr, inds)
            existsAndDefault('fname', fullfile(fileparts(tmr.tracker.filename), 'headloc.csv'));
            fid = fopen(fname, 'w');
            for j = 1:size(tmr.video.headloc, 2)
                if (all(isfinite(tmr.video.headloc(:,j))))
                    fprintf (fid, '%d,%.1f,%.1f\n', j, tmr.video.headloc(1,j), tmr.video.headloc(2,j));
                end
            end
            fclose(fid);
        end
        function tmr = readHeadLoc(tmr, fname)
            existsAndDefault('fname', fullfile(fileparts(tmr.tracker.filename), 'headloc.csv'));
            data = importdata(fname);
            tmr.video.headloc(:,data(:,1)) = data(:,2:3)';
        end
        function [path, xform, selected] = trackPoints(tmr,inds, updateInterval, rad, locToClick, constrainToSpot)
            existsAndDefault('locToClick', 'center');
            existsAndDefault('updateInterval', Inf);
            existsAndDefault('rad', 12);
            existsAndDefault('constrainToSpot', false);
            path = NaN(2, length(inds));
            selected = false(1, length(path));
            %im = tmr.getVidFrame(inds(1));
            im = adapthisteq(uint8(tmr.getVidFrame(inds(1))));
            clf;
            pcolor(im); shading flat; colormap gray; axis equal
            title (['please click on ' locToClick]);
            [x,y] = ginput(1); 
            ctr = [x(1) y(1)];
            if (constrainToSpot)
                sll = tmr.video.spotLocLaser;
                iii = find(all(isfinite(sll)));
                sll(:,~all(isfinite(sll))) = interp1(iii, sll(:,iii)', find(~all(isfinite(sll))), 'nearest', 'extrap')';
                dts = sqrt(sum((sll(:,inds(1))-[x;y]).^2));
            end
                
            selected(1) = true;
        
            r = [max(1,ctr(1)-rad) max(1,ctr(2)-rad) min(size(im,2)-ctr(1)+rad,2*rad) min(size(im,1)-ctr(2)+rad,2*rad)]; 
            points = detectMinEigenFeatures(im, 'ROI', r);
            points = double(points.Location);
            in = (points(:,1) - ctr(1)).^2 + (points(:,2) - ctr(2)).^2 < rad^2;
            points = points(in,:);
            hold on;
            plot(ctr(1) + rad*cos(linspace(0,2*pi,20)),ctr(2) + rad*sin(linspace(0,2*pi,20)) ,'m-',points(:,1), points(:,2),'g.');
            hold off;
            ii = 1;
            pointTracker = vision.PointTracker('MaxBidirectionalError', 1);
            path(:,1) = ctr';
            initialize(pointTracker, points, im);
            lastUpdate = 1;
            problem = false(size(inds));
            for j = inds(2:end)
               
                % Initialize the tracker with the initial point locations and the initial
                % video frame.
               
                
                im = adapthisteq(uint8(tmr.getVidFrame(j)));
                oldPoints = points;
                if (~isempty(points))
                    [points,isFound] = pointTracker.step(im);
                    validPoints = points(isFound,:);
                    validOldPoints = oldPoints(isFound,:);
                else
                    validOldPoints = [];
                    validPoints = [];
                end
                if (length(validPoints) < 3)
                    hold off;
                     pcolor(im); shading flat; colormap gray; axis equal;
                    title(['Lost points - click on' locToClick '; right click to indicate problem in previous segment']);
                    [ctr(1),ctr(2), button] = ginput(1);
                    if (constrainToSpot)
                        dts = sqrt(sum((sll(:,j)-[ctr(1);ctr(2)]).^2));
                    end

                     if (button > 1)
                        problem(lastUpdate:(ii-1)) = true;
                    end
                    lastUpdate = ii;                     
                    r = [max(1,ctr(1)-rad) max(1,ctr(2)-rad) min(size(im,2)-ctr(1)+rad,2*rad) min(size(im,1)-ctr(2)+rad,2*rad)]; 
                    points = detectMinEigenFeatures(im, 'ROI', r);
                    points = double(points.Location);
                    selected(ii) = true;
                   
                else
                     [xform(ii), ~,points] = estimateGeometricTransform(validOldPoints, validPoints, 'similarity','MaxDistance', 3); 
                     [ctr(1),ctr(2)] = xform(ii).transformPointsForward(ctr(1),ctr(2));
                end
                in = (points(:,1) - ctr(1)).^2 + (points(:,2) - ctr(2)).^2 < rad^2;
                if (constrainToSpot)
                    pdts = sqrt((points(:,1) - sll(1,j)).^2 + (points(:,2) - sll(2,j)).^2);
                    in = in & abs(dts - pdts) < 2*rad;
                end
                points = points(in,:);
                if (length(points) < 6)
                    r = [max(1,ctr(1)-rad) max(1,ctr(2)-rad) min(size(im,2)-ctr(1)+rad,2*rad) min(size(im,1)-ctr(2)+rad,2*rad)]; 
                    points = detectMinEigenFeatures(im, 'ROI', r);
                    points = double(points.Location);
                    in = (points(:,1) - ctr(1)).^2 + (points(:,2) - ctr(2)).^2 < rad^2;
                    if (constrainToSpot)
                        pdts = sqrt((points(:,1) - sll(1,j)).^2 + (points(:,2) - sll(2,j)).^2);
                        in = in & abs(dts - pdts) < 2*rad;
                    end
                    points = points(in,:);
                end
                

                ii = ii +1;
                path(:,ii) = ctr';
                if (~isempty(points))
                    pointTracker.setPoints(points);
                end
                if (ii - lastUpdate > updateInterval)
                    hold off
                    pcolor(im); shading flat; colormap gray; axis equal; 
                    title(['Routine update - click on' locToClick '; right click to indicate problem in previous track']);
                    [ctr(1),ctr(2), button] = ginput(1);
                    if (constrainToSpot)
                        dts = sqrt(sum((sll(:,j)-[ctr(1);ctr(2)]).^2));
                    end
                    if (button > 1)
                        problem(lastUpdate:(ii-1)) = true;
                    end
                    lastUpdate = ii;
                    selected(ii) = true;
                end
                
                if (mod(ii,20) == 0)
                    pcolor(im); shading flat; colormap gray; axis equal; hold on
                    %plot(x([1:end 1]),y([1:end 1]),'m-',points(:,1), points(:,2),'g.');
                    plot(ctr(1),ctr(2),'ro', ctr(1) + rad*cos(linspace(0,2*pi,20)),ctr(2) + rad*sin(linspace(0,2*pi,20)) ,'m-',points(:,1), points(:,2),'g.');
                    hold off

                    pause(0.001);
                end
                if (ctr(1) < rad/4 || ctr(1) > size(im,2) - rad/4 || ctr(2) < rad/4 || ctr(2) > size(im,2) - rad/4)
                    disp('ge -- terminating!');
                    break;
                end
                
            end
            
            if (any(problem))
                
                si = find(diff(problem) > 0); if (problem(1)), si = [1 si]; end
                ei = find(diff(problem) < 0); if (problem(end)), ei = [ei length(problem)]; end
                
                for k = 1:length(si)
                    ii = si(k):ei(k); %inds(si(k):ei(k));
                    ui = min(updateInterval/5, length(ii)/5);
                    [p,x,s] = tmr.trackPoints(inds(ii), ui, rad, locToClick);
                    path(:,ii) = p;
                    xform(ii(2:end)) = x;
                    selected(ii) = s;
                   % [path(:,ii), xform(ii), selected(ii)] = tmr.trackPoints(inds(ii), ui, rad);
                end
            end
                
        end
        
        function tmr = annotateBouts(tmr, maxBoutDuration, varargin)
            existsAndDefault('maxBoutDuration', median(diff(sort([tmr.neuron(1).forwardMotionPeaks tmr.neuron(1).backwardMotionPeaks]))));
            if (isfield(tmr.neuron, 'forwardBout'))
                tmr.neuron = rmfield(tmr.neuron, 'forwardBout');
            end
            if (isfield(tmr.neuron, 'backwardBout'))
                tmr.neuron = rmfield(tmr.neuron, 'backwardBout');
            end
            
           
            
            for j = 1:length(tmr.neuron)
                %clear isolated peaks
                fmp = tmr.neuron(j).forwardMotionPeaks;
                d2next = min([Inf abs(diff(fmp))], [abs(diff(fmp)) Inf]);
                
                nearest = NaN(size(tmr.tx));
                deltat = Inf(size(nearest));
                for k = 1:length(tmr.neuron(j).forwardMotionPeaks)
                    dtt = abs(tmr.tx - tmr.neuron(j).forwardMotionPeaks(k));
                    nearest(dtt < deltat) = k;
                    deltat = min(dtt,deltat);
                end
                for k = 1:length(tmr.neuron(j).backwardMotionPeaks)
                    dtt = abs(tmr.tx - tmr.neuron(j).backwardMotionPeaks(k));
                    nearest(dtt < deltat) = -k;
                    deltat = min(dtt,deltat);
                end
                fields = {'ratio_div_baseline', 'redrate', 'greenrate', 'loc_filt', 'speed_filt', 'vel_rel_ts'};
                
          
                
                
                for k = 1:length(tmr.neuron(j).forwardMotionPeaks)
                    dtt = abs(tmr.tx - tmr.neuron(j).forwardMotionPeaks(k));
                    bout.inds = find(dtt < maxBoutDuration/2 & nearest == k);
                    bout.tx = tmr.tx(bout.inds);
                    for f = 1:length(fields)
                        bout.(fields{f}) = tmr.neuron(j).(fields{f})(:,bout.inds);
                    end
                    tmr.neuron(j).forwardBout(k) = furtherProcessBout(bout, tmr.neuron(j).forwardMotionPeaks(k));
                end
                for k = 1:length(tmr.neuron(j).backwardMotionPeaks)
                    dtt = abs(tmr.tx - tmr.neuron(j).backwardMotionPeaks(k));
                    bout.inds = find(dtt < maxBoutDuration/2 & nearest == -k);
                    bout.tx = tmr.tx(bout.inds);
                    for f = 1:length(fields)
                        bout.(fields{f}) = tmr.neuron(j).(fields{f})(:,bout.inds);
                    end
                    tmr.neuron(j).backwardBout(k) = furtherProcessBout(bout, tmr.neuron(j).backwardMotionPeaks(k));
                end
                fwdinds = [tmr.neuron(j).forwardBout.peakVelInd];
                
                                %Mirna added this here
                if isempty(tmr.neuron(j).forwardMotionPeaks)
                    tmpbout = tmr.neuron(j).backwardBout(1);
                    f = fieldnames(tmpbout)';
                    f{2,1} = {};
                    emptybout = struct(f{:});
                    tmr.neuron(j).forwardBout = emptybout;
                end
%                 
                if isempty(tmr.neuron(j).backwardMotionPeaks)
                    tmpbout = tmr.neuron(j).forwardBout(1);
                    f = fieldnames(tmpbout)';
                    f{2,1} = {};
                    emptybout = struct(f{:});
                    tmr.neuron(j).backwardBout = emptybout;
                end     
                
                %Mirna added this
%                  if isfield(tmr.neuron(j),'backwardBout')
                    bakinds = [tmr.neuron(j).backwardBout.peakVelInd];
%                  else
%                     bakinds = [];
%                  end

                
                                 
                 
                 
                allinds = [fwdinds bakinds];
                fwd = [true(size(fwdinds)) false(size(bakinds))];
                
                
                [allinds,III] = sort(allinds);
                fwd = fwd(III);
                for k = 1:length(tmr.neuron(j).forwardBout)
                    ni = find(allinds > tmr.neuron(j).forwardBout(k).peakVelInd, 1, 'first');
                    pi = find(allinds < tmr.neuron(j).forwardBout(k).peakVelInd, 1, 'last');
                    if (~isempty(ni))
                        tmr.neuron(j).forwardBout(k).timeToNext = diff(tmr.tx([tmr.neuron(j).forwardBout(k).peakVelInd allinds(ni)]));
                        tmr.neuron(j).forwardBout(k).isNextFwd = fwd(ni);
                    else
                        tmr.neuron(j).forwardBout(k).timeToNext = Inf;
                        tmr.neuron(j).forwardBout(k).isNextFwd = false;
                    end
                    if (~isempty(pi))
                        tmr.neuron(j).forwardBout(k).timeFromLast = diff(tmr.tx([allinds(pi) tmr.neuron(j).forwardBout(k).peakVelInd]));
                        tmr.neuron(j).forwardBout(k).wasLastFwd = fwd(pi);
                    else
                        tmr.neuron(j).forwardBout(k).timeFromLast = Inf;
                        tmr.neuron(j).forwardBout(k).wasLastFwd = false;
                    end
                end
                
                
                %Mirna added this
%                 if isfield(tmr.neuron(j),'backwardBout')
                    
                for k = 1:length(tmr.neuron(j).backwardBout)
                    ni = find(allinds > tmr.neuron(j).backwardBout(k).peakVelInd, 1, 'first');
                    pi = find(allinds < tmr.neuron(j).backwardBout(k).peakVelInd, 1, 'last');
                    if (~isempty(ni))
                        tmr.neuron(j).backwardBout(k).timeToNext = diff(tmr.tx([tmr.neuron(j).backwardBout(k).peakVelInd allinds(ni)]));
                        tmr.neuron(j).backwardBout(k).isNextFwd = fwd(ni);
                    else
                        tmr.neuron(j).backwardBout(k).timeToNext = Inf;
                        tmr.neuron(j).backwardBout(k).isNextFwd = false;
                    end
                    if (~isempty(pi))
                        tmr.neuron(j).backwardBout(k).timeFromLast = diff(tmr.tx([allinds(pi) tmr.neuron(j).backwardBout(k).peakVelInd]));
                        tmr.neuron(j).backwardBout(k).wasLastFwd = fwd(pi);
                    else
                        tmr.neuron(j).backwardBout(k).timeFromLast = Inf;
                        tmr.neuron(j).backwardBout(k).wasLastFwd = false;
                    end
                end
                               
                
                
                
                for k = 1:length(tmr.neuron(j).forwardBout)
                    if (k == 1 || ~tmr.neuron(j).forwardBout(k).wasLastFwd)
                        tmr.neuron(j).forwardBout(k).boutNumFromStart = 0;
                    else
                         tmr.neuron(j).forwardBout(k).boutNumFromStart = tmr.neuron(j).forwardBout(k-1).boutNumFromStart + 1;
                    end
                end
                 for k = length(tmr.neuron(j).forwardBout):-1:1;
                    if (k == length(tmr.neuron(j).forwardBout) || ~tmr.neuron(j).forwardBout(k).isNextFwd)
                        tmr.neuron(j).forwardBout(k).boutNumFromEnd = 0;
                    else
                         tmr.neuron(j).forwardBout(k).boutNumFromEnd = tmr.neuron(j).forwardBout(k+1).boutNumFromEnd + 1;
                    end
                 end
                for k = 1:length(tmr.neuron(j).backwardBout)
                    if (k == 1 || tmr.neuron(j).backwardBout(k).wasLastFwd)
                        tmr.neuron(j).backwardBout(k).boutNumFromStart = 0;
                    else
                         tmr.neuron(j).backwardBout(k).boutNumFromStart = tmr.neuron(j).backwardBout(k-1).boutNumFromStart + 1;
                    end
                end
                
                
                 for k = length(tmr.neuron(j).backwardBout):-1:1;
                    if (k == length(tmr.neuron(j).backwardBout) || tmr.neuron(j).backwardBout(k).isNextFwd)
                        tmr.neuron(j).backwardBout(k).boutNumFromEnd = 0;
                    else
                         tmr.neuron(j).backwardBout(k).boutNumFromEnd = tmr.neuron(j).backwardBout(k+1).boutNumFromEnd + 1;
                    end
                end
            end
            function bout = furtherProcessBout(bout, peakVelTime)
                
                pvi = find(bout.tx >= peakVelTime, 1, 'first');
                bout.peakVelInd = bout.inds(pvi);
                [bout.maxRatio, I] = max(bout.ratio_div_baseline);
                bout.peakRatInd = bout.inds(I);
                bout.minRatio = min(bout.ratio_div_baseline);
                %en.wikipedia.org/wiki/Quantile
                bout.ratVentiles = percentile(bout.ratio_div_baseline, 0.05:0.05:0.95);
                bout.ratRange = diff(bout.ratVentiles([1 end]));
                bout.ratMean = mean(bout.ratio_div_baseline);
                 bout.redVentiles = percentile(bout.redrate, 0.05:0.05:0.95);
                bout.redRange = diff(bout.redVentiles([1 end]));
                bout.redMean = mean(bout.redrate);
                oldBaseline = Inf;
                inliers = bout.ratio_div_baseline < bout.ratVentiles(1);
                baseline = mean(bout.ratio_div_baseline(inliers));
                reps = 0;
                while (abs(baseline-oldBaseline) > 0.01*baseline)
                    inliers = bout.ratio_div_baseline < baseline + 2*std(bout.ratio_div_baseline(inliers));
                    baseline = mean(bout.ratio_div_baseline(inliers));
                    reps = reps + 1;
                    if (reps > 20)
                        break;
                    end
                end
                bout.ratBaseline = baseline;
                
                fwdprogress = cumsum(bout.vel_rel_ts);
                [~,I] = find(diff(sign(lowpass1D(bout.vel_rel_ts, 0.03/tmr.dt))) ~= 0);
                I1 = I(find(I < pvi, 1, 'last'));
                I2 = I(find(I > pvi, 1, 'first'));
                I = [I1 I2];
                [~,II] = min(fwdprogress(I));
                bout.farthestBackInd = bout.inds(I(II));
                
                
                
            end
        end
        function [datamatrix] = getNeuronFieldStimAligned(tmr, fieldname, txi, centering,rownum, varargin)
             spotnum = 1;       
            existsAndDefault('centering', 'start');
            existsAndDefault('rownum', 1);
            
            validTimes = {};
            varargin = assignApplicable(varargin);
            
            
            switch(lower(centering))
                case {'start', 'begin'}
                    tc = tmr.tx([false diff(tmr.tracker.status.stimOn)>0]);
                case {'stop', 'end'}
                    tc = tmr.tx([false diff(tmr.tracker.status.stimOn)<0]);
                otherwise
                    warning ('did not recognize chosen centering -- using start time');
                    tc = tmr.tx([false diff(tmr.tracker.status.stimOn)>0]);
            end
            
            if (~isempty(validTimes))
                if (~iscell(validTimes))
                    validTimes = {validTimes};
                end
                valid = false(size(tc));
                for j = 1:length(validTimes)
                    valid = valid | (tc > validTimes{j}(1) & tc < validTimes{j}(2));
                end
            else
                valid = true(size(tc));
            end
            tc = tc(valid);
            
            datamatrix = NaN(length(tc),length(txi));
            if (isfield(tmr.neuron, fieldname))
                ydat = tmr.neuron(spotnum).(fieldname);
            else
                if (isfield(tmr.tracker, fieldname))
                    ydat = tmr.tracker.(fieldname);
                else
                    error ('field not found');
                end
            end
            
            if (length(txi) > 1 && median(diff(txi)) > tmr.dt)
                ydat = lowpass1D(ydat, 1/6*(median(diff(txi)/tmr.dt)));
            end
            
            
            for j = 1:length(tc)
                datamatrix(j,:) = interp1(tmr.tx, ydat(rownum,:), tc(j) + txi, 'linear', NaN)';
            end
            
        
        end
        
        function [datamatrix_forward, datamatrix_backward] = getNeuronFieldCentered(tmr, fieldname, txi, centering, rownum, varargin)
        %function [datamatrix_forward, datamatrix_backward] = getNeuronFieldCentered(tmr, fieldname, txi, centering, rownum, varargin)
            existsAndDefault('centering', 'vel');
            existsAndDefault('rownum', 1);
            if (length(tmr) > 1)
                datamatrix_forward = zeros([0 length(txi)]);
                datamatrix_backward = zeros([0 length(txi)]);
                for j = 1:length(tmr)
                    [dmf,dmb] = tmr(j).getNeuronFieldCentered(fieldname, txi, centering, rownum, varargin{:});
                    datamatrix_forward = [datamatrix_forward;dmf];
                    datamatrix_backward = [datamatrix_backward;dmb];
                end
                return;
            end
            spotnum = 1;       
            
            taxis = 'time';
            validTimes = {};
            varargin = assignApplicable(varargin);
            
            
            switch(lower(centering))
                case {'vel', 'velocity'}
                    cn = 'peakVelInd';
                case {'rat', 'ratio'}
                    cn = 'peakRatInd';
                case {'farthestbackind', 'fbi', 'farthestback', 'fb'}
                    cn = 'farthestBackInd';
                otherwise
                    warning ('did not recognize chosen centering -- using peakVelInd');
                    cn = 'peakVelInd';
            end
            cycleaxis = (strcmpi(taxis, 'cycle'));
                
            
           
                    
            if (isempty(tmr.neuron(spotnum).forwardBout))
                tc = [];
            else
                tc = tmr.tx([tmr.neuron(spotnum).forwardBout.(cn)]);
            end
            %tc = tmr.tx([tmr.neuron(spotnum).forwardBout]);
            if (~isempty(validTimes))
                if (~iscell(validTimes))
                    validTimes = {validTimes};
                end
                valid = false(size(tc));
                for j = 1:length(validTimes)
                    valid = valid | (tc > validTimes{j}(1) & tc < validTimes{j}(2));
                end
            else
                valid = true(size(tc));
            end
            tc = tc(valid);
           
            datamatrix_forward = NaN(length(tc),length(txi));
            if (isfield(tmr.neuron, fieldname))
                ydat = tmr.neuron(spotnum).(fieldname);
            else
                if (isfield(tmr.tracker, fieldname))
                    ydat = tmr.tracker.(fieldname);
                else
                    error ('field not found');
                end
            end
            
            if (length(txi) > 1 && ~cycleaxis)
                ydat = lowpass1D(ydat, 1/6*(median(diff(txi)/tmr.dt)));
            end
%             if (length(txi) > 1 && cycleaxis)
%                 ydat = lowpass1D(ydat, median(diff(tc))/6*(median(diff(txi)/tmr.dt)));
%             end
            if (cycleaxis)
                deltaP = median(diff(txi));
                px = 1:deltaP:(length(tc));
                treg = interp1(1:length(tc), tc, px, 'spline');
                ydreg = interp1(tmr.tx, ydat, treg);
                for j = 1:length(tc)                   
                    datamatrix_forward(j,:) = interp1(px, ydreg, j + txi, 'linear', NaN);    
                end
            else
                for j = 1:length(tc)
                    datamatrix_forward(j,:) = interp1(tmr.tx, ydat(rownum,:), tc(j) + txi, 'linear', NaN)';
                end
            end
            if (isempty(tmr.neuron(spotnum).backwardBout))
                tc = [];
            else
                tc = tmr.tx([tmr.neuron(spotnum).backwardBout.(cn)]);
            end
            if (~isempty(validTimes))
                if (~iscell(validTimes))
                    validTimes = {validTimes};
                end
                valid = false(size(tc));
                for j = 1:length(validTimes)
                    valid = valid | (tc > validTimes{j}(1) & tc < validTimes{j}(2));
                end
            else
                valid = true(size(tc));
            end
            tc = tc(valid);
            
            datamatrix_backward = NaN(length(tc),length(txi));
             if (cycleaxis)
                v = find(valid);          
                for j = 1:length(v)
                    bb = tmr.neuron(spotnum).backwardBout(v(j));
                    ix = interp1([-1 0 1], [bb.inds(1) bb.(cn) bb.inds(end)], txi, 'linear', NaN);
                    datamatrix_backward(j,:) = interp1(ydat(rownum,:), ix, 'linear', NaN)';
                end
            else
                for j = 1:length(tc)
                    datamatrix_backward(j,:) = interp1(tmr.tx, ydat(rownum,:), tc(j) + txi, 'linear', NaN)';
                end
             end
        
        end 
        
        
        function tmr = findMotionBouts (tmr, varargin)
            % function tmr = findMotionBouts (tmr, varargin)
            % varargin = vthresh([] - use standard deviation), minPeakSep(1 second)
            vthresh = [];
            minPeakSep = 1;
            lptime = 0.3;
            varargin = assignApplicable(varargin);
            if (~isfield(tmr.neuron, 'vel_rel_ts'))
                disp ('need to know tail-spot displacement');
                return;
            end
            for j = 1:length(tmr.neuron)
                nv = lowpass1D(tmr.neuron(j).vel_rel_ts, lptime/tmr.dt);
                nvslow = lowpass1D(tmr.neuron(j).vel_rel_ts, 10*lptime/tmr.dt);
                fwd_certain = nvslow > 0.15 * max(abs(nvslow));
                rev_certain = nvslow < -0.15 * max(abs(nvslow));
                
                fwd = fwd_certain;
                rev = rev_certain;
                
                %fill in fwd gaps
                si = find(diff(fwd_certain) < 0);
                ei = find(diff(fwd_certain) > 0);
                if (ei(1) < si(1))
                    ei = ei(2:end);
                end
             %   plot (tmr.tx, nv, 'k', tmr.tx(fwd), nv(fwd), tmr.tx(rev), nv(rev)); pause
                for k = 1:min(length(si),length(ei))
                    testi = (si(k):ei(k));
                    if (~any(rev_certain(testi)))
                        fwd(testi) = nvslow(testi) > 0;
                  %      plot (tmr.tx, nv, 'k', tmr.tx(fwd), nv(fwd), tmr.tx(rev), nv(rev), tmr.tx(testi), nv(testi), 'b.'); 
                    end
                end
                
                %fill in rev gaps
                si = find(diff(rev_certain) < 0);
                ei = find(diff(rev_certain) > 0);
                if (~isempty(ei) && ei(1) < si(1))
                    ei = ei(2:end);
                end
                for k = 1:min(length(si),length(ei))
                    testi = (si(k):ei(k));
                    if (~any(fwd_certain(testi)))
                        rev(testi) = nvslow(testi) < 0;
                      %  plot (tmr.tx, nv,'k', tmr.tx(fwd), nv(fwd), tmr.tx(rev), nv(rev), tmr.tx(testi), nv(testi), 'r.'); pause

                    end
                end
                
                %fill in gaps between two segments with simple threshold
                fwd = fwd | (~rev & nv > 0);
                rev = rev | (~fwd & nv < 0);
                
                if (isempty(vthresh))
                    vthreshold = std(nv(isfinite(nv)));
                else
                    vthreshold = vthresh;
                end
                [~,locsf] = findpeaks(nv, 'MinPeakHeight', vthreshold, 'MinPeakProminence', vthreshold, 'MinPeakDistance', minPeakSep/tmr.dt);
                locsf = locsf(fwd(locsf)); %require that average motion be forward at peak time
                tmr.neuron(j).forwardMotionPeaks = tmr.tx(locsf);
                [~,locsb] = findpeaks(-nv, 'MinPeakHeight', vthreshold, 'MinPeakProminence', vthreshold,  'MinPeakDistance', minPeakSep/tmr.dt);
                locsb = locsb(rev(locsb));
                tmr.neuron(j).backwardMotionPeaks = tmr.tx(locsb);
                    
                
               
                timeToForward = Inf(size(tmr.tx));
                timeToBackward = timeToForward;
                for k = 1:length(tmr.neuron(j).forwardMotionPeaks);
                    timeToForward = min(timeToForward, abs(tmr.tx - tmr.neuron(j).forwardMotionPeaks(k)));
                end
                for k = 1:length(tmr.neuron(j).backwardMotionPeaks);
                    timeToBackward = min(timeToBackward, abs(tmr.tx - tmr.neuron(j).backwardMotionPeaks(k)));
                end
                
                deltat = median(diff(sort([tmr.neuron(j).forwardMotionPeaks tmr.neuron(j).backwardMotionPeaks])));
                
%                 if (length(tmr.neuron(j).forwardMotionPeaks) > 1)
%                     dtf = max(median(diff(tmr.neuron(j).forwardMotionPeaks)), 10);
%                 else
%                     dtf = 10;
%                 end
                tmr.neuron(j).forward = timeToForward < timeToBackward & timeToForward < deltat;
                 
                if (length(tmr.neuron(j).backwardMotionPeaks) > 1)
                    dtb = max(median(diff(tmr.neuron(j).backwardMotionPeaks)), 10);
                else
                    dtb = 10;
                end
                tmr.neuron(j).backward = timeToBackward < timeToForward & timeToBackward < deltat;
                
                 plot (tmr.tx, nv, 'k', tmr.tx(tmr.neuron(j).forward), nv(tmr.neuron(j).forward), tmr.tx(tmr.neuron(j).backward), nv(tmr.neuron(j).backward),...
                     tmr.tx(locsf), nv(locsf), 'g.', tmr.tx(locsb), nv(locsb), 'r.'); title (['vthresh = ' num2str(vthreshold)]);
                
                
            end
        end
        
        function tmr = manualAnnotateBoutCurvature(tmr, fwdboutinds, bakboutinds)
            existsAndDefault('fwdboutinds', []);
            existsAndDefault('bakboutinds', []);
            
            if (isempty(fwdboutinds) && isempty(bakboutinds))
                fwdboutinds = 1:length(tmr.neuron(1).forwardBout);
                bakboutinds = 1:length(tmr.neuron(1).backwardBout);
            end
            for j = fwdboutinds
                tmr.neuron(1).forwardBout(j).sidePoints = [];
                tmr.neuron(1).forwardBout(j) = manualAnnotateBoutCurv(tmr.neuron(1).forwardBout(j));
            end
            for j = bakboutinds
                tmr.neuron(1).backwardBout(j).sidePoints = [];
                tmr.neuron(1).backwardBout(j) = manualAnnotateBoutCurv(tmr.neuron(1).backwardBout(j));
            end
            function bout = manualAnnotateBoutCurv(bout)
                [~,vind] = min(abs(tmr.video.tind - bout.peakRatInd));
                for k = 1:3
                    im(:,:,k) = tmr.getVidFrame(vind + k - 1);
                end
                pcolor(adapthisteq(uint8(mean(im,3)))); title('click from tail to head along side opposite tracking spot');
                xlim(size(im,2)*[.25 .75]); ylim(size(im,1)*[.25 .75]);
                colormap gray; axis equal; shading flat
                bout.sidePoints = getpts();
            end    
        end
        
    end
    
    methods (Static)
        
        function tmr = loadProcessedResult(fname)
            if (isdir(fname))
                if (exist(fullfile(fname, 'processed_result.mat'), 'file'))
                    fname = fullfile(fname, 'processed_result.mat');
                else
                    d = dir(fullfile(fname, '*.mat'));
                    if (length(d) == 1)
                        fname = fullfile(fname, d(1).name);
                    else
                        error ('directory does not contain single valid mat file');
                    end
                end
            else
                if (~exist(fname, 'file'))
                    error ('file not found');
                end
            end
            bob = load(fname, 'tmr');
            tmr = bob.tmr;
            tmr = tmr.changeDirectory(fileparts(fname));
                        
        end
        
        function tmr = processDirectory(dirname, varargin)
            redoAnalysis = false;
            makeMovie = false;
            processVideo = true;
            vidoptions = {};
            varargin = assignApplicable(varargin);
            if (nargin < 1)
                disp ('need dirname');
                return;
            end
            resultname = fullfile(dirname, 'processed_result.mat');
            if (~exist(resultname, 'file') || redoAnalysis)
                tmr = TrackingMicroscopeResult(dirname);
                if (~any(tmr.tracker.status.tracking))
                    warning ('status flag for tracking is uniformly 0');
                    disp (dirname)
                else
                    if (processVideo || makeMovie)
                        tmr = tmr.addVideo(varargin{:});
                         if (abs(tmr.video.timewarp - 1) > 5e-3)
                             warning ('timewarp is off by at least 0.5%');
                             disp('if video recording has correct frame rate, try using tmr.enforceVideoDeltaT(actualVideoDeltaT)');
                         end
                    end
                end
                try
                    disp ('saving result struct ...'); tic
                    save (resultname, 'tmr', '-v7.3');
                    disp ('saved!'); toc
                catch me
                    disp (me.getReport());
                    return;
                end
            else
                disp('loading result struct ...'); tic
                load(resultname, 'tmr');
                disp ('loaded!'); toc;
                try
                    disp ('loading video');
                    if (tmr.hasVideo)%#ok<NODEF>
                       tmr = tmr.loadVideo(); 
                       disp ('video loaded!'); toc;
                    else
                        if (makeMovie)
                            tmr = tmr.addVideo();
                        end
                    end
                catch me
                    disp (me.getReport());
                    return;
                end
            end
            try
                [~,figname] = fileparts(fileparts(tmr.tracker.filename));
                tmr.overviewFigure(1, 'figTitle', figname);
                saveas(1, fullfile(dirname, 'overview_fig.pdf'), 'pdf');
                savefig(1, fullfile(dirname, 'overview_fig.fig'));
                tmr.pathAndActivity(1, 'figTitle', figname);
                saveas(1, fullfile(dirname, 'pathAndActivity.pdf'), 'pdf');
                savefig(1, fullfile(dirname, 'pathAndActivity.fig'));
                if (makeMovie && tmr.hasVideo && any(tmr.tracker.status.tracking))
                    ops = tmr.moviePlaybackOptions_default();
                    tmr.playMovie([],'record',ops{:}, vidoptions{:});
                end
            catch me
                disp (me.getReport());
            end
        end
        
        function tmr = annotateExperiment(dirname)
            %function tmr = annotateExperiment(dirname)
            %function tmr = annotateExperiment(fname)
            %function tmr = annotateExperiment(tmr)
            %function tmr = annotateExperiment()
            if (nargin == 0)
                [filename,pathname] = uigetfile('processed_result.mat');
                if (isequal(filename, 0))
                    tmr = [];
                    return;
                end
                dirname = fullfile(filename, pathname);
            end
            if (isa(dirname, 'TrackingMicroscopeResult'))
                tmr = dirname;
                dirname = fileparts(tmr.tracker.filename);
            else
                if (exist(dirname, 'dir'))
                    filename = fullfile(dirname, 'processed_result.mat');
                    if (~exist(filename, 'file'))
                        [filename,pathname] = uigetfile(dirname, 'select result mat file');
                        filename = fullfile(pathname, filename);
                    end
                else
                    if (exist(dirname, 'file'))
                        filename = dirname;
                    else
                        disp (['I don''t know what to do with ' dirname]);
                        tmr = [];
                        return;
                    end
                end
                disp ('loading tmr'); tic;
                load(filename, 'tmr');
                dirname = fileparts(filename);
                tmr = tmr.changeDirectory(dirname); %#ok<NODEF>
                disp('done, loading video'); toc
                tmr = tmr.loadVideo();
                toc
                if (~existsAndDefault('tmr', []))
                    disp ('failed to load TMR');
                    return;
                end
               
            end
            disp ('experiment loaded -- let''s create a fluorescence overlay -- hit cancel to skip');
            [filename,pathname] = uigetfile({'*.jpg;*.bmp;*.tif;*.tiff;*.png'}, 'select starting image', dirname);
            if (~isequal(filename, 0))
                filename = fullfile(pathname, filename);
                try 
                    fluorim = imread(filename);
                    figure(1); clf;
                    ah = tmr.overlayFluorescence(tmr.video.startFrame + 2, fluorim);
                    u = ah.ax(1).Units;
                    ah.ax(1).Units = 'pixels';
                    im = frame2im(getframe(1, ah.ax(1).Position));
                    ah.ax(1).Units = u;
                    [~,~,ext] = fileparts(filename);
                    imwrite(im, fullfile(dirname, ['start_overlay-' datestr(now,30) '.' ext]));
                catch
                    disp ('failed to create overlay');
                end
            end
            try
                disp ('ok, let''s calibrate the stage to video magnification -- click on high contrast fixed spots to track -- NOT spots on the larva and NOT laser artifacts');
                pcolor(adapthisteq(uint8(median(tmr.videocube(:,:,1:ceil(size(tmr.videocube,3)/1000):size(tmr.videocube,3)),3,'omitnan')))); shading flat; colormap gray; title ('avoid these spots -- click to continue'); axis equal; ginput(1);
                
                tmr = tmr.manualRegistration();
                disp (['magnification = ' num2str(tmr.micronsPerPixel) ' microns per pixel']);
            
                disp ('ok, let''s find the tail location -- click around the middle of A7');
                tmr = tmr.locateTail();
                
                tmr = tmr.addTailSpotVec();
                
                tmr = tmr.findMotionBouts();
                tmr = tmr.annotateBouts();
                
                disp ('saving tmr with updated information to disk');
                 tmr.writeTailLoc();
                save(fullfile(dirname, 'processed_result.mat'), 'tmr', '-v7.3');
                disp ('done!');             
            catch me
                disp(me.getReport());
            end
            
%             disp ('ok, let''s do the behavioral annotation; close the window when you''re done and I''ll save the results');
%             try
%                 if (isempty(tmr.bsl) || ~isa(tmr.bsl, 'BehavioralStateList') || isempty(tmr.bsl.stateVector))
%                     tmr = tmr.initBSL();
%                 end
%                 if (isempty(tmr.bsl.nframes) || tmr.bsl.nframes < length(tmr.video.et))
%                     tmr.bsl = tmr.bsl.setLength(length(tmr.video.et));
%                 end
%                 tmr.bsl = tmr.bsl.addState({'forward'  'backward'  'left_turn'  'right_turn'  'squished'  'hunch'  'other'});
%                 tmr = tmr.behaviorGUI();
%                 tmr.saveBSL(fullfile(dirname, ['behavioral-annotations-' datestr(now,30) '.csv']));
%                 disp ('saving tmr with updated annotations to disk');
%                 save(fullfile(dirname, 'processed_result.mat'), 'tmr', '-v7.3');
%                 disp ('done!');
%             catch me
%                 disp(me.getReport());
%             end
            
             
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
%             [testf, testg, testh] = problem.objective(problem.x0);
%             testf
%             problem.objective(problem.x0 - 1e-6*testg./sqrt(sum(testg.^2))) -  problem.objective(problem.x0)
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
        function [ratio, sigma_ratio, lambda, sigma_lambda, thetaOut, cmat_thetaOut] = MLERatesOld (tedge, nr, ng, D1, D2)
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
            problem.x0 = [lg-lr;lg+lr];
            problem.solver = 'fminunc';
            problem.options = optimoptions(problem.solver);
            problem.options.GradObj = 'on';
            problem.options.Hessian = 'on';
            problem.options.Algorithm = 'trust-region';
%             [testf, testg, testh] = problem.objective(problem.x0);
%             testf
%             problem.objective(problem.x0 - 1e-6*testg./sqrt(sum(testg.^2))) -  problem.objective(problem.x0)
            thetaOut = reshape(fminunc(problem), [], 2);
            ratio = exp(thetaOut(:,1));
            lambda = exp(thetaOut*[.5 -.5;.5 .5]);
            
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
            sigma_lambda = 0.5*lambda.*(cmat_thetaOut*[1 1;2 -2; 1 1]).^(0.5);
            
            
            function [negLogP, negLogPgrad, negLogPHess] = probfun (theta)
                theta = reshape(theta, [],2);
                lambdaR = exp(0.5*(-theta(:,1) + theta(:,2)));
                lambdaG = exp(0.5*(theta(:,1) + theta(:,2)));
                %deltat_k = t_k - t_k-1; deltat(2:end)_k = t_k+1 - t_k
                
                dtheta = [[0 0]; diff(theta)];
                
                logP = sum(nr.*(0.5*(-theta(:,1) + theta(:,2))) - lambdaR.*deltat ...
                    + ng.*(0.5*(theta(:,1) + theta(:,2))) - lambdaG.*deltat) ...
                    - sum(dtheta(:,1).^2./(4*D1*deltat) + dtheta(:,2).^2./(4*D2*deltat));
                
                %deltat_k = t_k - t_k-1; deltat(2:end)_k = t_k+1 - t_k
                logPgrad(:,1) = 0.5*(ng - nr + deltat.*(lambdaR - lambdaG)) + 1./(2*D1)*(dtheta([2:end 1],1)./deltat([2:end 1]) - dtheta(:,1)./deltat);
                logPgrad(:,2) = 0.5*(ng + nr - deltat.*(lambdaR + lambdaG)) + 1./(2*D2)*(dtheta([2:end 1],2)./deltat([2:end 1]) - dtheta(:,2)./deltat);
                
                H12 = 0.25*(lambdaR-lambdaG).*deltat;
                H11 = -0.25*(lambdaR+lambdaG).*deltat - 1./(2*D1).*([1./deltat(2:end); 0] + 1./deltat);
                H22 = -0.25*(lambdaR+lambdaG).*deltat - 1./(2*D2).*([1./deltat(2:end); 0] + 1./deltat);
                n = size(theta, 1);
                
                negLogP = -logP;
                negLogPgrad = -logPgrad(:);
                negLogPHess = spdiags(-[[H12;H12], [H11;H22], [H12;H12]], [-n 0 n], 2*n, 2*n);
            end
           
        end

    end
    
    
end

