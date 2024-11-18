% MODIFIED FROM AKIHIRO'S motion_alignment.m SCRIPT

classdef motion_alignmentHS
%     Typical Use: 
%           ma = motion_alignment
% 
%     properties(Constant) % Q: How to define a global color scheme?
%        colors = viridis(numel(tsr.neuron)); 
%     end
%%
 methods
        % -*-*-*-*-*-*-*-*- Get functions (Compute & Get data output) -*-*-*-*-*-*-*-*-

        function [t_idx_act] = get_time_index_act(ma, activity, ti, tf)
            % Get the index of time for a selected initial and final time
            % ti: initial time (in sec)
            % tf: final time (in sec)
            
            t = activity.txfine;
            t_adj = false;
            
            if t_adj
              %  disp('time range adjustment');
            else
               % disp('Using the default time range.');
            end
            try
                t1 = find(t < ti);
                t_idx_act(1) = t1(end); 
            catch
               % warning('ti is starting time');
                t_idx_act(1) = 1;
            end
            try
                t2 = find(t < tf);
                t_idx_act(2) = t2(end);
            catch
                %warning('tf is ending time');
                t_idx_act(2) = length(t);
            end
        end

        function [t_idx_tsr] = get_time_index_tsr(ma, tsr, ti, tf)
            % Get the index of time for a selected initial and final time
            % ti: initial time (in sec)
            % tf: final time (in sec)
            
            t = tsr.tx;
            t_adj = false;
            
            if t_adj
               % disp('time range adjustment');
            else
               % disp('Using the default time range.');
            end
            try
                t1 = find(t < ti);
                t_idx_tsr(1) = t1(end); 
            catch
                warning('ti is the start time');
                t_idx_tsr(1) = 1;
            end
            try
                t2 = find(t < tf);
                t_idx_tsr(2) = t2(end);
            catch
                warning('tf is the ending time')
                t_idx_tsr(2) = length(t);
            end
        end

        function [loc, pos, dpos, ddpos, t] = get_filt_loc(ma, tsr, ti, tf, pass1D)
            % Get the x, y, z position of the neuron (tracker), 
            % Apply a low-pass filter to xyz to smooth out
            % Take the first and second derivative of the positions and
            % normalize ddpos
            
            t = tsr.tx(ti:tf);
            trange = true(size(t)); 
            loc = (tsr.neuron.loc(1:3,ti:tf));

%             pass1D = 1.425;

            pos = lowpass1D(loc(trange), pass1D/tsr.dt);
            dpos = diff(pos)./diff(t);
            ddpos = diff(pos,2)./max(diff(pos,2));
        end

        function [pos, params] = get_parameters(ma, tsr, thresh, ti, tf)
            % Default values:
            params.thresh = thresh; % need to adjust this depending on the data
            params.sp = 1e15;
            params.l = 1;
            params.nBreg = 5000;
            params.nInner = 1;
            params.epsilon = 1e-5;
            
            % ti and tf are the indices that correspond to the time array.
            % *** NOT the time (in seconds) ***
            pos = tsr.neuron.loc(1:3,ti:tf);
        end

        function [] = checkCounterInds(tsr)
        figure('Position',[200 90 550 900]); 
                subplot(4,1,1);
                plot(t(1:end-1), ds); hold on;
                plot(t(inds), ds(inds), 'ro');
                grid(); xlabel('time (sec)'); ylabel('ds');
                axis tight;

                subplot(4,1,2);
                plot(t(1:end-1), dds); hold on;
                plot(t(inds), dds(inds), 'ro'); 
                yline(params.thresh,'r--');
                grid; xlabel('time (sec)'); ylabel('dds');
                axis tight;
                
                % Plot: sx and sy
%                 figure();
                subplot(4,1,3);
                plot(t, loc(1,:)); grid on; hold on;
                plot(t, sx); xlabel('time (sec)'); ylabel('sx');
                plot(t(inds), loc(1, inds), 'ro'); axis tight;

%                 figure(); 
                subplot(4,1,4);
                plot(t, loc(2,:), t, sy); grid on;
                hold on; plot(t(inds), loc(2, inds), 'ro');
                xlabel('time (sec)'); ylabel('sy'); axis tight;
        end

        function [x, y] = patchRegion(ma, x1, x2, y1, y2, xOffset)
            x = (x1:0.1:x2)-xOffset;
            y1 = y1*ones(size(x));
            y2 = y2*ones(size(x));
            x = [x fliplr(x)];
            y = [y1 fliplr(y2)];
        end

        % I think this is mostly obsolete given tsr.findCounterMovements
        function [inds, s, sx, sy] = get_L1timing(ma, tsr, params, ti, tf, varargin)
            % Inputs:
            %   loc: position vector of the neuron
            %   thresh: threshold for dds
            %   t: limited time range (array)
            
            % Parameters for l1spline:
            %   sp: smoothing parameter
            %   l: split bregman smoothing parameter
            %   nBreg: number of plit bergman iterations
            %   nInner: number of inner iterations
            %   epsilon: acceptable relative error
            % Output:
            %   inds: indices of the furthest back point
            %   s   : dot product of direction and velocity
            %   sx, sy: smoothed x and y trajectories
            fig_check = false;
            sigma = 300;
            neuron_id = 1;
            assignApplicable(varargin);
            
            t = tsr.tx(ti:tf);
            loc = tsr.neuron.loc(1:3,ti:tf);
            
            % L1 smoothing of x and y
            sx = l1spline(loc(1,:), params.sp, params.l, params.nBreg, params.nInner, params.epsilon);
            sy = l1spline(loc(2,:), params.sp, params.l, params.nBreg, params.nInner, params.epsilon);
            
            % Direction vector
            tx = diff(sx) ./ sqrt(diff(sx).^2 + diff(sy).^2);
            ty = diff(sy) ./ sqrt(diff(sx).^2 + diff(sy).^2);
            
            % Get the dot product of the direction & velocity
            s = cumsum(tx.*diff(loc(1,:)) + ty.*diff(loc(2,:)));
            ds = deriv(s, sigma); % deriv(x, sigma): sigma might need to be adjusted
            dds = deriv(ds, 10);
            
            % Find the indices of the furthest back point along the
            % direction. Use the threshold to adjust the no. of indices. 
            inds = find(diff(ds > 0) > 0 & dds(1:end-1) > params.thresh);
            
            %disp(['neuron:', num2str(neuron_id), ', sigma = ', num2str(sigma),...
            %    ', #inds: ', num2str(size(inds,2))]);

            % ==================================================            
            % Check the threshold and the timing positions
            % ==================================================
            if fig_check==true
                figure('Position',[200 90 550 900]); 
                subplot(4,1,1);
                plot(t(1:end-1), ds); hold on;
                plot(t(inds), ds(inds), 'ro');
                grid(); xlabel('time (sec)'); ylabel('ds');
                axis tight;

                subplot(4,1,2);
                plot(t(1:end-1), dds); hold on;
                plot(t(inds), dds(inds), 'ro'); 
                yline(params.thresh,'r--');
                grid; xlabel('time (sec)'); ylabel('dds');
                axis tight;
                
                % Plot: sx and sy
%                 figure();
                subplot(4,1,3);
                plot(t, loc(1,:)); grid on; hold on;
                plot(t, sx); xlabel('time (sec)'); ylabel('sx');
                plot(t(inds), loc(1, inds), 'ro'); axis tight;

%                 figure(); 
                subplot(4,1,4);
                plot(t, loc(2,:), t, sy); grid on;
                hold on; plot(t(inds), loc(2, inds), 'ro');
                xlabel('time (sec)'); ylabel('sy'); axis tight;
            end
        end

        function [time_act, r, ratio_mean] = get_traces_n(ma, activity, t1, t2, tx, vs, traceinds)
%             assignApplicable(varargin);
             n = length(traceinds);

             ratio = lowpass1D(activity.voisets(vs).lpr(traceinds,:),2)./activity.voisets(vs).r0(traceinds);

             inds = find(activity.motion.tcross >= t1 & activity.motion.tcross <= t2);

             dt = mean(diff(activity.txfine)); % dt for activity.txfine is not the same for each time step
             nt = floor(tx/dt); % Needs to be an integer, but is not for activity 
             r = nan(n, length(inds), 2*nt+1);

             time_act = linspace(-tx, tx, length((inds(1)-nt:+inds(1)+nt)));
 
             for i=1:length(inds); 
                 [~,idx] = min(abs(activity.txfine - activity.motion.tcross(inds(i)))); 
                 indx_act(i) = idx; 
             end

             for i=1:length(indx_act)
                 intv_t1 = indx_act(i)-nt;
                 intv_t2 = indx_act(i)+nt;
                 for j=1:n
                     r(j,i,:) = ratio(j, intv_t1:intv_t2); %- ratio(indx(i)-nt); 
                 end
             end

             for j=1:n
                 ratio_mean(j,:) = mean(squeeze(r(j,:,:)),1);
             end
        end

        function [time_tsr, x, z] = get_trace_loc(ma, tsr, tx, s, loc, indx_tsr)
            dt = tsr.dt; nt = tx/dt;
            time_tsr = linspace(-tx, tx, length((indx_tsr(1)-nt:+indx_tsr(1)+nt)));

             for i=1:length(indx_tsr)
                 x(i,:) = s(1, indx_tsr(i)-nt:indx_tsr(i)+nt) - s(1,indx_tsr(i)-nt);
                 z(i,:) = loc(3, indx_tsr(i)-nt:indx_tsr(i)+nt) - loc(3,indx_tsr(i)-nt);
             end
        end

 end
end