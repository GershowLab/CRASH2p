classdef RatioAnalysis
    %{
        Typical Use: 
            ra = RatioAnalysis;
    %}
    
    methods
        function ra = plot_ratio(ra, n_id, n_neuron, srcdir, t, trange, tsr, varargin)
            % NAME-VALUE inputs:
            %   'n_id': neuron id
            %   'srcdir': 
            %   't':
            %   'trange': [ti:tf]
            %   'tsr':
            %   'savefigure', true/[false]: whether to save a .png figure
            %   to disk
            
            % Default saving setting:
            dstdir = srcdir; %#ok<NASGU> % allow for more flexible saving - Rui 20201006
            savefigure = false;
            assignApplicable(varargin);
            colors = viridis(n_neuron);
            
            % Get the location of the neuron
            % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            loc = (tsr.neuron(n_id).loc(1:3,:)); % time is saved in tx         

            % Get photon counts and ratio
            % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            grn_rate = tsr.neuron(n_id).greenrate;
            red_rate = tsr.neuron(n_id).redrate;
            ratio = tsr.neuron(n_id).ratio_div_baseline;
            
            % --------------------------------------------
            % Plots of the photon counts, ratio, and xyz positions: 
            % --------------------------------------------
            fig = figure('Position', [50 50 800 600]); % [x y width height]

            subplot(3,1,1); 
            yyaxis right;
            plot(t(trange), red_rate(trange), 'r-', 'LineWidth', 1, 'DisplayName', 'red');
            ylabel('Photon Counts (red)'); hold on;
            yyaxis left;
            plot(t(trange), grn_rate(trange), 'g-', 'LineWidth', 1, 'DisplayName', 'green');
            legend(); xlabel('time (sec)'); ylabel('Photon Counts (green)'); grid on; hold off;

            subplot(3,1,2);
            plot(t(trange), ratio(trange), '-','Color', colors(n_id,:), 'LineWidth', 1, 'DisplayName', 'ratio');
            legend(); xlabel('time (sec)'); ylabel('Ratio'); grid on;

            subplot(3,1,3);
            % plot(t(trange), loc(1,trange), 'LineWidth', 2, 'DisplayName', 'x'); hold on;
            % plot(t(trange), loc(2,trange), 'LineWidth', 2, 'DisplayName', 'y');
            % plot(t(trange), loc(3,trange), 'LineWidth', 2, 'DisplayName', 'z');
            locx = loc(1,trange); locy = loc(2,trange); locz = loc(3,trange);
            yyaxis left;
            plot(t(trange), locx-locx(1),'b-', 'LineWidth', 2, 'DisplayName', 'x'); hold on;
            plot(t(trange), locy-locy(1),'k-', 'LineWidth', 2, 'DisplayName', 'y');
            ylabel('position (x,y)');
            
            yyaxis right;
            plot(t(trange), locz-locz(1), 'LineWidth', 2, 'DisplayName', 'z');
            ylabel('position (z)');
            
            legend(); xlabel('time (s)'); grid on;
            
            if strcmpi(srcdir(1),'G') 
                % only shorten data folder path when working directly on G:
                % drive and data structure is predictable - Rui 20201006
                title_name = strsplit(srcdir,'\');
                t_name = strcat(title_name(5),'\',title_name(6),'\',title_name(end),': n=',num2str(n_id));
            else
                t_name = [srcdir ': n=' num2str(n_id)];
            end
            sgtitle(t_name,'FontSize',12,'Interpreter','none');
            % set(titl, 'FontSize', 12);

            % --------------------------------------------
            % Save figure
            % --------------------------------------------
            if savefigure
%                 if strcmpi(srcdir(1),'G')
%                     % same reasoning as above - Rui 20201006
%                     save_name = strrep(t_name, '\', '_');
%                     save_name = strrep(save_name, ': ', '_');
%                     fpath = fullfile(dstdir, strcat(string(save_name), '.png'));
%                 else
                    % just save to dstdir with naming convention
                    % TIMESTAMP_FIGURETYPE.png regardless of whether working local or remote - Rui 20201006
                    timestamp = split(srcdir,filesep); %#ok<*UNRCH>
                    timestamp = timestamp{end};
                    fpath = fullfile(dstdir,[timestamp '_n=' num2str(n_id) '.png']);
%                 end
                saveas(fig,fpath);
                close(fig) % close figure window outside this function instead - Rui 20201006
            end
            % ======================================================
        end

        function ra = plot_one_trajectory(ra, n_id, srcdir, trange, tsr, varargin)
            % NAME-VALUE inputs:
            %   'n_id': neuron id
            %   'srcdir': 
            %   't':
            %   'trange': 
            %   'tsr':
            %   'savefigure', true/[false]: whether to save a .png figure
            %   to disk
            
            % Default saving setting:
            dstdir = srcdir; %#ok<NASGU>
            savefigure = false;
            assignApplicable(varargin);
            
            % Get the location of the neuron
            % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            loc = (tsr.neuron(n_id).loc(1:3,:)); % time is saved in tx         

            % --------------------------------------------
            % 2D plot of the trajectory
            % --------------------------------------------
            fig1 = figure('Position', [50 50 800 600]);
            plot(loc(1,trange), loc(2,trange)); hold on;
            plot(loc(1,1), loc(2,1), 'ro', 'LineWidth', 2);
            plot(loc(1,length(loc)), loc(2,length(loc)), 'rX', 'LineWidth', 2);
            grid on; xlabel('x'); ylabel('y');
            legend('trajectory', 'start', 'end');
            axis equal;
            
            if strcmpi(srcdir(1),'G')
            title_name = strsplit(srcdir,'\');
            t_name = strcat(title_name(5),'\',title_name(6),'\',title_name(end),': n=',num2str(n_id));
            else
                t_name = [srcdir ': n=' num2str(n_id)];
            end
            sgtitle(t_name,'FontSize',12,'Interpreter','none');

            % --------------------------------------------
            % Save figure
            % --------------------------------------------
            if savefigure
%                 if strcmpi(srcdir(1),'G')
%                     save_name = strrep(t_name, '\', '_');
%                     save_name = strrep(save_name, ': ', '_');
%                     fpath = fullfile(dstdir, strcat(string(save_name), '_map.png'));
%                 else
                    timestamp = split(srcdir,filesep); %#ok<*UNRCH>
                    timestamp = timestamp{end};
                    fpath = fullfile(dstdir,[timestamp '_n=' num2str(n_id) '.png']);
%                 end
                saveas(fig1,fpath);
                close(fig1)
            end
        end

        function ra = plot_ntrajectory(ra, srcdir, n_neuron, trange, tsr, varargin)
            % NAME-VALUE inputs:
            %   'n_id': neuron id
            %   'srcdir': 
            %   'trange': 
            %   'tsr':
            %   'savefigure', true/[false]: whether to save a .png figure
            %   to disk
            
            % Default saving setting:
            dstdir = srcdir; %#ok<NASGU>
            savefigure = false;
            assignApplicable(varargin);
            colors = viridis(n_neuron);
            
            fig1 = figure('Position', [50 50 800 600]);

            for n_id = 1:n_neuron
                % Get the location of the neuron
                % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
                loc = (tsr.neuron(n_id).loc(1:3,:)); % time is saved in tx         

                % --------------------------------------------
                % 2D plot of the trajectory
                % --------------------------------------------
                plot(loc(1,trange), loc(2,trange),'Color', colors(n_id,:), 'DisplayName', strcat('n=',num2str(n_id))); hold on;
                if n_id == n_neuron
                    plot(loc(1,1), loc(2,1), 'ro', 'LineWidth', 2, 'DisplayName', 'start');
                    plot(loc(1,length(loc)), loc(2,length(loc)), 'rX', 'LineWidth', 2, 'DisplayName', 'end');
                end
            end
            grid on; xlabel('x'); ylabel('y');
            legend('show');
            axis equal;
            
            
            % --------------------------------------------
            % Save figure
            % --------------------------------------------
            if strcmpi(srcdir(1),'G')
                title_name = strsplit(srcdir,'\');
                t_name = strcat(title_name(5),'\',title_name(6),'\',title_name(end));
            else
                t_name = srcdir;
            end

            if savefigure
%                 if strcmpi(srcdir(1),'G')
%                     save_name = strrep(t_name, '\', '_');
%                     save_name = strrep(save_name, ': ', '_');
%                     fpath = fullfile(dstdir, strcat(string(save_name), '_map.png'));
%                 else
                    timestamp = split(srcdir,filesep); %#ok<*UNRCH>
                    timestamp = timestamp{end};
                    fpath = fullfile(dstdir,[timestamp '_map.png']);
%                 end
                saveas(fig1, fpath);
                close(fig1)
            end

        end

        function ra = plot_nactivity(ra, srcdir, n_neuron, t, trange, tsr, varargin)
            % NAME-VALUE inputs:
            %   'n_id': neuron id
            %   'srcdir': 
            %   'trange': 
            %   'tsr':
            %   'savefigure', true/[false]: whether to save a .png figure
            %   to disk
            
            % Default saving setting:
            dstdir = srcdir; %#ok<NASGU>
            savefigure = false;
            c_range = false;
            assignApplicable(varargin);
            if c_range
                fig1 = figure('Position', [50 50 600 600]);
            else
                fig1 = figure('Position', [50 50 1200 600]);
            end

            colors = viridis(n_neuron);
            for n_id = 1:n_neuron
                % Get the ratio
                % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
                ratio = tsr.neuron(n_id).ratio_div_baseline;
      
                % --------------------------------------------
                % Plot activity
                % --------------------------------------------
                plot(t(trange), ratio(trange), '-', 'Color', colors(n_id,:), 'LineWidth', 1, 'DisplayName', strcat('n=',num2str(n_id)));
                hold on;
            end
            xlabel('time (sec)'); 
            ylabel('Ratio'); 
            grid on;
            legend('show');
            if c_range
                axis tight;
            end
            
            % --------------------------------------------
            % Save figure
            % --------------------------------------------
            if strcmpi(srcdir(1),'G')
                title_name = strsplit(srcdir,'\');
                t_name = strcat(title_name(5),'\',title_name(6),'\',title_name(end));
            else
                t_name = srcdir;
            end
            title(t_name,'Interpreter','none');
            
            if savefigure
%                 if strcmpi(srcdir(1),'G')
%                     save_name = strrep(t_name, '\', '_');
%                     save_name = strrep(save_name, ': ', '_');
%                     if c_range
%                         fpath = fullfile(dstdir, strcat(string(save_name), '_ratio_crange.png'));
%                     else
%                         fpath = fullfile(dstdir, strcat(string(save_name), '_ratio.png'));
%                     end
%                 else
                    timestamp = split(srcdir,filesep); %#ok<*UNRCH>
                    timestamp = timestamp{end};
                    if c_range
                        fpath = fullfile(dstdir, [timestamp '_ratio_crange.png']);
                    else
                        fpath = fullfile(dstdir, [timestamp '_ratio.png']);
                    end
%                 end
                saveas(fig1,fpath);
                close(fig1)
            end

        end

        function ra = rel_locs(ra, srcdir, n_neuron, t, trange, tsr, varargin)
            % NAME-VALUE inputs:
            %   'n_id': neuron id
            %   'srcdir': 
            %   'trange': 
            %   'tsr':
            %   'savefigure', true/[false]: whether to save a .png figure
            %   to disk
            
            % Default saving setting:
            dstdir = srcdir; %#ok<NASGU>
            savefigure = false;
            assignApplicable(varargin);
            
            loc_all = NaN(3,length(t),n_neuron); % (axis-->x,y,z, value, neuron#)
            
            for n_id = 1:n_neuron
                % Get the location of neurons
                loc_all(:,:,n_id) = tsr.neuron(n_id).loc(:,trange);
            end
            
            % Compute the center of mass (mean)
            loc_com = mean(loc_all,3);
            
            % --------------------------------------------
            % Plot stuff
            % --------------------------------------------
            
            fig = figure('Position', [50 50 1200 200+200*n_neuron]); % [x y width height]
    
            n_subplot = n_neuron + 1;
            
            for i =1:n_neuron
                loc_i = tsr.neuron(i).loc(:,:);
                
                subplot(n_subplot,1,i);
                title(['neuron ', num2str(i)], 'FontSize', 10);

%                 title('neuron ', num2str(i));
%                 yyaxis left;
%                 disp(['t(trange:)', size(t(trange))]);
%                 disp(['loc:)', size(loc_i(1,trange)-loc_com(1,trange))]);
                plot(t(trange), loc_i(1,trange)-loc_com(1,trange), 'b-', 'LineWidth', 1, 'DisplayName', '$x_i$-$\bar{x}$');
                hold on;
                plot(t(trange), loc_i(2,trange)-loc_com(2,trange), 'k-', 'LineWidth', 1, 'DisplayName', '$y_i$-$\bar{y}$');
%                 ylabel('position(s) (x,y)', 'FontSize', 12);
                
%                 yyaxis right;
                plot(t(trange), loc_i(3,trange)-loc_com(3,trange), '-', 'LineWidth', 1, 'DisplayName', '$z_i$-$\bar{z}$');
                ylabel('position', 'FontSize', 12); grid on; hold off;
                
                leg = legend('show', 'Location', 'northwest', 'FontSize', 12); set(leg, 'Interpreter', 'LaTeX');
                clear leg
            end
            
            subplot(n_subplot,1,n_subplot);
            if n_neuron==2
%                 yyaxis left;
                plot(t(trange), loc_all(1,trange,1)-loc_all(1,trange,2), 'b-', 'LineWidth', 1, 'DisplayName', '\Deltax'); hold on;
                plot(t(trange), loc_all(2,trange,1)-loc_all(2,trange,2), 'k-', 'LineWidth', 1, 'DisplayName', '\Deltay');
%                 ylabel('positions (x,y)', 'FontSize', 12);

%                 yyaxis right;
                plot(t(trange), loc_all(3,trange,1)-loc_all(3,trange,2), '-', 'LineWidth', 1, 'DisplayName', '\Deltaz');
                ylabel('position', 'FontSize', 12);
                title(['Distance between ', num2str(n_neuron), ' neurons'], 'FontSize', 12);
            else 
%                 yyaxis left;
                plot(t(trange), loc_com(1,trange), 'b-', 'LineWidth', 1, 'DisplayName', 'x_{COM}'); hold on;
                plot(t(trange), loc_com(2,trange), 'k-', 'LineWidth', 1, 'DisplayName', 'y_{COM}');
%                 ylabel('positions (x,y)', 'FontSize', 12);

%                 yyaxis right;
                plot(t(trange), loc_com(3,trange), '-', 'LineWidth', 1, 'DisplayName', 'z_{COM}');
                ylabel('position', 'FontSize', 12);
                title(['Center of mass of ', num2str(n_neuron), ' neurons'], 'FontSize', 12);
            end 
            xlabel('time', 'FontSize', 12); 
            grid on;
            legend('show', 'Location', 'northwest', 'FontSize', 10);
            
            % --------------------------------------------
            % Save figure
            % --------------------------------------------
            if strcmpi(srcdir(1),'G')
                title_name = strsplit(srcdir,'\');
                t_name = strcat(title_name(5),'\',title_name(6),'\',title_name(end));
            else
                t_name = srcdir;
            end
            sgtitle(t_name,'Interpreter','none');
            
            if savefigure
%                 if strcmpi(srcdir(1),'G')
%                     save_name = strrep(t_name, '\', '_');
%                     save_name = strrep(save_name, ': ', '_');
%                     fpath = fullfile(dstdir, strcat(string(save_name), '_locs.png'));
%                 else
                    timestamp = split(srcdir,filesep); %#ok<*UNRCH>
                    timestamp = timestamp{end};
                    fpath = fullfile(dstdir,[timestamp '_locs.png']);
%                 end
                saveas(fig,fpath);
                close(fig)
            end
        end
        
        function ra = plot_piezo_stage(ra, srcdir, t, trange, tsr, varargin)
            % NAME-VALUE inputs:
            %   'n_id': neuron id
            %   'srcdir': 
            %   'trange': 
            %   'tsr':
            %   'savefigure', true/[false]: whether to save a .png figure
            %   to disk
            
            % Default saving setting:
            dstdir = srcdir; %#ok<NASGU>
            savefigure = false;
            assignApplicable(varargin);
            
            % --------------------------------------------
            % Plot stuff
            % --------------------------------------------
            
            fig = figure('Position', [50 50 800 600]); % [x y width height]
            
            loc_piezo = tsr.tracker.piezo(trange);
            loc_stage = tsr.tracker.stageloc(:,trange);

            subplot(2,1,1);
            title(['Stage loc'], 'FontSize', 12);

            yyaxis left;
            plot(t(trange), loc_stage(1,:), 'b-', 'LineWidth', 1, 'DisplayName', '$x$');
            hold on;
            plot(t(trange), loc_stage(2,:), 'k-', 'LineWidth', 1, 'DisplayName', '$y$');
            ylabel('position(s) (x,y)', 'FontSize', 12);

            yyaxis right;
            plot(t(trange), loc_stage(3,:), '-', 'LineWidth', 1, 'DisplayName', '$z$');
            ylabel('position (z)', 'FontSize', 12); grid on; hold off;

            leg = legend('show', 'Location', 'southwest', 'FontSize', 12); set(leg, 'Interpreter', 'LaTeX'); clear leg;
            xlabel('time');
            
            subplot(2,1,2);
            title(['Piezo and Stage loc (z)'], 'FontSize', 12);
            yyaxis left;
            plot(t(trange), loc_piezo, 'b-', 'LineWidth', 1, 'DisplayName', 'piezo'); hold on;
            ylabel('piezo loc (z)', 'FontSize', 12);

            yyaxis right;
            plot(t(trange), loc_stage(3,:), '-', 'LineWidth', 1, 'DisplayName', 'stage');
            ylabel('stage loc (z)', 'FontSize', 12);

            xlabel('time', 'FontSize', 12); 
            grid on;
            leg = legend('show', 'Location', 'southwest', 'FontSize', 10); set(leg, 'Interpreter', 'LaTeX'); clear leg;
            
            % --------------------------------------------
            % Save figure
            % --------------------------------------------
            if strcmpi(srcdir(1),'G')
                title_name = strsplit(srcdir,'\');
                t_name = strcat(title_name(5),'\',title_name(6),'\',title_name(end));
            else
                t_name = srcdir;
            end
            sgtitle(t_name,'Interpreter','none');
            
            if savefigure
%                 if strcmpi(srcdir(1),'G')
%                     save_name = strrep(t_name, '\', '_');
%                     save_name = strrep(save_name, ': ', '_');
%                     fpath = fullfile(dstdir, strcat(string(save_name), '_piezo.png'));
%                 else
                    timestamp = split(srcdir,filesep); %#ok<*UNRCH>
                    timestamp = timestamp{end};
                    fpath = fullfile(dstdir,[timestamp '_piezo.png']);
%                 end
                saveas(fig,fpath);
                close(fig)
            end
        end
        
    end
end