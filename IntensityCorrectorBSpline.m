classdef IntensityCorrectorBSpline 
    % Calculates image values, derivatives, and Hessians 
    % at specified locations where photon was observed
    
    %alpha multiplies dwell or divides counts
    %corrected_rate = original_rate/alpha 

    properties (Constant)
     

        

    end
    
    properties
        
        npts = [6,6,4]; %num control points in x,y,z
        rate; %nx x ny x nz
        orig_rate; %nx x ny x nz
        
        step_size_termination = 1e-4;
%         max_dp_step = 1;
        dwell_min = 12.5e-9; %minimum dwell time to count towards correction
        remove_min_rate_voxels = true; %if template has a non-zero minimum - don't use voxels that are close to the minimum rate in analysis
        
        rateboxinds; %rate valid box -- used to clip template if it includes non-valid edges

        ratevalid; %where the clipped rate is valid -- used to remove non-valid template after box application
        validrate; %only valid rates (rate(ratevlid))
    end
    
    properties (Access = public)
        A; %NxM 
        N; %numel(ratevalid)
        M; %prod(npts)
        S; %s axes for spline representation - run from 0 to npts(j)-1
        J; %control point locations for spline interpolation - run from 1 to npts(j)
        Agrid; %to map fit values back to full grid (including where rate is not valid)
    end
    
    
    methods
       
        

        function ic = IntensityCorrectorBSpline(rate, npts) %immodel, photonloc, photonphase, photontime,  z_scale)

            ic.orig_rate = rate;
            if (nargin > 1)
                ic.npts = npts;
            end
            ic = ic.calculateBSPlineParameters();
        end
        

        %calculates parameters required for fitting
        function ic = calculateBSPlineParameters(ic)
            
            nd = ndims(ic.orig_rate);
        
            if (ic.remove_min_rate_voxels)
                mr = 1.01*min(ic.orig_rate(:),[],'omitnan');
            else
                mr = 0;
            end
            orvalid = ic.orig_rate > mr & isfinite(ic.orig_rate);
                
            for j = 1:nd
                ii = setdiff(1:nd,j);                 
                ic.rateboxinds{j} = find(squeeze(any(orvalid,ii)));
            end
            ic.rate = ic.orig_rate(ic.rateboxinds{:});
            ic.ratevalid = orvalid(ic.rateboxinds{:});
            ic.validrate = ic.rate(ic.ratevalid);

            ic.N = nnz(ic.ratevalid);
            ic.M = prod(ic.npts);
            ic.A = zeros(ic.N, ic.M); %num valid points
            ic.Agrid = zeros(numel(ic.rate), ic.M); %total points including invalids
            ic.S = zeros(numel(ic.rate), nd);
            ic.J = zeros(ic.M, nd);
            
            %create s grids 
            mygrid = cell([1, nd]);
            args = cell([1, nd]);
            for j = 1:ndims(ic.rate)
                args{j} = linspace(0,ic.npts(j)-1,size(ic.rate,j));
            end
            [mygrid{:}] = ndgrid(args{:});
            for j = 1:length(mygrid)
                ic.S(:,j) = mygrid{j}(:);
            end
            
            mygrid = cell([1, length(ic.npts)]);
            args = cell([1, length(ic.npts)]);
            for j = 1:length(ic.npts)
                args{j} = 1:ic.npts(j);
            end
            [mygrid{:}] = ndgrid(args{:});
            for j = 1:length(mygrid)
                ic.J(:,j) = mygrid{j}(:);
            end
            
            for d = 1:nd
                for j = 1:ic.npts(d)
                    BS{j,d} = IntensityCorrectorBSpline.chanCubicPolynomial(ic.S(:,d), j, ic.npts(d));
                end
            end

            %calculate spline weights
            for k = 1:ic.M %number of control points
                ic.Agrid(:,k) = 1;
                for j = 1:nd %number of dimensions
                    %ic.Agrid(:,k) = ic.Agrid(:,k).*IntensityCorrectorBSpline.chanCubicPolynomial(ic.S(:,j), ic.J(k,j), ic.npts(j));
                    ic.Agrid(:,k) = ic.Agrid(:,k).*BS{ic.J(k,j),j};
                end
            end
            ic.A = ic.Agrid(ic.ratevalid(:),:);
            
        end

        
        
         function [alpha,p] = correction_iterative_log(ic, counts, dwell)
             p = zeros([ic.M 1]);
             maxsteps = 100;
             %first clip to rate box
             counts = counts(ic.rateboxinds{:});
             dwell = dwell(ic.rateboxinds{:});

             %now remove inds with bad rate values
             counts = counts(ic.ratevalid);
             dwell = dwell(ic.ratevalid);


             %now only work with valid dwell times
              valid = (dwell(:) > ic.dwell_min) & isfinite(counts(:)) & isfinite(dwell(:));
              n = counts(valid);
              AA = ic.A(valid,:);
              sAA = sum(AA,1);

              ppvalid = sAA > min(sum(ic.Agrid,1))/10; %talks to less than 10% of the weight of the smallest element of the full grid
              
              AAA = AA(:,ppvalid);
%               if (any(~ppvalid))
%                   fprintf('%d/%d control points ignored as not sufficiently impacting sampled data\n', nnz(~ppvalid), ic.M);
%               end
              nbar0 = ic.validrate(valid).*dwell(valid);
              nn = nnz(valid);

              if rcond(AAA'*spdiags(nbar0,0,nn,nn)*AAA) < 1e-9
                ppvalid = sAA > min(sum(ic.Agrid,1))/4; %talks to less than 25% of the weight of the smallest element of the full grid
                fprintf('further restricted control points: %d/%d control points ignored as not sufficiently impacting sampled data\n', nnz(~ppvalid), ic.M);
              end
              AAA = AA(:,ppvalid);
              nbar0 = ic.validrate(valid).*dwell(valid);
              nn = nnz(valid);

              nmin = min(nbar0,[],'all')/100; 

             for j = 1:maxsteps
                 dp = zeros([ic.M 1]);
                 alpha = exp(AA*p);
                 nbar = nbar0.*alpha;
                 if (any(nbar < nmin))
                     fprintf('corrected photon rate too low - aborting\n');
                     break;
                 end
                 dp(ppvalid) = (AAA'*spdiags(nbar,0,nn,nn)*AAA)\AAA'*(n-nbar);
                 p = p + dp;
                 if (all(abs(dp) < ic.step_size_termination))
                     break;
                 end
                % max(abs(dp)) %for display
             end
             alpha = ones(size(ic.orig_rate));
             alpha(ic.rateboxinds{:}) = reshape(exp(ic.Agrid*p), size(ic.rate));
             p = reshape(p,ic.npts);
         end

         function [alpha,p] = calculate_correction(ic,counts,dwell)
             [alpha,p] = ic.correction_iterative_log(counts,dwell);
         end

         function [alpha,p] = calculate_correction_parallel(ic, count_arr, dwell_arr)
             %count_arr, dwell_arr are T x X x Y x Z (or T x X x Y etc) -
                     
                     D = parallel.pool.DataQueue;
            h = waitbar(0, 'Calculating corrections ...');
            afterEach(D, @nUpdateWaitbar);
            NN = size(count_arr,1); 
            q = 1;
            function nUpdateWaitbar(~)
                waitbar(q/NN, h);
                q = q + 1;
            end

             alpha = zeros(size(count_arr));
             p = zeros([size(count_arr,1) ic.npts]);
             tic
             parfor j = 1:size(count_arr,1)
                 [aa,pp] = ic.calculate_correction(squeeze(count_arr(j,:,:,:,:,:)), squeeze(dwell_arr(j,:,:,:,:,:)));
                 alpha(j,:,:,:,:,:) = aa;
                 p(j,:,:,:,:,:) = pp;
                 send(D,j);
             end
                 
         end

         function alpha = calculateAlpha(ic, p, maxcorrection)
             p = reshape(p, [ic.M 1]);
             alpha = ones(size(ic.orig_rate));
             alpha(ic.rateboxinds{:}) = reshape(exp(ic.Agrid*p), size(ic.rate));
             if (nargin >=3)
                 alpha(alpha > maxcorrection) = maxcorrection;
                 alpha(alpha < 1/maxcorrection) = 1/maxcorrection;
             end
         end
        
        
         function corrections = calc_frame_corrections(ic, framestub, frame_range)

            if (isfolder(framestub))
                dirname = framestub;
                d = dir(fullfile(dirname, 'registered_frames', '*_frame*.mat'));
            else
                dirname = '';
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
            fullframerange = [min(ind,[],'omitnan') max(ind,[],'omitnan')];

            existsAndDefault('frame_range',fullframerange);

            frame_range = [min(frame_range) max(frame_range)];


            frame = load(sprintf('%s%d.mat', framestub, frame_range(1)));



            corrections.frame_range = frame_range;
            corrections.xaxis = frame.xaxis;
            corrections.yaxis = frame.yaxis;
            corrections.ic = ic;

            corrections_tx = zeros([diff(frame_range)+1 1]);
            corrections_p = ones(diff(frame_range)+1,ic.M);

            frameinds = frame_range(1):frame_range(end);
            ts1 = tic;


            D = parallel.pool.DataQueue;
            h = waitbar(0, 'Calculating corrections ...');
            afterEach(D, @nUpdateWaitbar);
            NN = length(frameinds);
            p = 1;
            function nUpdateWaitbar(~)
                waitbar(p/NN, h);
                p = p + 1;
            end


            parfor j = 1:NN
                % tind = (1:4) + (j-1)*4;
                frame = load(sprintf('%s%d.mat', framestub, frameinds(j)));%, varsToLoad{:});
                corrections_tx(j) = mean(frame.time_edges, 'omitnan');
                [~,p] = ic.calculate_correction(sum(frame.non_rigid_red_counts,4), sum(frame.non_rigid_red_dwell,4));
                corrections_p(j,:) = p(:);
                send(D, j);
            end
            corrections.tx = corrections_tx;
            corrections.p = corrections_p;
            toc(ts1);
            close(h)

         end

    end
    
    methods(Static)

        
        
        function b = chanCubicPolynomial(s, index, n)
            if (nargin < 3)
                n = ceil(s(end)) + 1;
            end
            
            if (index == 1)
                b0 = IntensityCorrectorBSpline.chanCubicPolynomial(s,0);
            else
                if (index == n)
                    b0 = IntensityCorrectorBSpline.chanCubicPolynomial(s,n+1);
                else
                    b0 = zeros(size(s));
                end
            end
            s = s - index + 1;
            
            b = b0+ (2+s).^3/6 .*(s >= -2 & s < -1) + (2/3 - s.^2 - s.^3/2).*(s >= -1 & s < 0) + (2/3 - s.^2 + s.^3/2).*(s >= 0 & s < 1) + (2-s).^3/6 .*(s >= 1 & s < 2);
        end

        function corrections  = correct_directory (srcdir, register, saveresults, framerange, min_ctrl_pts,ctrl_pt_spacing)

            existsAndDefault('saveresults', false);
            existsAndDefault ('min_ctrl_pts', 3);
            existsAndDefault('ctrl_pt_spacing',20);
            existsAndDefault('framerange',[]);

            [~,timestamp,~] = fileparts(srcdir);
            load(fullfile(srcdir, register, [timestamp '_template.mat']), 'final_template');
            npts = max(ceil([diff(final_template.u([1 end])) diff(final_template.v([1 end])) diff(final_template.w([1 end]))]/ctrl_pt_spacing),min_ctrl_pts);
            tic
            ic = IntensityCorrectorBSpline(final_template.F,npts);
            toc
            corrections = ic.calc_frame_corrections(fullfile(srcdir, register), framerange);
            if (saveresults)
                save(fullfile(srcdir, register, [timestamp '_intensity_corrections.mat']), 'corrections', '-v7.3');
            end
        end

    end
    
end



