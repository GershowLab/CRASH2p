classdef VRFileAccess
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
         galvoDeltaT = 0; %added to galvo time to fine tune sampling time (in seconds, typical values probably <~ us)
    end
    
    properties %(Access = protected)
        vrfc
        filename
        rollover;
        maxBytes = 2E9;
        counterData = [];
        hasCounterData = false;
    end
    properties (Transient = true)
        fid = -1;
    end
    
    methods
        function vrfa = setMaxBytes(vrfa, maxBytes)
            vrfa.maxBytes = maxBytes;
        end
        
        function vrfa = VRFileAccess(filename, maxBytes)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            vrfa.filename = filename;
            if (nargin < 2)
                maxBytes = 2E9;
            end
            vrfa = vrfa.openFile(); %keep file open
            vrfa.maxBytes = maxBytes;
            vrfa = vrfa.locateRollovers(maxBytes/16);
            vrfa = vrfa.initalizeFileChunks;
            ctrfile = strrep(vrfa.filename, '_photontiming.bin', '_str_ctrs.txt');
            if (exist(ctrfile, 'file'))
                 vrfa.hasCounterData = true;
                try
                    vrfa.counterData = parseCounterFile(ctrfile);
                catch me
                     vrfa.hasCounterData = false;
                end
               
            end
            
        end
        
        function rollover = getRollover(vrfa)
            rollover = vrfa.rollover;
        end
        
        function [hasErrors, msg] = consistencyCheck (vrfa)
            hasErrors = false;
            if (any(diff(vrfa.rollover.counter) ~= 1))
                hasErrors = true;
                msg = {'missing rollover counter'};
            else
                msg = {};
            end
            if (any(vrfa.rollover.ncounters(2:end-1) ~= 2^16))
                hasErrors = true;
                msg = [msg; sprintf('%d rollovers missing at least one timer', nnz(vrfa.rollover.ncounters(2:end-1) ~= 2^16))];
            end
            if (~vrfa.hasCounterData)
                msg = [msg; 'streaming counter file (_str_ctr.txt) not present, cannot do further checking'];
                return;
            end
            if (vrfa.counterData.hasErrors)
                hasErrors = true;
                msg = [msg; 'streaming counters reports the following errors in data transmission on fpga or to computer'];
                fn = fieldnames(vrfa.counterData.errors);
                for j = 1:length(fn)
                    if (vrfa.counterData.errors.(fn{j}))
                        msg = [msg; fn{j}]; %#ok<AGROW>
                    end
                end
            end
            ro = vrfa.counterData.rollover;
            ii = interp1(double(ro.cnum), 1:length(ro.cnum), vrfa.rollover.counter, '*nearest', 'extrap');
            
            if (any(ro.cnum(ii) ~= vrfa.rollover.counter))
                msg = [msg; ['counters missing from streaming counter text file:', sprintf(' %d', vrfa.rollover.counter(ro.cnum(ii) ~= vrfa.rollover.counter))]];
            end
            valid = ro.cnum(ii) == vrfa.rollover.counter;
            ii = ii(valid);
            threshold = 32; %allow up to 32 bytes mismatch due to odd synching up
            deltaBytes = ro.bytesInCounter(ii) ~= vrfa.rollover.bytes(valid);        
            if (any(abs(deltaBytes) > threshold))
                hasErrors = true;
                inds = find(abs(deltaBytes) > threshold);
                str = 'bytes in chunk on fpga do not match bytes in chunk on disk (counter / deltaBytes:';
                for j = 1:length(inds)
                    str = [str sprintf(' %d/%g', ro.cnum(ii(inds(j))), delta(inds(j)))]; %#ok<AGROW>
                end
                msg = [msg; str];
            end
            if (~all([vrfa.vrfc.parsed]))
                fprintf('%d / %d vrfcs parsed. run vrfa.parseChunks() first to get most accurate report', nnz([vrfa.vrfc.parsed]), length(vrfa.vrfc));
            end
            if (~any([vrfa.vrfc.parsed]))
                return;
            end
            msg = [msg; sprintf('%d/%d vrfc chunks parsed; comparing these to streaming counters', nnz([vrfa.vrfc.parsed]), length(vrfa.vrfc))];
            ii = interp1(double(ro.cnum), 1:length(ro.cnum), [vrfa.vrfc.counter], '*nearest', 'extrap');            
            valid = ro.cnum(ii) == [vrfa.vrfc.counter] & [vrfa.vrfc.parsed];
            ii = ii(valid);
            
            ntag = [vrfa.vrfc(valid).ntag];
            ctr = [vrfa.vrfc(valid).counter];
            nphoton = [vrfa.vrfc(valid).nphotons];
            
            ntag2 = [ro.tags_b1(ii);ro.tags_b2(ii)];
            nphoton2 = [ro.photons_b1(ii);ro.photons_b2(ii)];
            if (any(ntag ~= ntag2))
                hasErrors = true;
                str = 'tag triggers mismatch (counter,disk/fpga):';
                for k = 1:2
                    if (any(ntag(k,:) ~= ntag2(k,:)))
                        str = [str sprintf(' beam %d:',k)];
                        for j = find(ntag(k,:) ~= ntag2(k,:))
                            str = [str sprintf(' %d, %d/%d', ctr(j), ntag(k,j), ntag2(k,j))];
                        end
                    end
                end
                msg = [msg; str];
            end
            if (any(nphoton ~= nphoton2))
                hasErrors = true;
                str = 'photon counts mismatch (counter,disk/fpga):';
                for k = 1:2
                    if (any(nphoton(k,:) ~= nphoton2(k,:)))
                        str = [str sprintf(' beam %d:',k)];
                        for j = find(nphoton(k,:) ~= nphoton2(k,:))
                            str = [str sprintf(' %d, %d/%d', ctr(j), nphoton(k,j), nphoton2(k,j))];
                        end
                    end
                end
                msg = [msg; str];
            end
            if (~hasErrors)
                msg = [msg; 'number of bytes, photons, and tags are consistent between fpga counters and disk values'];
            end
                   
        end
        
        function [vrfa, cf] = checkFileOpen(vrfa)
           cf = false;
           try
               cf = (ftell(vrfa.fid) < 0);
           catch
               vrfa = vrfa.closeFile;
               cf = true;
           end
        end
        
        function vrfa = parseChunks(vrfa)
            %function vrfa = parseChunks(vrfa)
           [vrfa, cf] = checkFileOpen(vrfa);
           vrfa = vrfa.openFile;
           ts1 = tic;
           vrfa.vrfc.parse(vrfa.fid);
           if (cf && nargout < 1)
               vrfa = vrfa.closeFile();
           end
           elapsed = toc(ts1);
           fprintf('%.1f seconds to parse %.1f megabytes = %.1f MB/sec\n', elapsed, vrfa.rollover.totalBytes/(1024^2), vrfa.rollover.totalBytes/(1024^2)/elapsed);
        end
        
        function [timer, beam, aod, vrfa] = getData(vrfa, beamnum, timerange)
          [vrfa, cf] = checkFileOpen(vrfa);
           vrfa = vrfa.openFile;
           if (isempty(vrfa.galvoDeltaT))
               warning ('galvoDeltaT not set - using 0');
               vrfa.galvoDeltaT = 0;
           end
           gdelta = ones(size(timerange)); gdelta(1) = -1; gdelta = gdelta*vrfa.galvoDeltaT;
           [timer, beam, aod] = vrfa.vrfc.getData(vrfa.fid, beamnum, timerange+gdelta);
           vrfa.manageMemory();
           if (isempty(timer.t))
               return;
           end
           tr = timer.t([1 end]);
           tr(1) = max(tr(1),timerange(1));
           tr(2) = min(tr(2),timerange(2));
           valid = beam.gt >= tr(1) & beam.gt <= tr(2);
           beam.gt = beam.gt(valid);
           beam.gPhase = beam.gPhase(valid);
           valid = beam.rt >= tr(1) & beam.rt <= tr(2);
           beam.rt = beam.rt(valid);
           beam.rPhase = beam.rPhase(valid);
           valid = beam.tagt >= tr(1) & beam.tagt <= tr(2);
           beam.tagt = beam.tagt(valid);
           beam.tagPeriod = beam.tagPeriod(valid);
           
           if (isempty(aod.t))
               aod.t = reshape(timerange,1,[]);
               aod.loc = [0 0;0 0];
           end
           if (cf && nargout < 4)
               vrfa = vrfa.closeFile();
           end
        end
        
        function [pos, phase, time, vrfa] = getPhotonLocations(vrfa, beamnum, color, timerange)
            %[pos, phase, time, vrfa] = getPhotonLocations(vrfa, beamnum, color, timerange)
            %pos(1:2) = x,y pos(3) = piezo [1/5/24 code examination] - for
            %both beams (same value)
           
            if (nargout > 3)
                c;
            else
                [timer, beam, aod] = vrfa.getData(beamnum, timerange);
            end
            if (color == 'a')
                time = [beam.rt, beam.gt];
            else
                time = beam.([color 't']);
            end
            if (numel(timer.t) < 2 || isempty(time))
                pos = zeros([3 0]);
                time = zeros([1 0]);
                phase = zeros([1 0]);
                return;
            end
            if (color == 'a')
                phase = [beam.rPhase, beam.gPhase];
                [time,I] = sort(time);
                phase = phase(I);
            else
                phase = beam.([color 'Phase']);
            end
            aodt = [-1E9 aod.t 1E9];
            aodloc = aod.loc(1:2, [1 1:end end]);
            %have to pad galvo location too, because otherwise galvoDeltaT
            %can throw out of range
            galvot = [-1E9 timer.t+vrfa.galvoDeltaT 1E9];
            galvopos = timer.pos(1:3, [1 1:end end]);
            
            %galvopos = interp1(timer.t + vrfa.galvoDeltaT, timer.pos', time, '*linear')';
            galvopos = interp1(galvot, galvopos', time, '*linear')';
            
            aodpos = interp1(aodt, aodloc', time, '*previous')';
            aodpos(3,:) = 0;
            pos = galvopos + aodpos;
        end
        
          function [pos, phase, tvec, dwelltime, vrfa] = getDwellLocationSubsampled(vrfa, beamnum, timerange, deltaT, isaod)
              %tvec is interpolated between timerange([1 2]) with
              %approximate interval deltaT (updated to make deltaT an even
              %multiple of the period); jitter is added to prevent aliasing
              
              if (numel(timerange) > 2 || numel(deltaT) > 1)
                  error ('revised code wants timerange, deltaT, not tvec')
              end
              if (nargout > 3)
                  [timer, beam, aod, vrfa] = vrfa.getData(beamnum, timerange);
              else
                  [timer, beam, aod] = vrfa.getData(beamnum, timerange);
              end
              if (nargin < 5)
                  isaod = mean(aod.loc(:) ~= 0) > 0.1;
              end
              if (numel(timer.t) < 2 || isempty(beam.tagt))
                  pos = zeros([3 0]);
                  tvec = zeros([1 0]);
                  phase = zeros([1 0]);
                  dwelltime = zeros([1 0]);
                  return;
              end
              
              
              time = beam.tagt;
              period = beam.tagPeriod;
              aodt = [-1E9 aod.t 1e9];
              aodloc = aod.loc(1:2, [1 1:end end]);
              
              
              
              if (isaod) %only sample the first 1/2 of each period
                  expansion = ceil(0.5*mean(period, 'omitnan')/deltaT); %samples per half period
                  deltaT = 0.5*mean(period, 'omitnan')/expansion;
                  poff =  0.25;         
              else
                  expansion = ceil(mean(period, 'omitnan')/deltaT);
                  deltaT = mean(period, 'omitnan')/expansion;
                  poff = 0.5;
              end
              
              time = time + poff*period;
              
              deltaTime = deltaT*repmat((((1-expansion)/2):((expansion-1)/2))', [1 length(time)]);
              jitter = repmat((rand(size(time))-0.5)*deltaT/2, [expansion 1]);
              toff = deltaTime + jitter;
              tvec = reshape(repmat(time,  [expansion 1]) + toff, 1,[]);
              phase = 2*pi*mod(reshape((toff)./repmat(period, [expansion 1]), 1, [])+poff, 1);

              galvopos = interp1(timer.t + vrfa.galvoDeltaT, timer.pos', time, 'linear', 'extrap')';
              galvovel = deriv(galvopos, 0.001)./deriv(time, 0.001); %creates a kernel [-0.5 0 0.5]
              aodpos = interp1(aodt, aodloc', time, 'previous')';
              aodpos(3,:) = 0;
              pp = galvopos + aodpos;
              x = reshape(repmat(pp(1,:), [expansion 1]) + repmat(galvovel(1,:), [expansion 1]).*toff, 1, []);
              y = reshape(repmat(pp(2,:), [expansion 1]) + repmat(galvovel(2,:), [expansion 1]).*toff, 1, []);
              z = reshape(repmat(pp(3,:), [expansion 1]) + repmat(galvovel(3,:), [expansion 1]).*toff, 1, []);
              pos = [x;y;z];
              dwelltime = ones(size(tvec))*deltaT;
          end
         
         function [pos, time, period, vrfa] = getDwellLocations(vrfa, beamnum, timerange, isaod)
         %function [pos, time, vrfa] = getDwellLocations(vrfa, beamnum, timerange, isaod)
         %
         %reports the central location of each tag sweep. If isaod is true,
         %only reports the forward (phase = pi/4) portion of the sweep
         %otherwise, reports position at phase = pi/4, 3pi/4
            if (nargout > 3)
                [timer, beam, aod, vrfa] = vrfa.getData(beamnum, timerange);
            else
                [timer, beam, aod] = vrfa.getData(beamnum, timerange);
            end
            if (nargin < 4)
                isaod = mean(aod.loc(:) ~= 0) > 0.1;
            end
            if (numel(timer.t) < 2 || isempty(beam.tagt))
                pos = zeros([3 0]);
                time = zeros([1 0]);
                period = zeros([1 0]);
                return;
            end
            
            
            time = beam.tagt;
            period = beam.tagPeriod;
            aodt = [-1E9 aod.t 1e9];
            aodloc = aod.loc(1:2, [1 1:end end]);
            
            if (isaod)
                time = time + 0.25*period;
            else
                t1 = time + 0.25*period;
                t2 = time + 0.75*period;
                time = zeros([1 length(time)*2]);
                time(1:2:end) = t1;
                time(2:2:end) = t2;
                pp = period;
                period = reshape([pp;pp], [], 2*length(pp));
            end
            
            galvopos = interp1(timer.t + vrfa.galvoDeltaT, timer.pos', time, 'linear', 'extrap')';
            aodpos = interp1(aodt, aodloc', time, 'previous', 'extrap')';
            aodpos(3,:) = 0;
            pos = galvopos + aodpos;
            
         end
        
         function [phase, deltaT] = getPhaseFromTime(vrfa, beamnum, tvec)
             %phase = tag lens phase at tvec
             %deltaT = time since last tag trigger at tvec
            [~,beam] = vrfa.getData(beamnum, tvec([1 end]));
            tdat = interp1(beam.tagt, [beam.tagt' beam.tagPeriod'], tvec, 'previous')';
            deltaT = tvec - tdat(1,:);
            phase = 2*pi*deltaT./tdat(2,:);
        end
        
        
        function bytes = memoryUsed(vrfa)
            bytes = sum(vrfa.vrfc.memoryUsed);
        end
        
        function manageMemory(vrfa, maxBytes)
            if (nargin < 2)
                maxBytes = vrfa.maxBytes;
            end
            vrfa.vrfc.manageMemory(maxBytes);
        end
        
        function unload(vrfa)
            vrfa.vrfc.unload();
        end
        
        function tr = getTimeRange(vrfa)
            tr = vrfa.vrfc(1:end-1).getTimeRange();
        end
        
        function tx = getMemoryLimitedTimeSteps(vrfa, tr, chunksize)
            t = vrfa.rollover.counter * VRFileChunk.ticksPerCounterIncrement/VRFileChunk.FPGA_CLOCK;
            m = cumsum(vrfa.rollover.bytes);
            ms = interp1(t, m, tr(1), 'linear', min(m));
            me = interp1(t,m,tr(end), 'linear', max(m));
            nchunks = ceil((me-ms)/chunksize);
            mx = linspace(ms,me,nchunks + 1);
            tx = interp1(m, t, mx);
        end
        
         function vrfa = openFile(vrfa)
             %do nothing if file is already open
            try
                if (ftell(vrfa.fid) < 0)
                    vrfa = vrfa.closeFile;
                end
            catch
                vrfa.fid = -1;
            end
            if (vrfa.fid <= 0 || ftell(vrfa.fid) < 0)
                vrfa.fid = fopen(vrfa.filename, 'r');              
            end
            if (vrfa.fid <= 0)
                error (['failed to open file ' vrfa.filename ' for reading']);
            end
        end
        function vrfa = closeFile(vrfa)
            if (vrfa.fid > 0)
                try
                    fclose(vrfa.fid);
                catch
                end
            end
            vrfa.fid = -1;
        end
    end
    
    methods (Access = protected)
        function vrfa = locateRollovers(vrfa, chunkSize)
            vrfa = vrfa.openFile();
            fseek(vrfa.fid, 0,1); 
            totalBytes = ftell(vrfa.fid); 
            fseek(vrfa.fid,0,-1);
            vrfa.rollover.totalBytes = totalBytes;
            
            if (nargin < 2)
                chunkSize = 2.5E8; %500 MB chunks
            end
            floc = 0;
            iter = 1;
            counterval = [];
            counterloc = [];
            while (true)
                fseek(vrfa.fid, floc, -1);
                nread = min((totalBytes-floc)/2, chunkSize);
                chunk = uint16(fread(vrfa.fid, [1 nread], 'uint16'));
                [timer, lastInd] = findTimers(chunk);
                ro = timer.rollover;
                counterval = [counterval timer.counter]; %#ok<AGROW>
                counterloc = [counterloc (floc+timer.inds - 1)*2];  %#ok<AGROW>
                fullCounter(iter).counter = ro.counters(ro.ntimers == 2^16);%#ok<AGROW>
                fullCounter(iter).location = floc + (ro.inds(ro.ntimers == 2^16)-1)*2;%#ok<AGROW>
                [partialCounter(iter).counter,I] = setdiff(ro.allCounters, fullCounter(iter).counter);%#ok<AGROW>
                partialCounter(iter).location = floc + (ro.allInds(I)-1)*2;%#ok<AGROW>
                if (nread < chunkSize) %read to end of file
                    break;
                end
                floc = floc + (lastInd-1)*2;
                iter = iter+1;
            end

            
            fullcounters = [fullCounter.counter];
            fullLoc = [fullCounter.location];
            
            [fullcounters, IA, IC] = unique(fullcounters);
            if (length(fullcounters) < IC)
                %check that duplicate counters have duplicate locations
                for j = 1:length(fullcounters)
                    if (any(fullLoc(IC == j) ~= fullLoc(IA(j))))
                        error ('did not parse counters correctly');
                    end
                end
            end
            
            partcounters = [partialCounter.counter];
            pL = [partialCounter.location];
            
            partLoc = accumarray(partcounters', pL, [max(partcounters) 1], @min)';
            
            partcounters = setdiff(partcounters, fullcounters);
            partLoc = partLoc(partcounters);
            
            ac = [partcounters, fullcounters];
            al = [partLoc, fullLoc];
            [vrfa.rollover.counter,I] = sort(ac, 'ascend');
            vrfa.rollover.loc = al(I);
            vrfa.rollover.bytes = diff([vrfa.rollover.loc vrfa.rollover.totalBytes]);
            vrfa.rollover.rawtimer = timer;
            [~,I] = unique(counterloc);
            counterval = counterval(I);
            vrfa.rollover.cx = min(counterval):max(counterval);
            vrfa.rollover.ncounters = hist(counterval, vrfa.rollover.cx);
            while (true)
                I = find(diff(vrfa.rollover.counter) < 0,1,'first');
                if (isempty(I))
                    break;
                end
                warning ('main counter rollover detected; this may affect other data synching if rollovers are not compensated there as well');
                if (vrfa.rollover.counter(I) ~= 4095)
                    warning ('reset to 0 at non-rollover point - expect counter to go to 4095 and reset');
                end
                vrfa.rollover.counter((I+1):end) = vrfa.rollover.counter((I+1):end) + 4096;
            end
            
  
        end
        
        function vrfa = initalizeFileChunks(vrfa)
           
            if (isempty(vrfa.rollover) || ~isstruct(vrfa.rollover) || ~isfield(vrfa.rollover, 'counter'))
                return;
            end
            loc = [vrfa.rollover.loc vrfa.rollover.totalBytes-16];
            vrfa.vrfc = repmat(VRFileChunk([0 0], [0 0], 0), [1 length(vrfa.rollover.counter)]);
            for j = 1:length(vrfa.rollover.counter) 
                vrfa.vrfc(j) = VRFileChunk([loc(j) loc(j+1)+16], VRFileChunk.ticksPerCounterIncrement/VRFileChunk.FPGA_CLOCK * (vrfa.rollover.counter(j) + [0 1]), 4096*floor(vrfa.rollover.counter(j)/4096));
            end
        end
        
        
        
       
        
    end
    
end

