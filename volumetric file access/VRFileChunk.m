classdef VRFileChunk < handle
    % VRFileChunk
    %   accesses portion of a file that contains all counters for one
    %   rollover
    
    properties (Constant)
        ticksPerCounterIncrement = 2^28;
        FPGA_CLOCK = 1.6E8;
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        
        timerange
        fileLoc
        counterOffset = 0;
        bytesUsed = 0;
    end
    properties (Transient = true, Access = protected)
        timer
        beam
        aod
        inmemory = false;
        timeLastAccessed = Inf;
    end
    properties (GetAccess = public, SetAccess = protected)
        parsed = false; %whether this section of the data has been loaded from file at least once
        counter
        unprocessedBytes
        ntag
        ncounter
        u_deltacounter %mean change in ticks between counters
        s_deltacounter %std change in ticks between counters
        bytesOnDisk
        nphotons
        parseIssues
    end
    
    methods
        function vrfc = VRFileChunk(fileLoc, timerange,counterOffset)
            %vrfc = VRFileChunk(fileLoc, timerange,tickOffset)
            %   
            vrfc.fileLoc = fileLoc;
            vrfc.timerange = timerange;
            vrfc.counterOffset = counterOffset;
        end
        
        function parse(vrfc, fid)
            %function parse(vrfc, fid)
            %loads and immediately unloads vrfc
            %use if you want to iterate through without accessing data for
            %consistency checks
            ts1 = tic;
            for j = 1:length(vrfc)
                vrfc(j).load(fid);
                vrfc(j).unload();
                if (toc(ts1) > 30)
                    fprintf('%d / %d chunks parsed\n', j, length(vrfc));
                    ts1 = tic;
                end
            end
        end
        
        function load(vrfc,fid)
        %loads data from file into memory
            if (length(vrfc)> 1)
                error ('require specifying a single vrfc to load at a time to prevent accidental memory overfill');
            end
            if (vrfc.inmemory)
                return;
            end
            fseek(fid, vrfc.fileLoc(1),-1);
            chunk = uint16(fread(fid, [1 diff(vrfc.fileLoc)/2], 'uint16'));    
            [vrfc.timer, vrfc.beam, vrfc.aod, lastInd, vrfc.parseIssues] = parseTemporalChunkMinimal(chunk);
            vrfc.unprocessedBytes = vrfc.bytesOnDisk - lastInd*2;
            vrfc.postLoadProcess();
            vrfc.inmemory = true;
            vrfc.timeLastAccessed = tic();
            vrfc.calcMemoryUsage();
        end
        function unload(vrfc)
            if (length(vrfc) > 1)
                for j = 1:length(vrfc)
                    vrfc(j).unload;
                end
                return;
            end
            vrfc.timer = [];
            vrfc.beam = [];
            vrfc.aod = [];
            vrfc.inmemory = false;
            vrfc.timeLastAccessed = Inf;
            vrfc.calcMemoryUsage();
        end
       
        function tr = getTimeRange(vrfc)
            tr = [vrfc(1).timerange(1) vrfc(end).timerange(end)];
        end

        function trs = getTimeRangeVectors(vrfc)
            trs = cat(1,vrfc.timerange);
        end
                
        
        function bytes = memoryUsed(vrfc)
            bytes = sum([vrfc.bytesUsed]);
            
            
        end
        
        function [timer, beam, aod] = getData(vrfc, fid, beamnum, timerange, bufferTime)
            %bufferTime is the time to expand time range on AOD and timer
            %(but not beam) to allow for later interpolation
            %default = 30 microseconds (capture at least 1 additional timer
            %on each end)
            
            if (nargin < 5)
                bufferTime = 3E-5;
            end
            [timer,beam,aod] = vrfc.emptyResponse();
            if (length(vrfc) > 1)
                for j = 1:length(vrfc)
                    [tt(j), bb(j), aa(j)] = vrfc(j).getData(fid, beamnum, timerange, bufferTime); %#ok<AGROW>
                end
                fn = fieldnames(timer);
                for j = 1:length(fn)
                    timer.(fn{j}) = cat(2, tt.(fn{j}));
                end
                fn = fieldnames(beam);
                for j = 1:length(fn)
                    beam.(fn{j}) = cat(2, bb.(fn{j}));
                end
                 fn = fieldnames(aod);
                for j = 1:length(fn)
                    aod.(fn{j}) = cat(2, aa.(fn{j}));
                end
                return;
            end
            if (~vrfc.checkTimeRange(timerange, fid))
                return;
            end
            vrfc.timeLastAccessed = tic();
            
            bufftr = timerange([1 end]) + [-bufferTime bufferTime];
            valid = vrfc.timer.tickCount >= bufftr(1)*vrfc.FPGA_CLOCK & vrfc.timer.tickCount <= bufftr(end)*vrfc.FPGA_CLOCK;
            timer.t = vrfc.timer.tickCount(valid)/vrfc.FPGA_CLOCK;
            timer.pos = vrfc.timer.pos6((1:3) + (beamnum-1)*3, valid);
            
            valid = vrfc.aod(beamnum).tickCount >= bufftr(1)*vrfc.FPGA_CLOCK & vrfc.aod(beamnum).tickCount <= bufftr(end)*vrfc.FPGA_CLOCK;
            aod.t = vrfc.aod(beamnum).tickCount(valid)/vrfc.FPGA_CLOCK;
            aod.loc = vrfc.aod(beamnum).loc(:,valid);
            
            valid = vrfc.beam(beamnum).gelapsed >= timerange(1)*vrfc.FPGA_CLOCK & vrfc.beam(beamnum).gelapsed <= timerange(end)*vrfc.FPGA_CLOCK;
            beam.gt = vrfc.beam(beamnum).gelapsed(valid)/vrfc.FPGA_CLOCK;
            beam.gPhase = vrfc.beam(beamnum).gPhase(valid);
            
            valid = vrfc.beam(beamnum).gelapsed >= timerange(1)*vrfc.FPGA_CLOCK & vrfc.beam(beamnum).gelapsed <= timerange(end)*vrfc.FPGA_CLOCK;
            beam.gt = vrfc.beam(beamnum).gelapsed(valid)/vrfc.FPGA_CLOCK;
            beam.gPhase = vrfc.beam(beamnum).gPhase(valid);
            
            valid = vrfc.beam(beamnum).relapsed >= timerange(1)*vrfc.FPGA_CLOCK & vrfc.beam(beamnum).relapsed <= timerange(end)*vrfc.FPGA_CLOCK;
            beam.rt = vrfc.beam(beamnum).relapsed(valid)/vrfc.FPGA_CLOCK;
            beam.rPhase = vrfc.beam(beamnum).rPhase(valid);
            
            valid = vrfc.beam(beamnum).elapsed >= timerange(1)*vrfc.FPGA_CLOCK & vrfc.beam(beamnum).elapsed <= timerange(end)*vrfc.FPGA_CLOCK;
            beam.tagt = vrfc.beam(beamnum).elapsed(valid)/vrfc.FPGA_CLOCK;
            beam.tagPeriod = vrfc.beam(beamnum).period(valid)/vrfc.FPGA_CLOCK;
        end
        
        
        
        function manageMemory(vrfc, maxBytes)
            %given vector of vrfc, unloads those who have not been accessed
            %in the longest time to keep total memory usage under maxBytes
            bytes = [vrfc.bytesUsed];
            if (sum(bytes) <= maxBytes)
                return;
            end
            [~,I] = sort([vrfc.timeLastAccessed], 'descend');
            ind = find(cumsum(bytes(I)) > maxBytes, 1, 'first');
            vrfc(I(ind:end)).unload;
        end
       
    end
    
    methods (Access = protected)
        function postLoadProcess(vrfc)
            vrfc.timer.counter = vrfc.timer.counter + vrfc.counterOffset;
            vrfc.timer.validCounters = vrfc.timer.validCounters + vrfc.counterOffset;
            deltaT = vrfc.counterOffset*vrfc.ticksPerCounterIncrement;
            vrfc.timer.tickCount = vrfc.timer.tickCount + deltaT;
            for j = 1:2
                if (isempty(vrfc.aod(j).inds))
                    vrfc.aod(j).tickCount = [min(vrfc.timer.tickCount)-1 max(vrfc.timer.tickCount)+1];
                    vrfc.aod(j).loc = [0 0;0 0];
                    vrfc.aod(j).timermismatch = [false false];
                end
                vrfc.aod(j).tickCount = vrfc.aod(j).tickCount + deltaT;
                vrfc.beam(j).elapsed = vrfc.beam(j).elapsed + deltaT;
                vrfc.beam(j).gelapsed = vrfc.beam(j).gelapsed + deltaT;
                vrfc.beam(j).relapsed = vrfc.beam(j).relapsed + deltaT;
            end
            vrfc.timerange = vrfc.timer.tickCount([1 end])/vrfc.FPGA_CLOCK;   
            
            vrfc.counter = mode(vrfc.timer.counter);
            valid = vrfc.timer.counter == vrfc.counter;
            if (nnz(~valid) > 1)
                warning('vrfc is supposed to contain all timers for one counter + one timer for next counter only');
            end
            fn = {'inds', 'counter', 'validTimer', 'tickCount', 'pos6'};
            for j = 1:length(fn)
                vrfc.timer.(fn{j}) = vrfc.timer.(fn{j})(:,valid);
            end
            
            vrfc.pareToMinimum();
            
            vrfc.ncounter = length(vrfc.timer.tickCount);
            vrfc.u_deltacounter = mean(diff(vrfc.timer.tickCount));
            vrfc.s_deltacounter = std(diff(vrfc.timer.tickCount));
            vrfc.ntag = [length(vrfc.beam(1).elapsed);length(vrfc.beam(2).elapsed)];
            vrfc.nphotons = [length(vrfc.beam(1).relapsed)+length(vrfc.beam(1).gelapsed);length(vrfc.beam(2).relapsed)+length(vrfc.beam(2).gelapsed)];
            vrfc.bytesOnDisk = diff(vrfc.fileLoc);
            vrfc.parsed = true;
        end
        
        function pareToMinimum(vrfc)
            vrfc.timer = rmfield(vrfc.timer, setdiff(fieldnames(vrfc.timer), {'tickCount', 'pos6'}));
            vrfc.aod = rmfield(vrfc.aod, setdiff(fieldnames(vrfc.aod), {'tickCount', 'loc'}));
            vrfc.beam = rmfield(vrfc.beam, setdiff(fieldnames(vrfc.beam), {'gelapsed', 'gPhase', 'relapsed', 'rPhase', 'elapsed', 'period'}));
                       
        end
        
        
        function inrange = checkTimeRange(vrfc, times, fid)
            inrange = vrfc.timerange(1) < max(times) && vrfc.timerange(2) >= min(times);
            if (inrange && ~vrfc.inmemory && nargin > 2 && fid > 0)
                vrfc.load(fid);
            end
        end
        
        function calcMemoryUsage(vrfc)
            props = union(properties(vrfc),{'timer', 'beam', 'aod', 'counter','timerange','fileLoc','inmemory','timeLastAccessed','counterOffset'});
            vrfc.bytesUsed = 0;
            
            for ii=1:length(props)
                currentProperty = vrfc.(char(props(ii))); %#ok<NASGU>
                s = whos('currentProperty');
                 vrfc.bytesUsed =  vrfc.bytesUsed + s.bytes;
            end
            
        end
        
        function vrfc = saveobj(vrfc)
            vrfc.unload();
        end
    end
    
    methods (Static, Access = protected)
        function [timer, beam, aod] = emptyResponse()
            timer.t = zeros([1 0]);
            timer.pos = zeros([3 0]);
            beam.gt = zeros([1 0]);
            beam.gPhase = zeros([1 0]);
            beam.rt = zeros([1 0]);
            beam.rPhase = zeros([1 0]);
            beam.tagt = zeros([1 0]);
            beam.tagPeriod = zeros([1 0]);
            aod.loc = zeros([2 0]);
            aod.t = zeros([1 0]);
        end        
    end
    
end

