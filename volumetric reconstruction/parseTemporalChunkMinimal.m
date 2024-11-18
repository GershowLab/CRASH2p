function [timer, beam, aod, lastInd, parseIssue] = parseTemporalChunkMinimal(chunk)
% [timer, beam, aod, diagnostics, lastInd] = parseTemporalChunk(chunk)
% parses a chunk of data loaded from a file - minimal processing; does not
% do any extra extrapolation or calculation
%
% new encoding (beginning 2/25/19; AOD beginning 1/8/2020) is
% 4 bit code ; 12 bit data
% 1 = ch0 photon (green), beam 1 ; ticks since last trigger 1
% 2 = ch1 photon (red), beam 1; ticks since last trigger 1
% 3 = ch0+ch1, beam 1; ticks since last trigger 1
% 7 = tag trigger 1; ticks since last trigger 1
% 9,10,11,15 are same as 1,2,3,7 except for beam 2, tag trigger 2
% 8 = TIMER; # of 160 MHz ticks on counter / 2^28 (increments every 1.6777
% seconds)
% following code 8, this sequence of 7 16 bit integers is output
% xsum, x1, y1, z1, x2, y2, z2, where xsum makes the XOR of all 6 integers
% 0
% 12 = AOD update beam 1 = TIMER/16 (increments every 16 ticks/ 0.1 ns) 
% 13 = AOD update beam 2 = TIMER/16 (increments every 16 ticks/ 0.1 ns)
% following 12,13 this sequence of 3 16 bit integers is output
% xsum, x, y where xsum is the xor of x,y
% all locations are specified in (best estimated) microns / 2^6. IE a
% reading of 64 is 1 micron
%
% note that although this is set up to parse only a chunk of the file,
% the parsed data takes significantly (~8.4x) more space (because it is a double,
% and because more data is stored) than the original file. If the whole parsed data will fit in memory,
% there is no reason not to load the whole file 
%
% if the parsed data will not fit in memory, then analysis will have to be
% carried out piecemail
 
   
    parseIssue = uint32(0);
    
    numShortsRead = length(chunk);
    code = bitshift(bitand(chunk,uint16(61440)),-12); %get top 4 bits
    data = double(bitand(chunk,uint16(4095))); %bottom 12 bits

    
    
    timerinds = find(code == 8);
    timerinds = timerinds(timerinds <= numShortsRead - 7); %in case no more to be read from chunk
    locinds = repmat(timerinds, [7 1]) + repmat((1:7)', [1 size(timerinds,2)]);
    
    aodinds = find(code == 12 | code == 13);
    aodinds = aodinds(aodinds <= numShortsRead - 3); %in case no more to be read from chunk
    aodlocinds = repmat(aodinds, [3 1]) + repmat((1:3)', [1 size(aodinds,2)]);
    
    timerarr = chunk(locinds);
    aodarr = reshape(chunk(aodlocinds),size(aodlocinds));
    
    %possible ways you can get code 8, 12, 13
    %due to overflow, due to aod trigger, in location for either of these
   
    %first let's check xsums
    xs = timerarr(1,:);
    for q = 2:7
        xs = bitxor(xs, timerarr(q,:));
    end
    locvalid = (xs == 0);
    
     %first let's check xsums
    xs = aodarr(1,:);
    for q = 2:3
        xs = bitxor(xs, aodarr(q,:));
    end
    aodvalid = (xs == 0);
    
    badtimerinds = timerinds(~locvalid);
    timerinds = timerinds(locvalid);
    badaodinds = aodinds(~aodvalid);
    aodinds = aodinds(aodvalid);
    
    %now extract the putatively valid location indices
    locinds = repmat(timerinds, [7 1]) + repmat((1:7)', [1 size(timerinds,2)]);
    if (isempty(aodinds))
        aodlocinds = zeros([3 0]);
    else
        aodlocinds = repmat(aodinds, [3 1]) + repmat((1:3)', [1 size(aodinds,2)]);
    end
    allloc = [locinds(:); aodlocinds(:)];
    %first, let's verify that all bad checksums are in locations somewhere
    
    badaodinds = setdiff(badaodinds, allloc);
    badtimerinds = setdiff(badtimerinds, allloc);
    badinds = union(badaodinds, badtimerinds);
    if ~isempty(badinds)
        parseIssue = bitor(parseIssue, StreamingDataParseError.invalid_checksums);
        warning ([num2str(length(badinds)) ' invalid checksums discovered']);
        fprintf(2, '%d bad code 8, %d bad code 12, %d bad code 13\n', nnz(code(badinds) == 8), nnz(code(badinds) == 12), nnz(code(badinds) == 13));
    end
    
    %now let's check that there aren't any overlapping hits

    aodintimer = intersect(aodinds, locinds);
    if (~isempty(aodintimer))
        %a possible reason there might be an aod hit inside a timer hit is
        %the pattern xs, x,y,z, 0, 0, 0 would trigger an AOD hit if the top
        %4 bits of z are 12 or 13
        %so let's take care of this case and not raise any alarms
        allzero = false(size(aodintimer));
        for j = 1:length(aodintimer)
            allzero(j) = all(chunk(aodintimer(j) + (1:3)) == 0);
        end
        
        aodinds = setdiff(aodinds, aodintimer);
        aodlocinds = repmat(aodinds, [3 1]) + repmat((1:3)', [1 size(aodinds,2)]);
        if (any(~allzero))
            parseIssue = bitor(parseIssue, StreamingDataParseError.aod_in_timer);
            warning(['found ' num2str(nnz(~allzero)) ' aod valid checksums inside timer location block ']);
        end
    end
    
    timerinaod = intersect(timerinds, aodlocinds);
    if (~isempty(timerinaod))
        %this should always be an error
        parseIssue = bitor(parseIssue, StreamingDataParseError.timer_in_aod);
        timerinds = setdiff(timerinds, timerinaod);
        locinds = repmat(timerinds, [7 1]) + repmat((1:7)', [1 size(timerinds,2)]);
      %  if (any(~allzero))
            warning(['found ' num2str(length(timerinaod)) ' timer valid checksums inside aod location block ']);
       % end
    end
        
    timerarr = chunk(locinds);
    aodarr = chunk(aodlocinds);
    
    for j = 1:2
        av = code(aodinds) == 11 + j;        
        aod(j).inds = aodinds(av);
        aod(j).ticks = double(data( aod(j).inds))*16;
        ap = aodarr(2:3,av);
        aod(j).loc = double(reshape(typecast(ap(:), 'int16'), size(ap)))/64;
    end
 
    timer.badInds = badtimerinds;
    timer.inds = timerinds;
    timer.counter = data(timerinds);
    while (true)
        I = find(diff(timer.counter) < 0,1,'first');
        if (isempty(I))
            break;
        end
        warning ('main counter rollover detected; this may affect other data synching if rollovers are not compensated there as well');
        if (timer.counter(I) ~= 4095)
            parseIssue = bitor(parseIssue, StreamingDataParseError.unexpected_rollover);
            warning ('reset to 0 at non-rollover point - expect counter to go to 4095 and reset');
        end
        timer.counter((I+1):end) = timer.counter((I+1):end) + 4096;
    end
  
    tcedges = min(timer.counter):(max(timer.counter) + 1);
   
    ncounters = histcounts(timer.counter, tcedges); tcedges = tcedges(1:end-1);
    problem = ncounters ~= 2^16;% & tcedges > min(timer.counter) & tcedges < max(timer.counter); %we will just discard the end counters anyway
   % problem([1 end]) = false;
    timer.validCounters = tcedges(~problem);
    rolloverInds = find([true diff(timer.counter)>0]);
    rolloverCounters = timer.counter(rolloverInds);
    
    timer.rollover.allCounters = rolloverCounters;
    timer.rollover.allInds = timer.inds(rolloverInds);
    
    % [C,IA,IB] = intersect(A,B)
    [timer.rollover.counters,I] = intersect(rolloverCounters, timer.validCounters);
    timer.rollover.whichtimer = rolloverInds(I);
    timer.rollover.inds = timer.inds(timer.rollover.whichtimer);
    for j = 1:length(timer.rollover.counters)
        timer.rollover.ntimers(j) = nnz(timer.counter == timer.rollover.counters(j));
    end
    timer.validTimer = ismember(timer.counter, timer.rollover.counters);
    if (numel(timer.rollover.whichtimer) > 1)
        timer.tickCount = round(interp1(timer.rollover.whichtimer, 2^28*timer.rollover.counters, 1:length(timer.inds), '*linear', 'extrap'));
    else
        if (isempty(timer.rollover.whichtimer))
            if (isempty(timer.rollover.allCounters))
                error ('cannot determine timing without a counter increment (need to load at least 1.7 sec of data');
            end
            refcounter = timer.rollover.allCounters(end);
            refind = rolloverInds(end);
        else
            refcounter = timer.rollover.counters;
            refind = timer.rollover.whichtimer;
        end
        timer.tickCount = ((1:length(timer.inds)) - refind)*4096 + 2^28*refcounter;
    end
    if (any(problem(2:end-1)))
        warning('missing timers');
        for k = find(problem(2:end-1))
            disp(['counter = ', num2str(tcedges(k+1)), '; num timers = ', num2str(ncounters(k+1))]);
        end
    end
    
    %get absolute timer for aods
    %aod time rolls over every 2^16 ticks
    %each counter should be 2^12 ticks
    %so bits 12-15 should be the same for both timers if everything is
    %correct
    %depending on when increment happens and order of AOD/timer update,
    %could be that only next or previous is right - so check against both &
    %use that
    
    % timer put into timer FIFO 1 cycle after rollover
    % aod put into AOD FIFO on change detection
    % timer dequeues with priority over AOD
    % possible sequences - tp = timer <= AOD timem, tn = timer after AOD
    % time - tp has correct bits
    % tp queue/dequeue AOD queue/dequeue tn queue/dequeue [correct timer is
    % 1 before]
    % tp queue/dequeue (AOD queue tn queue) tn dequeue AOD deqeue [correct
    % timer is 2 before]
    % (AOD queue tp queue with 1 count delay) AOD dequeue tp dequeue ([correct timer is 1 after]) - this is unlikeley but might as well check 
    
    %timer ticks are found using double interpolation,
    %so delta > 1 when timer > 2^51 ticks, or about 167 days -- we can take
    %that risk
    for j = 1:2
        if (isempty(aod(j).inds))
            continue;
        end
        tiprev = interp1(timer.inds, 1:length(timer.inds), aod(j).inds, '*prev', 'extrap'); %index of the timer before the aod in the file
        tiprev(~isfinite(tiprev)) = 1;
        timerPrev = timer.tickCount(tiprev);
        timerPrev2 = timer.tickCount(max(1, tiprev-1));
        timerNext = timer.tickCount(min(tiprev+1, length(timer.tickCount)));
%         timerPrev = interp1(timer.inds, timer.tickCount, aod(j).inds, '*prev', 'extrap');
%         timerNext = interp1(timer.inds, timer.tickCount,  aod(j).inds, '*next', 'extrap');
        prevGood =  bitand(uint64(timerPrev), 61440) ==  bitand(uint64(aod(j).ticks), 61440);
        nextGood =  bitand(uint64(timerNext), 61440) ==  bitand(uint64(aod(j).ticks), 61440);
        prev2Good =  bitand(uint64(timerPrev2), 61440) ==  bitand(uint64(aod(j).ticks), 61440);
        
        oneGood = prevGood | nextGood | prev2Good;
        if (any(~oneGood))
            parseIssue = bitor(parseIssue, StreamingDataParseError.aods_do_not_match);
            warning ([num2str(nnz(~oneGood)) ' aods do not match their forward or backward timers']);
        end
        aod(j).tickCount = double(bitand(uint64(timerPrev),  bitcmp(65535, 'uint64'))) + aod(j).ticks;
        aod(j).tickCount(nextGood) = double(bitand(uint64(timerNext(nextGood)),  bitcmp(65535, 'uint64'))) + aod(j).ticks(nextGood);
        aod(j).tickCount(prev2Good) = double(bitand(uint64(timerNext(prev2Good)),  bitcmp(65535, 'uint64'))) + aod(j).ticks(prev2Good);
        ii = find(oneGood);
         aod(j).timermismatch = ~oneGood;
        if (numel(ii) < 2)
            warning ('aod and timer really do not add up');
            continue;
        end
        aod(j).tickCount(~oneGood) = interp1(ii, aod(j).tickCount(ii), find(~oneGood), 'linear');
       
        
    end
    

    pos6 = timerarr(2:7,:);
    timer.pos6 = double(reshape(typecast(pos6(:), 'int16'), size(pos6)))/64;
   
    if (isempty(aodinds))
        lastInd = max(timerinds) + 7;
    else
        lastInd = max(max(aodinds + 3), max(timerinds)+7);
    end
    chunk = chunk(1:lastInd);
    code = code(1:lastInd);
    data = data(1:lastInd);
    valid = true(size(chunk));
    valid(timerinds) = false; valid(aodinds) = false; valid(aodlocinds) = false; valid (locinds) = false;
    

    codes = 0:16;
    ncodes = histcounts(code(valid), codes); codes = 0:15;
    unused = [0 4 5 6 8 12 13 14];
    if (any(ncodes(unused + 1)))
        warning ('bad code detected');
        disp(['bad codes: ', num2str(unused)]);
        bar(axes(figure(13)), codes, ncodes);
        for k = find(ncodes)
            disp(['code = ' num2str(codes(k)) '; n = ' num2str(ncodes(k))]);
        end
    end
    delta = [0 8];
    for j = 1:2
        beam(j).tagInd = double(find(code == (7+delta(j)) & valid)); %#ok<*AGROW>        
        beam(j).period = data(beam(j).tagInd);
        if (isempty(beam(j).tagInd))
            continue;
        end
        if (std(beam(j).period,'omitnan') > 0.2 * mean(beam(j).period,'omitnan'))
            %new timing scheme - record tick count mod 4 
            %so tick count * 4 is a 14 bit number; the highest 2 bits
            %should match the same 2 bits of the timer tick count before
            %sometimes it might match the 2 bits of the timer 2 before,
            %because the fpga software prioritizes timers and aods to keep
            %the bytes together; thus it will hold a tag for a timer but
            %never the other way around
            %
            %however, if we are missing timers, then there is some
            %uncertainty in which timer is correct, so we should just take
            %the closest that matches
            %
            maskForAgreement = bitset(bitset(uint64(0), 14),13);
            maskForAddition =  bitcmp(16383, 'uint64');
            [tickCount, offset, oneGood, validTimer] = alignTimes (timer, beam(j).tagInd, 4*beam(j).period, maskForAgreement, maskForAddition);
            
            beam(j).elapsed = double(tickCount);
            problem = (offset ~= 0 & offset ~= -1 & validTimer) | ~oneGood; %problem if the timer is valid (i.e. there are the correct number of timers with the same counter) and it isn't the first or second before the tag trigger
            if (any(problem))
                 parseIssue = bitor(parseIssue, StreamingDataParseError.tags_do_not_match);
                warning (['counter: ' num2str(ceil(max(timerPrev/2^28))) ', beam ' num2str(j) ': ' num2str(nnz(problem)) ' tag triggers do not match either of the two previous timers']);
            end
            ii = find(oneGood);
            beam(j).timermismatch = problem;
            if (numel(ii) < 2 && length(ii) > 2)
                %the way alignTimes is written, this should never happen,
                %because we test against a huge number of forward and
                %backward timers
                warning ('huge parse problem - tag timers do not match counter timers at all');                
            else
                beam(j).elapsed(~oneGood) = interp1(ii, beam(j).elapsed(ii), find(~oneGood), 'linear');
            end
            beam(j).period = [diff(beam(j).elapsed(1:2)) diff(beam(j).elapsed)];
%             
%             ticks = 4*beam(j).period;
%             tiprev = interp1(timer.inds, 1:length(timer.inds), beam(j).tagInd, '*prev', 'extrap'); %index before the timer in the file
%             tiprev(~isfinite(tiprev)) = 1;
%             timerPrev = timer.tickCount(tiprev);
%             timerPrev2 = timer.tickCount(max(1, tiprev-1));
% %             timerPrev = interp1(timer.inds, timer.tickCount, beam(j).tagInd, '*prev', 'extrap');
% %             timerPrev2 = interp1(timer.inds, timer.tickCount,  beam(j).tagInd, '*prev', 'extrap');
%             mask = bitset(bitset(uint64(0), 14),13);
%             prevGood =  bitand(uint64(timerPrev), mask) ==  bitand(uint64(ticks), mask);
%             prev2good =  bitand(uint64(timerPrev2), mask) ==  bitand(uint64(ticks), mask);
%             oneGood = prevGood | prev2good;
%             if (any(~oneGood))
%                  parseIssue = bitor(parseIssue, StreamingDataParseError.tags_do_not_match);
%                 warning (['counter: ' num2str(ceil(max(timerPrev/2^28))) ', beam ' num2str(j) ': ' num2str(nnz(~oneGood)) ' tag triggers do not match either of the two previous timers']);
%             end
% %             bothgood = prevGood & nextGood;
% %             if (any(timerNext(bothgood) ~= timerPrev(bothgood)))
% %                 warning ('wtf');
% %             end
%             beam(j).elapsed = double(bitand(uint64(timerPrev),  bitcmp(16383, 'uint64'))) + ticks;
%             beam(j).elapsed(prev2good) = double(bitand(uint64(timerPrev2(prev2good)),  bitcmp(16383, 'uint64'))) + ticks(prev2good);
%             ii = find(oneGood);
%             beam(j).timermismatch = ~oneGood;
%             if (numel(ii) < 2 && length(ii) > 2)
%                 warning ('huge parse problem - tag timers do not match counter timers at all');
%                 continue;
%             end
%             beam(j).elapsed(~oneGood) = interp1(ii, beam(j).elapsed(ii), find(~oneGood), 'linear');
%             beam(j).period = [diff(beam(j).elapsed(1:2)) diff(beam(j).elapsed)];
        else
            beam(j).elapsed = cumsum(beam(j).period + 1);%1 extra tick for the rollover
            updateinds = [];
            for k = 1:length(timer.rollover.counters)
                %a valid rollover has exactly 2^16 counters
                if (timer.rollover.ntimers(k) ~= 2^16)
                    continue;
                end
                timer_ind = timer.inds(timer.rollover.whichtimer(k) + (0:65535));
                timer_ticks = timer.tickCount(timer.rollover.whichtimer(k) + (0:65535));
                si = bsearch(beam(j).tagInd, timer_ind(1)); %find(beam(j).tagInd <= timer_ind(1), 1, 'last');
                ei = min(length(beam(j).tagInd), bsearch(beam(j).tagInd, timer_ind(end)) + 1); %find(beam(j).tagInd >= timer_ind(end), 1, 'first');
                matching_tag = interp1(beam(j).tagInd(si:ei), si:ei, timer_ind, '*nearest', NaN);
                %             if (any(~isfinite(matching_tag)))
                %                 warning (['something odd happened with counter interval ' num2str(k)]);
                %             end
                
                timer_ticks = timer_ticks(isfinite(matching_tag));
                elapsed = beam(j).elapsed(matching_tag(isfinite(matching_tag)));
                tt0 = mean(timer_ticks);
                
                el0 = mean(elapsed);
                p = polyfit(elapsed - el0, timer_ticks - tt0, 1);
                timer.rollover.scalingError(j,k) = 1 - p(1);
                timer.rollover.offset(j,k) = p(2) - el0 + tt0;
                if (k == 1 || isempty(updateinds))
                    updateinds = 1:ei;
                else
                    updateinds = (updateinds(end)+1):ei;
                end
                if (k == length(timer.rollover.counters))
                    updateinds = updateinds(1):length(beam(j).elapsed);
                end
                beam(j).elapsed(updateinds) =  beam(j).elapsed(updateinds) + timer.rollover.offset(j,k);
            end
        end
        %extend tag indexing to 1 previous so that elapsed time can be
        %approximated for pre-tag photons - this is mainly a fix to avoid
        %NaNs 
        tI = [-1 beam(j).tagInd];
        te = [beam(j).elapsed(2)-2*beam(j).elapsed(1), beam(j).elapsed]; 
        
        green{j} = find((code == 1+delta(j) | code == 3+delta(j)) & valid);
        beam(j).gticks = data(green{j});
        beam(j).gelapsed = interp1(tI, te, green{j}, '*previous', 'extrap') + beam(j).gticks;
        beam(j).gPhase = 2*pi*beam(j).gticks ./ interp1( beam(j).tagInd, beam(j).period, green{j}, '*nearest', 'extrap');
        if (any(~isfinite(beam(j).gPhase)))
            warning('something foobar');
        end
        
        red{j} = find((code == 2+delta(j) | code == 3+delta(j)) & valid);
        beam(j).rticks = data(red{j});
        beam(j).relapsed = interp1( tI, te, red{j}, '*previous', 'extrap') + beam(j).rticks;
        beam(j).rPhase = 2*pi*beam(j).rticks ./ interp1( beam(j).tagInd, beam(j).period, red{j}, '*nearest', 'extrap');
        
            
        

      
    end
   
end

function [tickCount, offset, oneGood, validTimer] = alignTimes (timer, inds, tickVal, maskForAgreement, maskForAddition)
    %function alignTimes (timer, inds, tickVal, maskForAgreement, maskForAddition)
    %timerStruct = timer structure which contains fields inds, validTimer, tickCount
    %
    % inds = index in chunk locating timer/aod
    % tickVal = value of lower bits of global tick counter at time of
    %   timer/aod [for aod it's the lowest 16 bits of tick clock, for timer
    %   it's lowest 14 bits but with lowest 2 bits always 0] 
    % maskForAgreement = mask indicating bits that should be the same in
    %   timer.tickCounts and tickVal
    % maskForAddition = mask eliminating bits in timer that are covered by
    %   tickVal -- tickCount = timer.tickCount&mask + tickVal
    %
    % tickCount = best value of tick counter for each timer/aod
    % offset = how far ahead or behind of the the previous timer the
    %   matching timer was
    % oneGood = whether we were able to find a matching timer
    % validTimer = whether the matching timer was "valid," ie one where we
    %   were sure the tick count of the timer is correct
            tiprev = interp1(timer.inds, 1:length(timer.inds), inds, '*prev', 'extrap'); %index before the timer in the file
            tiprev(~isfinite(tiprev)) = 1;
            
            offsettrials = [0 reshape([-1:-1:-8;1:8], 1,[])]; %previous, 2 before, 1 after, 3 before, 2 after etc.
            
            oneGood = false(size(inds));
            tickCount = NaN(size(inds));
            offset = NaN(size(inds));
            validTimer = true(size(inds));
            for j = 1:length(offsettrials)
                ti = min(max(tiprev + offsettrials(j), 1), length(timer.inds));
                timerVal = timer.tickCount(ti);                
                timerGood =  bitand(uint64(timerVal), maskForAgreement) ==  bitand(uint64(tickVal), maskForAgreement);
                update = ~oneGood & timerGood;
                oneGood(update) = true;
                tickCount(update) = double(bitand(uint64(timerVal(update)),  maskForAddition)) + tickVal(update);
                offset(update) = offsettrials(j);
                validTimer(update) = timer.validTimer(ti(update));
                if (all(oneGood))
                    break;
                end
            end
            
            %use the previous timer even though it doesn't match
             %leave offset as NaN and oneGood as false
            update = ~oneGood;
            ti = tiprev;
            timerVal = timer.tickCount(ti);                
            tickCount(update) = double(bitand(uint64(timerVal(update)),  maskForAddition)) + tickVal(update);
            validTimer(update) = timer.validTimer(ti(update));
            
            
    
end



function bxor = bitxorarray(arr)
bxor = arr(1);
for j = 2:length(arr)
    bxor = bitxor(bxor, arr(j));
end
end
