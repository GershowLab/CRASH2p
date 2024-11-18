function [timer, lastInd] = findTimers(chunk)
% [timer, beam, aod, diagnostics, lastInd] = findTimers(chunk)
% locates timers only for purposes of later file parsing 
%
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
    aodarr = chunk(aodlocinds);
    
    %possible ways you can get code 8, 12, 13
    %due to overflow, due to aod trigger, in location for either of these
   
    %first let's check xsums
    xs = timerarr(1,:);
    for q = 2:7
        xs = bitxor(xs, timerarr(q,:));
    end
    locvalid = (xs == 0);
    
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
    aodlocinds = repmat(aodinds, [3 1]) + repmat((1:3)', [1 size(aodinds,2)]);
   
    allloc = [locinds(:); aodlocinds(:)];
    %first, let's verify that all bad checksums are in locations somewhere
    
    badaodinds = setdiff(badaodinds, allloc);
    badtimerinds = setdiff(badtimerinds, allloc);
    badinds = union(badaodinds, badtimerinds);
    if ~isempty(badinds)
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
            warning(['found ' num2str(nnz(~allzero)) ' aod valid checksums inside timer location block; this will happend randomly 1/65535 times so may be ok ']);
        end
    end
    
    timerinaod = intersect(timerinds, aodlocinds);
    if (~isempty(timerinaod))
        %this should always be an error ?
        timerinds = setdiff(timerinds, timerinaod);
        locinds = repmat(timerinds, [7 1]) + repmat((1:7)', [1 size(timerinds,2)]);
        warning(['found ' num2str(length(timerinaod)) ' timer valid checksums inside aod location block ']);
    end
        
    timerarr = chunk(locinds);
 
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
            warning ('reset to 0 at non-rollover point - expect counter to go to 4095 and reset');
        end
        timer.counter((I+1):end) = timer.counter((I+1):end) + 4096;
    end
  
    tcedges = min(timer.counter):(max(timer.counter) + 1);
   
    ncounters = histcounts(timer.counter, tcedges); tcedges = tcedges(1:end-1);
    problem = ncounters ~= 2^16;% & tcedges > min(timer.counter) & tcedges < max(timer.counter); %we will just discard the end counters anyway
    problem([1 end]) = false;
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
    if (~isempty(aodinds))
        lastInd = max(max(aodinds + 3), max(timerinds)+7,'omitnan');
    else
        lastInd =  max(timerinds,[],'omitnan') +7;
    end
    if (isempty(lastInd))
        lastInd = length(chunk);
    end

end