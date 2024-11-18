function [seg,nseg] = getSegment(inds,varargin)
% get start and end inds from 1D array containing segments of continuous inds
% 
% transposed output - Rui 20200424

minlength = [];
maxskip = 0;

assignApplicable(varargin);

if ~isrow(inds)
    inds = inds';
end

startind = [];
endind = [];

try
    % check start
    if inds(2)-inds(1)>1
        startind = [startind inds(1)];
        endind = [endind inds(1)];
        inds(1) = [];
    end
    % check end
    if inds(end)-inds(end-1)>1
        startind = [startind inds(end)];
        endind = [endind inds(end)];
        inds(end) = [];
    end
catch
    seg = [];
    return;
end

dinds = diff(inds);

% find single-pt "segment" in the middle of inds
% check start
dinds1 = circshift(dinds,1);
% figure;
% subplot(3,1,1); plot(inds);
% subplot(3,1,2); plot(dinds);
% subplot(3,1,3); plot(dinds1);
found = find(dinds>1&dinds1>1);
startind = [startind inds(found)];
endind = [endind inds(found)];
inds(found) = [];
dinds(found) = [];

% find regular segments
% figure;
% subplot(2,1,1); plot(inds);
% subplot(2,1,2); plot(dinds);
found = find(dinds>maxskip+1);
startind = [startind inds(1) inds(found+1)];
endind = [endind inds(found) inds(end)];

seg = [sort(startind)' sort(endind)'];

% filter by minlength
if ~isempty(minlength)
    dseg = diff(seg,2);
    seg(dseg<minlength,:) = [];
end

nseg = size(seg,1);

end