function [counts,dwell] = binWithCDMFFDShifts(plocvec,pphase,dlocvec,dweight,fa,shifts)

if iscell(plocvec)
    counts = zeros([[fa.getImSize] length(plocvec)]);
    dwell = zeros([[fa.getImSize] length(plocvec)]);
    for i=1:length(plocvec)
        [ct,dw] = binWithShifts(plocvec{i},pphase{i},dlocvec{i},dweight{i},fa,shifts);
        counts(:,:,:,i) = ct;
        dwell(:,:,:,i) = dw;
    end
else
    [counts,dwell] = binWithShifts(plocvec,pphase,dlocvec,dweight,fa,shifts);
end

end

function [counts,dwell] = binWithShifts(plocvec,pphase,dlocvec,dweight,fa,shifts)

% get counts
plocvec1 = applyCDMFFDShiftsToPts(plocvec,fa,shifts);
x = plocvec1(1,:)';
y = plocvec1(2,:)';
z = plocvec1(3,:)';
counts = fa.accumVals(fa.depthCorrection(pphase),x,y,z);

% get dwell
dlocvec1 = applyCDMFFDShiftsToPts(dlocvec,fa,shifts);
x = dlocvec1(1,:)';
y = dlocvec1(2,:)';
z = dlocvec1(3,:)';
dwell = fa.accumVals(dweight,x,y,z);

end