function locvec1 = applyCDMFFDShiftsToPts(locvec0,fa,shifts)

x = locvec0(1,:)';
y = locvec0(2,:)';
z = locvec0(3,:)';

Tx = shifts.Tx;
Ty = shifts.Ty;
Tz = shifts.Tz;

[~,subs0] = fa.accumVals(ones(size(x)),x,y,z);
locvec1 = locvec0;
subs1 = subs0;
for i=1:size(locvec1,2)
    if any(~isfinite(subs1(i,:)))
        continue;
    end
    thisploc = locvec1(:,i);
    ix = subs1(i,1);
    iy = subs1(i,2);
    iz = subs1(i,3);
    thisT = [Ty(ix,iy,iz);-Tx(ix,iy,iz);-Tz(ix,iy,iz)];
    locvec1(:,i) = thisploc+thisT;
end

end