function [at, pvec] = orientPoints_xyrot(x,y,z,w,thetaz)
%function [at] = orientPoints3D(x,y,z,w)
%function [at] = orientPoints3D(xaxis,yaxis,zaxis,im3d)
%pvec = [theta_z, theta_y, theta_x, 0, 0, 0]
% transform order is rot about z, rot about x, rot about y
% or (in homogenous coordinates) at = ry*rx*rz 
% keeps origin fixed

if (nargin < 4)
    w = ones(size(x));
end

pvec = zeros([1 6]);
pvec180 = pvec;


if (ndims(w) == 3)
    [xx,yy,zz] = ndgrid(x,y,z);
    x = xx(:)';
    y = yy(:)';
    z = zz(:)';
    w = w(:)';
end

pinitial = [x;y;z;ones(size(x))];

w = w/mean(w);

x2 = x - mean(w.*x);
y2 = y - mean(w.*y);
z2 = z - mean(w.*z);


uz = [cos(thetaz), -sin(thetaz); sin(thetaz), cos(thetaz)];


pts2 = uz*[x;y];
dx = -mean(w.*pts2(1,:));
dy = -mean(w.*pts2(2,:));


atz = zeros(4,4);
atz(1:2,1:2) = uz;
atz(3,3) = 1;
atz(:,4) = [dx;dy;0;1];

pvec(1) = mod(atan2(uz(2,1), uz(1,1)), 2*pi);

% repeat for rotation around x axis

pts2 = atz*pinitial;
x = pts2(1,:);
y = pts2(2,:);
z = pts2(3,:);
x2 = x - mean(w.*x);
y2 = y - mean(w.*y);
z2 = z - mean(w.*z);

m1 = [mean(w.*(y2.^2)) mean(w.*y2.*z2);mean(w.*y2.*z2) mean(w.*(z2.^2))];
[ux,~,~] = svd(m1);
if (det(ux) < 0)
    ux = [1 0;0 -1]*ux;
end
%keep rotation in the range -pi/2 to pi/2
if (ux(1) < 0)
    ux = -ux;
end

pts2 = ux*[y;z];
dy = -mean(w.*pts2(1,:));
dz = -mean(w.*pts2(2,:));
atx = zeros(4,4);
atx(2:3,2:3) = ux;
atx(1,1) = 1;
atx(:,4) = [0;dy;dz;1];

pvec(3) = atan2(ux(2,1), ux(1,1));


% finally for y axis


pts2 = atx*atz*pinitial;
x = pts2(1,:);
y = pts2(2,:);
z = pts2(3,:);
x2 = x - mean(w.*x);
y2 = y - mean(w.*y);
z2 = z - mean(w.*z);

m1 = [mean(w.*(x2.^2)) mean(w.*z2.*x2);mean(w.*z2.*x2) mean(w.*(z2.^2))];
[uy,~,~] = svd(m1);
if (det(uy) < 0)
    uy = [1 0;0 -1]*uy;
end
%keep rotation in the range -pi/2 to pi/2

if (uy(1,1) < 0)
    uy = -uy;
end


pvec(2) = atan2(uy(1,2), uy(1,1));

at = RigidAligner3D.rigidTransform(pvec);

pvec180 = pvec;
pvec180(1) = mod(pvec180(1) + pi, 2*pi);
pvec180(2:3) = -pvec180(2:3);
at180 = RigidAligner3D.rigidTransform(pvec180);
