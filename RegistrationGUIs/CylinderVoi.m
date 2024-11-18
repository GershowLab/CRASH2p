classdef CylinderVoi < Voi
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        pt1
        pt2
        radius
        inclEndCaps = true;
    end

    methods
        function bwim = makeMask(obj)
            
            bwim = obj.volumeWithinDistance(obj.radius, obj.inclEndCaps);

        end


        function [obj,r] = autoSetRadius(obj, im, r_range)
            if (nargin < 3)
                r_range = [1 10];
            end
            r_range = max(r_range,0);
            r = fminbnd(@(r) obj.segmentObjective(r,im), min(r_range), max(r_range));
            obj.radius = r;
        end

        function val = segmentObjective(obj, r, im)
            bwim1 = obj.volumeWithinDistance(r,obj.inclEndCaps);
            bwim2 =  obj.volumeWithinDistance(r*sqrt(2),obj.inclEndCaps) & ~bwim1;

            if (~any(bwim2) | ~any(bwim1))
                val = 0;
                return;
            end

            val = sum(im(bwim2))/nnz(bwim2) - sum(im(bwim1))/nnz(bwim1);
        end



        function bwim = volumeWithinDistance(obj, r, inclEndCaps)
            xa = obj.pt1(1);
            xb = obj.pt2(1);
            xc = obj.xx;
            ya = obj.pt1(2);
            yb = obj.pt2(2);
            yc = obj.yy;
            za = obj.pt1(3);
            zb = obj.pt2(3);
            zc = obj.zz;
            a = 0.5*sqrt(((xa - xc).*(yb - ya) - (xa-xb).*(yc-ya)).^2 + ...
                ((za - zc).*(yb - ya) - (za-zb).*(yc-ya)).^2 +...
                ((xa - xc).*(zb - za) - (xa-xb).*(zc-za)).^2);
            d = 2*a ./ sqrt((xa-xb).^2 + (ya - yb).^2 + (za - zb).^2);
            
            betweenPoints = (xc-xa).*(xb-xa) + (yc-ya).*(yb - ya) + (zc-za).*(zb-za) > 0 & ...
                            (xc-xb).*(xa-xb) + (yc-yb).*(ya - yb) + (zc-zb).*(za-zb) > 0;

            bwim = (d < r) & betweenPoints;

            if (inclEndCaps)
                bwim = bwim | (((xc-xa).^2 + (yc-ya).^2 + (zc-za).^2) < r^2) | (((xc-xb).^2 + (yc-yb).^2 + (zc-zb).^2) < r^2);

            end




        end

    end

    

end