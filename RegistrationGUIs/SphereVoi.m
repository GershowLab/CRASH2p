classdef SphereVoi < Voi
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        pt
        radius        
    end

    methods
        function bwim = makeMask(obj)
            
            bwim = obj.volumeWithinDistance(obj.radius);

        end


        function [obj,r] = autoSetRadius(obj, im, r_range)
            if (nargin < 3)
                r_range = [1 10];
            end
            r_range = max(r_range,0);
            r = fminbnd(@(r) obj.segmentObjective(r,im), min(r_range), max(r_range));
            obj.radius = r;
        end

        function obj = transform(obj, at, invert)

            %at is 3x4 (or 4x4 if made invertible)
            %pts is Nx3
            
            newpt = obj.pt(:); %(3x1)
            newpt(4,:) = 1;  %4x1
            at(4,:) = [0 0 0 1];
            if (nargin > 2 && invert)
                npt = (at\newpt); 
            else
                npt = (at*newpt);
            end
            obj.pt = npt(1:3);
        end

        function tf = inShape3D(obj, ptlist)
            %whether points are inside the 3D object
            if (size(ptlist,2) ~=3)
                ptlist = ptlist';
            end
            if (size(ptlist,2) ~=3)
                error ('need a N x 3 list of points')
            end
            tf = (ptlist(:,1)-obj.pt(1)).^2 + (ptlist(:,2)-obj.pt(2)).^2 + (ptlist(:,2)-obj.pt(2)).^2 < obj.radius.^2;
        end

        function h = height(obj, x, y)
            h = 2*sqrt(max(obj.radius.^2 - (obj.pt(1)-x).^2 - (obj.pt(2)-y).^2,0));
        end


        function val = segmentObjective(obj, r, im)
            %minimize average outside minus average inside; outside being a
            %volume twice as large as inside
            bwim1 = obj.volumeWithinDistance(r);
            bwim2 =  obj.volumeWithinDistance(r*1.2599) & ~bwim1;

            if (~any(bwim2,'all') || ~any(bwim1,'all'))
                val = 0;
                return;
            end

            val = sum(im(bwim2))/nnz(bwim2) - sum(im(bwim1))/nnz(bwim1);
        end



        function bwim = volumeWithinDistance(obj, r)
            xa = obj.pt(1);
            xc = obj.xx;
            ya = obj.pt(2);
            yc = obj.yy;
            za = obj.pt(3);
            zc = obj.zz;
            bwim = ((xc-xa).^2 + (yc-ya).^2 + (zc-za).^2) < r^2;
        end

    end

    

end

%                 valid = cv(j).inShape3D(ploc);
%                 counts(j) = sum(abs(sin(pphase(valid))))/0.7846; %correction used in fast assembler
%                 valid = cv(j).inShape3D(ploc_alt);
%                 altcounts(j) = sum(abs(sin(pphase_alt(valid))))/0.7846; %correction used in fast assembler
%             end
% 
%             if (nargout < 3)
%                 return;
%             end
% 
%            [pos, ~, period] = vra.getDwellLocations(timerange, false);
%            [pos_alt, ~, period_alt] = vra.getDwellLocations(timerange, true);
% 
%            tau_per_um = period/(4*vra.z_scale);
%            tau_per_um_alt = period_alt/(4*vra.z_scale);
%            
%            for j = 1:numel(cv)
%                %h = cv(j).getZHeight(pos(1:2,:)');
%                h = cv(j).height(pos(1,:), pos(2,:));
%                dwell(j) = sum(tau_per_um(:).*h(:));
%                %h = cv(j).getZHeight(pos_alt(1:2,:)');
%                h = cv(j).height(pos_alt(1,:), pos_alt(2,:));
%                altdwell(j) = sum(tau_per_um_alt(:).*h(:));
%            end
%         end