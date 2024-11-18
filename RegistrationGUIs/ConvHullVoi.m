classdef ConvHullVoi < Voi
    %uses matlab alphaShape 

    properties
       pts %list of points on convex hull, Nx3 (surely works) or Nx2 (may work)
       as %shape object
       dt %delaunay triangulation of 
       ch %triangulation of convex hull
       chxy %list of points on convex hull projected onto xy plane
       trxy %triangulation of 3D convex hull projected onto plane (i.e. how to locate facets from x,y point) 
       height %interpolant that finds the thickness of the convex hull 
    end

    methods

        function obj = transform(obj, at, invert)
            %applies an affine transform to the convex hull and returns a
            %new object

            %at is 3x4 (or 4x4 if made invertible)
            %pts is Nx3
            
            newpts = obj.pts'; %(3xN)
            newpts(4,:) = 1;  %4xN
            if (nargin > 2 && invert)
                newpts = (at\newpts)'; %Nx4 or Nx3 depending on at 
            else
                newpts = (at*newpts)';
            end
            obj = obj.setPoints(newpts(:,1:3));

        end

        function obj = setPoints(obj, ptlist)
            if (size(ptlist,2) ~= 2 && size(ptlist,2) ~=3)
                ptlist = ptlist';
            end
            if (size(ptlist,2) ~= 2 && size(ptlist,2) ~=3)
                error ('need a N x D list of points, where D = 2 or 3')
            end
            
            k = convhull(ptlist);
            obj.pts = ptlist(unique(k(:)),:);%all points that appear on the convex hull at least onc
            obj.as = alphaShape(obj.pts,Inf); %Inf = convex hull; Specify a = criticalAlpha(shp,'one-region') to use the smallest alpha radius that produces an alpha shape with only one region. HoleThreshold=100 means fill in interior holes
            obj.dt = delaunayTriangulation(obj.pts);
           
            k = convhull(ptlist(:,1:2));
            obj.chxy = ptlist(k,1:2);

            obj.ch = triangulation(obj.dt.convexHull, obj.pts);
            obj.trxy = triangulation(obj.dt.convexHull,obj.dt.Points(:,1:2));

            [xx,yy] = ndgrid(obj.xaxis, obj.yaxis);
            zz = reshape(obj.getZHeight([xx(:) yy(:)]),size(xx));
            obj.height = griddedInterpolant(xx,yy,zz,'linear','nearest');

        end

        function obj = addPoints(obj, ptlist)
            if (size(ptlist,2) ~= size(obj.pts,2))
                ptlist = ptlist';
            end
            if (size(ptlist,2) ~= size(obj.pts,2))
                error ('added points must have same dimension as existing points')
            end
            obj = obj.setPoints([obj.pts;ptlist]);
        end

        function tf = inShape3D(obj, ptlist)
            %whether points are inside the 3D object
            if (size(ptlist,2) ~=3)
                ptlist = ptlist';
            end
            if (size(ptlist,2) ~=3)
                error ('need a N x 3 list of points')
            end
            tf = obj.as.inShape(ptlist(:,1), ptlist(:,2), ptlist(:,3));
        end

        function tf = inShape2D(obj, ptlist)
            %whether an infinite z-line passing through a point intersects
            %the object
            if (size(ptlist,2) ~= 2 )
                ptlist = ptlist';
            end
            if (size(ptlist,2) ~= 2 )
                error ('need a N x 2 list of points')
            end
            tf = inpolygon(ptlist(:,1), ptlist(:,2), obj.chxy(:,1), obj.chxy(:,2));
        end

        function z = getZ(obj,ID,pt2d)
            if (length(ID) > 1)
                z = 0*ID;
                for j = 1:length(ID)
                    z(j) = obj.getZ(ID(j),pt2d);
                end
                return;
            end
            fn = obj.ch.faceNormal(ID);
            pt1 = obj.ch.Points(obj.ch.ConnectivityList(ID,1),:);
            %(p - p1).fn = 0
            %x fn_x + y fn_y + z fn_z = p1.fn
            %z = (p1.fn - xfn_x - yfn_y)/fn_z
            z = (dot(pt1,fn) - dot(fn(1:2),pt2d(1:2)))/fn(3);
        end

        function h = getZHeight(obj, pt2d)
            h = diff(obj.getZBounds(pt2d),[],2);
        end

        function zbounds = getZBounds(obj, pt2d)
            if(size(pt2d,1) > 1)
                if (size(pt2d,2) == 1)
                    pt2d = pt2d';
                else
                    zbounds = 0*pt2d;
                    valid = obj.inShape2D(pt2d);
                    for j = find(valid(:)')
                        zbounds(j,:) = obj.getZBounds(pt2d(j,:));
                    end
                    return;
                end
            end

%             if (~inpolygon(pt2d(1), pt2d(2), obj.chxy(:,1), obj.chxy(:,2)))
%                 zbounds = [0 0];
%                 return;
%             end
            ntri = size(obj.trxy.ConnectivityList,1);
            x = obj.trxy.Points(:,1);
            x = x(obj.trxy.ConnectivityList);
            y = obj.trxy.Points(:,2);
            y = y(obj.trxy.ConnectivityList);
            delta = 1e-12;
            inbox = find(min(x,[],2)-delta <= pt2d(1) & max(x,[],2)+delta >= pt2d(1) &   min(y,[],2)-delta <= pt2d(2) & max(y,[],2)+delta >= pt2d(2) );

            bc = obj.trxy.cartesianToBarycentric(inbox(:), repmat(pt2d,[numel(inbox) 1]));
            intri = inbox(all(bc >= -delta,2)); %account for numerical imprecision
%             if (numel(intri) < 2)
%                 intri = inbox(all(bc >= -delta,2)); %account for numerical imprecision
%             end
            z = obj.getZ(intri,pt2d);
            zbounds = [min(z(isfinite(z))) max(z(isfinite(z)))];
        end




    end
    methods 
        function bwim = makeMask(obj)
            bwim = reshape(obj.as.inShape([obj.xx(:),obj.yy(:),obj.zz(:)]), size(obj.xx));
        end
    end

    methods (Static)
        
        function cv = convertVoi(voi)
            %converts any voi object to convex hull voi
            bwim = voi.makeMask();
            cv = ConvHullVoi(voi.xaxis, voi.yaxis, voi.zaxis);
            cv = cv.setPoints([voi.xx(bwim), voi.yy(bwim), voi.zz(bwim)]);
        end



    end

end