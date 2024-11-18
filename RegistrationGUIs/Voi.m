classdef (Abstract) Voi
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        xaxis
        yaxis
        zaxis

        xx
        yy
        zz
    end

    methods
        function obj = Voi(varargin)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.xaxis = varargin{1};
            obj.yaxis = varargin{2};
            obj.zaxis = varargin{3};
            if ((length(varargin) < 4) || varargin{4})
                [obj.xx,obj.yy,obj.zz] = ndgrid(obj.xaxis, obj.yaxis, obj.zaxis);
            end
        end
        
        function im = labelImage(obj, val, im)
            if nargin < 3
                im = zeros(size(obj(1).makeMask));
            end
            for j = 1:length(obj)
                im(obj(j).makeMask) = val(j);
            end
        end

        function pt = centerOfMass(obj)
            if (length(obj) > 1)
                pt = zeros(3,length(obj));
                if (iscell(obj))
                    for j = 1:length(obj)
                        pt(:,j) = obj{j}.centerOfMass;
                    end
                else
                    for j = 1:length(obj)
                        pt(:,j) = obj(j).centerOfMass;
                    end
                end
                return
            end

            bwim = obj.makeMask();


            if (isempty(obj.xx) || any(size(obj.xx) ~= size(bwim)))
                [obj.xx,obj.yy,obj.zz] = ndgrid(obj.xaxis, obj.yaxis, obj.zaxis);
            end
            pt = [mean(obj.xx(bwim));mean(obj.yy(bwim));mean(obj.zz(bwim))];

        end

    end
    methods (Abstract)
        bwim = makeMask(obj);
    end

    methods (Static)
        function bwim = makeLabeledMask(voiarr)
            if (iscell(voiarr))
                bwim = 0.0*voiarr{1}.makeMask();
                for j = 1:length(voiarr)
                    bwim(voiarr{j}.makeMask()) = j;
                end
            else
                bwim = double(voiarr(1).makeMask());
                for j = 1:length(voiarr)
                    bwim(voiarr(j).makeMask()) = j;
                end
            end
    end



    end

end