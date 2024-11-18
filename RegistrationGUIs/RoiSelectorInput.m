classdef RoiSelectorInput < handle
     % helper class for RoiSelector3D app
    %   contains input data and space for return values
    %   input: xedges,yedges,zedges - edges for volume data (size
    %           nx+1,ny+1,nz+1)
    %          volume - 1 volume or cell of volumes, dimension of volume
    %          matches xaxis,yaxis,zaxis sizes (nx,ny,nz)
    %           label - if more than 1 volume, cell of labels that describe
    %           volume (e.g. best 10, evenly spaced 30) 
    %   output: roi (fields x,y,z): the limits in microns of the 3D rectangle roi from
    %       which to make the template and proceed with non-rigid registration
    %   finished: a boolean flag the TemplateCreator app sets to true when
    %       output fields are filled in and we are ready to proceed
    %
    %  typical use (see getRoi function):
    %       rsi = RoiSelectorInput(xaxis,yaxis,zaxis,volume,label);
    %       roi = rsi.getRoiGui(); % includes gui initialization and
    %                                waitfor command
    
    properties
        xedges
        yedges
        zedges
        volume %if more than one volume, make a cell {best100, random100, best5} etc.
        label %if more than one volume, then these are their names, e.g. best, random, etc.
        roi

        finished = false;
    end
    
    methods
        function obj = RoiSelectorInput(xedges, yedges, zedges, volume, label)
            %ROISELECTORINPUT Construct an instance of this class
            %   Detailed explanation goes here
            obj.xedges = xedges;
            obj.yedges = yedges;
            obj.zedges = zedges;
            if (~iscell(volume))
                obj.volume = padarray(volume,[1 1 1],0,'post');
            else
                for j = 1:length(volume)
                    obj.volume{j} = padarray(volume{j},[1 1 1],0,'post');
                end
            end
            obj.finished = false;
            if (~existsAndDefault('label',{}) || isempty(label))
                obj.label = {};
            else
                obj.label = label;
            end

        end

        function roi = snapRoiToGrid(rsi, roi)
            roi.x = RoiSelectorInput.snapRangeToGrid(rsi.xedges, roi.x);
            roi.y = RoiSelectorInput.snapRangeToGrid(rsi.yedges, roi.y);
            roi.z = RoiSelectorInput.snapRangeToGrid(rsi.zedges, roi.z);
            
        end
        
        function roi = getRoiGui(rsi)
            rsi.finished = false;
            RoiSelector3D(rsi);
            waitfor(rsi, 'finished', true)
            roi = rsi.roi;
        end
    end
    methods(Static)
        function roi = getRoi(xedges, yedges, zedges, volume, label)
            existsAndDefault('label',{});
            rsi = RoiSelectorInput(xedges,yedges,zedges,volume,label);
            roi = rsi.getRoiGui();
        end

        function rng = snapRangeToGrid(grid,rng)
            %last grid element <= rng(1) and first grid element >= rng(2),
            %taking care of out of range values

            gi = griddedInterpolant(grid,grid,'previous','nearest');
            rng(1) = gi(rng(1));
            gi = griddedInterpolant(grid,grid,'next','nearest');
            rng(2) = gi(rng(2));
        end



    end
end

