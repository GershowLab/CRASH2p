classdef VOISelectorObject < handle
     % helper class for VoiSelector app
    %   contains input data and space for return values
    %   input: xaxis,yaxis,zaxis - axes for volume data 
    %          volume - 1 volume or cell of volumes, dimension of volume
    %          matches xaxis,yaxis,zaxis sizes (nx,ny,nz)
    %           label - if more than 1 volume, cell of labels that describe
    %           volume (e.g. red, green) 
    %   output: vois (fields x,y,z): Voi objects that define volumes of
    %   interest
    %   finished: a boolean flag the VoiSelector app sets to true when
    %       output fields are filled in and we are ready to proceed
    
    
    properties
        xaxis
        yaxis
        zaxis
        volume %if more than one volume, make a cell {best100, random100, best5} etc.
        label %if more than one volume, then these are their names, e.g. best, random, etc.
        vois

        finished = false;
        imsize;

        srcdir = '';
        dstfile = '';

    end
    
    methods
        
        function obj = VOISelectorObject(xaxis,yaxis,zaxis,volume,label)
            if nargin==1 % either given template variable or path to template.mat
                try
                    if ischar(xaxis) % given path, try loading from disk
                        if isfile(xaxis)
                            obj.srcdir = fileparts(xaxis);
                            load(xaxis,'final_template');
                        elseif isfolder(xaxis)
                            obj.srcdir = xaxis;
                            d = dir(fullfile(obj.srcdir, '*_template.mat'));
                            load(fullfile(d.folder, d.name), 'final_template');
                        end
                    else % given template variable
                        if isfield(xaxis,'final_template')
                            final_template = xaxis.final_template;
                        else
                            final_template = xaxis;
                        end
                    end
                    xaxis = final_template.u;
                    yaxis = final_template.v;
                    zaxis = final_template.w;
                    volume = final_template.F;
                catch me
                    disp('VOISelectorObject(dir containing template.mat, file containing template.mat, or template');
                    rethrow(me);
                end
            end
            
            obj.xaxis = xaxis;
            obj.yaxis = yaxis;
            obj.zaxis = zaxis;
            obj.volume = volume;
            if (~iscell(obj.volume))
                obj.imsize = size(obj.volume);
            else
                obj.imsize = size(obj.volume{1});
            end
            obj.finished = false;
            if (~existsAndDefault('label',{}) || isempty(label))
                obj.label = {};
            else
                obj.label = label;
            end

            obj.vois = {};
            obj.dstfile = fullfile(obj.srcdir, 'voi_selections.mat');
        end
        
        function obj = VOISelectorObjectOLD(xaxis, yaxis, zaxis, volume, label)
            if (nargin == 1)
                try
                    if (ischar(xaxis))
                        if exist(xaxis, 'dir')
                            obj.srcdir = xaxis;
                            d = dir(fullfile(obj.srcdir, '*_tm.mat')); %template maker
                            load(fullfile(d.folder, d.name), 'tm');
                        else
                            if exist(xaxis, 'file')
                                obj.srcdir = fileparts(xaxis);
                                load(xaxis, 'tm');
                            end
                        end
                    else
                        if (is(xaxis, 'TemplateMaker'))
                            tm = xaxis;
                        end
                    end
                catch me
                    disp('VOISelectorObject(dir containing template maker, file containing template maker, or templatemaker')
                    rethrow(me);
                end
                xaxis = tm.template.u;
                yaxis = tm.template.v;
                zaxis = tm.template.w;
                volume = tm.template.F;

            end
               
            



            %ROISELECTORINPUT Construct an instance of this class
            %   Detailed explanation goes here
            obj.xaxis = xaxis;
            obj.yaxis = yaxis;
            obj.zaxis = zaxis;
            obj.volume = volume;
            if (~iscell(obj.volume))
                obj.imsize = size(obj.volume);
            else
                obj.imsize = size(obj.volume{1});
            end
            obj.finished = false;
            if (~existsAndDefault('label',{}) || isempty(label))
                obj.label = {};
            else
                obj.label = label;
            end

            obj.vois = {};
            obj.dstfile = fullfile(obj.srcdir, 'voi_selections.mat');

        end

        function addVoi (obj, voi)
            obj.vois{end+1} = voi;
        end

        function deleteVoi(obj, ind)
            if (ind < 1 || ind > length(obj.vois))
                return;
            end
            obj.vois = obj.vois([1:(ind-1) ind+1:end]);

        end

        function labelim = makeLabelMap (obj)
            %if vois overlap, map uses largest value
            labelim = zeros(obj.imsize);
            for j = 1:length(obj.vois)
                labelim(obj.vois{j}.makeMask) = j;
            end
        end

        function bwim = getMask(obj, ind)
            bwim = obj.vois{ind}.makeMask;
            if (isempty(bwim))
                bwim = false(obj.imsize);
            end
        end

        function saveVoiSet(obj, nickname)
            voiset.nickname = nickname;
            voiset.vois = obj.vois;
            matobj = matfile(obj.dstfile,'Writable',true);
            if (~isempty(who(matobj, 'voisets')))
                matobj.voisets = [matobj.voisets voiset];
            else
                matobj.voisets = voiset;
            end
        end


    end
    
end

