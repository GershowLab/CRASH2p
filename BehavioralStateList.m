classdef BehavioralStateList
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stateNames = {};
        stateVector;
        nframes = 0;
    end
    
    methods
        function bsl = BehavioralStateList(varargin)
            bsl.stateNames = regexprep(varargin, '\s+', '_');
        end
        
        function valid = getValidStates(bsl, trueStateNames, falseStateNames)
        %function valid = getValidStates(bsl, trueStateNames, falseStateNames)
        %combines trueStates with AND and falseStates with AND NOT
            if ~iscell(trueStateNames)
                trueStateNames = {trueStateNames};
            end
            existsAndDefault('falseStateNames', {});
            if (~iscell(falseStateNames))
                falseStateNames = {falseStateNames};
            end
            valid = true([1 bsl.nframes]);
            for j = 1:length(trueStateNames)
                valid = valid & [bsl.stateVector.(trueStateNames{j})];
            end
            for j = 1:length(falseStateNames)
                valid = valid & ~[bsl.stateVector.(falseStateNames{j})];
            end
            
        end
            
        
        function bsl = addState(bsl, stateName)
        %function bsl = addState(bsl, stateName)
            if (iscell(stateName))
                for j = 1:length(stateName)
                    bsl = bsl.addState(stateName{j});
                end
                return;
            end
        
            stateName = regexprep(lower(stateName), '\s+', '_');
            if (bsl.nframes <= 0)
                return;
            end
            if(~any(strcmpi(bsl.stateNames, stateName)))
                bsl.stateVector(bsl.nframes).(stateName) = false;
                [bsl.stateVector.(stateName)] = deal(false);
                bsl.stateNames = union(bsl.stateNames, stateName);
            end
            
        end
        
        function bsl = processCellList(bsl, stateCellList, merge)
        %function bsl = processCellList(bsl, stateCellList, merge)
        %stateCellList is of form {{index1, flag1, flag2}, {index1, flag2,
        %flag3}} 
        %if index is excluded, then assumes jth element is for jth frame
        %if merege is true, then old state vector is not overwritten
        %but each assigned state is overwritten
            existsAndDefault('merge', false);
            scl = [stateCellList{:}];
            names = scl(cellfun(@ischar, scl));
            names = unique(lower(names));
            names = regexprep(names, '\s+', '_');
            bsl.stateNames = union(bsl.stateNames, names);
            nf = max(cell2mat(scl(cellfun(@isnumeric, scl))));
            
            if (isempty(nf))
                nf = length(stateCellList);
            else
                nf = max(length(stateCellList), nf);
            end
            if (merge)
                bsl.nframes = max(nf, bsl.nframes);
            else
                bsl.nframes = nf;
            end
            
            for j = 1:length(bsl.stateNames)
                s.(bsl.stateNames{j}) = false;
            end
            if (~merge)
                bsl.stateVector = repmat(s, [bsl.nframes 1]);
            end
                
            for j = 1:length(stateCellList)
                scl = stateCellList{j};
                if (isnumeric(scl{1}))
                    ind = scl{1};
                    scl = scl(2:end);
                else
                    ind = j;
                end
                bsl = bsl.setStateFrame(ind, scl);
            end
            bsl.stateNames = union(bsl.stateNames, fieldnames(bsl.stateVector));        
        end
        
        function bsl = setStateFrame(bsl, ind, stateCell)
            %function bsl = setStateFrame(bsl, ind, stateCell)
            %ind = frame#
            %stateCell = {flag1, flag2} is list of flags to associate with that
            %state (in text)
            %also OK bsl = bsl.setStateFrame({ind, flag, flag});
            if (iscell(ind))
                if length(ind) > 1
                    stateCell = ind{2:end};
                else
                    stateCell = {};
                end
                ind = ind{1};
            end
            if (~iscell(stateCell))
                stateCell = {stateCell};
            end
            stateCell =  regexprep(stateCell, '\s+', '_');
            fn = bsl.stateNames;
            for k = 1:length(fn)
                if (~isfield(bsl.stateVector, fn{k}))
                     bsl.stateVector(ind).(fn{k}) = false;
                    [bsl.stateVector.(fn{k})] = deal(false);
                end
                bsl.stateVector(ind).(fn{k}) = false;
            end
            for k = 1:length(stateCell)
                if (~isfield(bsl.stateVector, stateCell{k}))
                    [bsl.stateVector.(stateCell{k})] = deal(false);
                    bsl.stateNames = union(bsl.stateNames, stateCell{k});
                end
                bsl.stateVector(ind).(stateCell{k}) = true;
            end
            bsl.nframes = max(bsl.nframes, ind);
        end
        
        function stateCell = getStateFrame(bsl, ind, name)
        %function stateCell = getStateFrame(bsl, ind)
        %function stateValue = getStateFrame(bsl, ind,stateName)
        
        
            if (length(ind) == 1)
                
                if (ind > length(bsl.stateVector) || ind < 1)
                    stateCell = {};
                    return;
                end
                if (nargin >=3)
                    stateCell = bsl.stateVector(ind).(name);
                    return;
                end
                fn = fieldnames(bsl.stateVector);
                stateCell = fn(logical(cell2mat(struct2cell(bsl.stateVector(ind))))); 
            else
                for j = 1:length(ind)
                    stateCell{j} = bsl.getStateFrame(ind(j));%#ok<AGROW>
                end
            end
        end
        
        function scl = getStateCellList(bsl)
            %function scl = getStateCellList(bsl)
            scl = bsl.getStateFrame(1:length(bsl.stateVector));
        end
        
        function toCSV(bsl, filename, condensed)
        % function toCSV(bsl, filename, condensed)
        % condensed means do not write out frames where no behavior is
        % flagged
            existsAndDefault('condensed', false);
            scl = bsl.getStateCellList();
            fid = fopen(filename, 'w');
            try
                for j = 1:length(scl)
                    if (~isempty(scl{j}))
                        fprintf(fid, '%d', j);
                        fprintf(fid, ', %s', scl{j}{:});
                        fprintf(fid, '\n');
                    else
                        if (~condensed)
                            fprintf(fid, '%d\n', j);
                        end
                    end
                end
            catch me
                fclose(fid);
                rethrow(me);
            end
            fclose(fid);
        end
        function bsl = addCSV(bsl, filename, merge)
        %function bsl = addCSV(bsl, filename, merge)
            fid = fopen(filename, 'r');
            try
                j = 1;
                while true
                    tline = fgetl(fid);
                    if (~ischar(tline))
                        break;
                    end
                    scl{j} = regexp(tline, ',\s*', 'split'); %#ok<AGROW>
                    if ~isempty(regexp(scl{j}{1}, '\d', 'ONCE'))
                        scl{j}{1} = str2double(scl{j}{1}); %#ok<AGROW>
                    end
                    j = j+1;
                end
                bsl = bsl.processCellList(scl, merge);
                
            catch me
                fclose(fid);
                rethrow(me);
            end
            fclose(fid);
        end
        
        function bsl = setLength(bsl, nframes)
            bsl.nframes = nframes;
            if (length(bsl.stateVector) > nframes)
                bsl.stateVector = bsl.stateVector(1:nframes);
            else
                bsl = bsl.setStateFrame(nframes, {});
            end
        end
    end
    
    methods (Static)
        function bsl = fromCSV(filename)
            %function bsl = fromCSV(filename)
            
            bsl = BehavioralStateList();
            bsl = bsl.addCSV(filename, false);
            
                    
        end
    end
    
end

