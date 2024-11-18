function [S_merged,mergeinds] = detectAndMergeShortStrides(tsr,S,varargin)
%S_merged = tsr.detectAndMergeShortStrides(tsr.com.strides,varargin)
%tsr.com.strides = S_merged;
%
% default parameters:
%   nmean = 5           : compare each stride length to local mean of 5 strides
%                         (actually 2 before and 2 after, excluding itself)
%   mergethresh = 0.5   : fraction of local mean stride length; if stride
%                         length < half of local mean, then not a real stride
%   diffthresh = 0.25   : fraction of vfwd range within a stride; if end-start
%                         is smaller than this value then not sure whether
%                         to merge to previous or next stride

% default parameters
nmean = 5; % compare each stride length to local mean of 5 strides
mergethresh = 0.5; % fraction of local mean stride length; if stride length < half of local mean, then not a real stride
diffthresh = 0.25; % fraction of vfwd range within a stride; if end-start is smaller than this value then not sure whether to merge to previous or next stride

assignApplicable(varargin);

disp('detecting too-short strides and merging into its appropriate neighbor...');

behav_str = {'forward','backward'};

% add continuous crawl number and original stride number to struct
if ~isfield(S,'crawl_num')
    nstrides = length(S);
    crawl_num = 1;
    S(1).crawl_num = crawl_num;
    S(1).stride_num = 1;
    for j=2:nstrides
        S(j).stride_num = j;
        if S(j).t(1)-S(j-1).t(end)>tsr.dt+range(diff(tsr.tx)) || S(j).behav_label~=S(j-1).behav_label
            crawl_num  = crawl_num+1;
        end
        S(j).crawl_num = crawl_num;
    end
else
    crawl_num = max([S.crawl_num]);
end
fprintf('%d crawling bouts detected\n',crawl_num);

% detect and merge strides for each crawl
nphasebins = length(S(1).tt);
fields = fieldnames(S);
emptyc = cell(length(fields),1);
S_merged = S([]);
mergepairs = [];
for k=1:crawl_num
    valid_strides = find([S.crawl_num]==k);
    nvalidstrides = length(valid_strides);
    % 1) identify too-short strides
    S2 = S(valid_strides);
    q = mean([S2.behav_label]);
    if mod(q,1)~=0
        warning('crawling bout %d contains more than one type of behavior, shouldn''t happen',k);
        S_merged = cat(1,S_merged,S2);
        continue;
    end
    slen = NaN(1,nvalidstrides);
    for j=1:nvalidstrides
        slen(j) = length(S2(j).t)*tsr.dt;
    end
    slen_mean = movmean(slen,nmean);
    slen_mean = (slen_mean*nmean-slen)/(nmean-1); % excluding each stride itself from the moving mean
    mergeinds = find(slen<slen_mean*mergethresh);
    mergeinds(diff(mergeinds)==1) = []; % remove adjacent ones
    if isempty(mergeinds)
        fprintf('no stride meets merging criteria in crawling bout %d (%s)\n',k,behav_str{q});
        S_merged = cat(1,S_merged,S2);
        continue;
    end
    % 2) figure out whether to merge into previous or next stride
    mergeinds = cat(1,mergeinds,NaN(size(mergeinds)));
    for j=1:size(mergeinds,2)
        jj = mergeinds(1,j);
        if q==1 % forward: stride starts at vfwd peak
            if diff(S2(jj).vfwd([1 end]))<-range(S2(jj).vfwd)*diffthresh % end is not a real peak, merge into next
                mergeinds(2,j) = mergeinds(1,j)+1;
            elseif diff(S2(jj).vfwd([1 end]))>range(S2(jj).vfwd)*diffthresh % start is not a real peak, merge into previous
                mergeinds(2,j) = mergeinds(1,j)-1;
            end
        elseif q==2 % backward: stride starts at vfwd trough
            if diff(S2(jj).vfwd([1 end]))>range(S2(jj).vfwd)*diffthresh % end is not a real trough, merge into next
                mergeinds(2,j) = mergeinds(1,j)+1;
            elseif diff(S2(jj).vfwd([1 end]))<-range(S2(jj).vfwd)*diffthresh % start is not a real trough, merge into previous
                mergeinds(2,j) = mergeinds(1,j)-1;
            end
        end
    end
    mergeinds(:,isnan(mergeinds(2,:))) = []; % if not sure, then don't merge yet
    if size(mergeinds,2)>1
        mergeinds(:,any(mergeinds<1,2)) = []; % don't merge first stride into previous
        mergeinds(:,any(mergeinds>nvalidstrides,2)) = []; % don't merge last stride into next
    else
        if any(mergeinds<1 | mergeinds>nvalidstrides) % same as above
            mergeinds = [];
        end
    end
    if isempty(mergeinds)
        fprintf('no stride meets merging criteria in crawling bout %d (%s)\n',k,behav_str{q});
        S_merged = cat(1,S_merged,S2);
        continue;
    end
    mergeinds = sort(mergeinds,1);
%     stridenums = [S2.stride_num];
    stridenums = valid_strides;
    fprintf('merging %d pairs of strides in crawling bout %d (%s): ',size(mergeinds,2),k,behav_str{q});
    if size(mergeinds,2)>1
        fprintf('%d-%d,',stridenums(mergeinds(:,1:end-1)));
    end
    fprintf('%d-%d\n',stridenums(mergeinds(:,end)));
    % 3) merge strides
    % add merged strides to end of S2, then clear out to-merge strides
    for j=1:size(mergeinds,2)
        s1 = S2(mergeinds(1,j));
        s2 = S2(mergeinds(2,j));
        s = cell2struct(emptyc,fields);
        s.txinds = [s1.txinds s2.txinds];
        s.t = [s1.t s2.t];
        s.vfwd = [s1.vfwd s2.vfwd];
        s.tt = linspace(s.t(1),s.t(end),nphasebins);
        s.vfwd_tt = interp1(s.t,s.vfwd,s.tt);
        s.behav_label = s1.behav_label;
        s.stride_num = mean([s1.stride_num s2.stride_num]);
        s.crawl_num = s1.crawl_num;
        S2 = cat(1,S2,s);
    end
    S2(mergeinds) = [];
    [~,sortorder] = sort([S2.stride_num]);
    S2 = S2(sortorder);
    S_merged = cat(1,S_merged,S2);
    mergepairs = cat(2,mergepairs,reshape(stridenums(mergeinds),size(mergeinds)));
end

S_merged = rmfield(S_merged,'stride_num');

end