function s = estetime(seconds,formatSpec)
%ESTETIME convert a time period from seconds to largest appropriate time unit

seconds = double(seconds);
if nargin==1
    formatSpec = '%0.3g';
end
if strcmp(formatSpec,'%d')
    formatSpec = '%.0g';
end

if seconds>=1
    timeunits = {' hr',' min',' sec'};
    t = mod(seconds,[0 3600 60])./[3600 60 1];
    [~,ind,~] = find(t>=1,1);
    s = [sprintf(formatSpec,t(ind)) timeunits{ind}];
else
    s = [sprintf(formatSpec,seconds) ' sec'];
end

end