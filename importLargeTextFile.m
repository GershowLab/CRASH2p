function out = importLargeTextFile(filename,varargin)

% default params
delimiter = '\t';
chunksize = 1E7; % ~500MB

assignApplicable(varargin);

% open .txt file and read column headers
fid = fopen(filename);
% frewind(fid);
headerstr = fgetl(fid); % read pointer is now at next line (real data begin)
headercell = strsplit(headerstr,delimiter);
ncol = numel(headercell);

% create formatSpec string
formatstr = '%ld';
for i=1:ncol-1
    formatstr = [formatstr delimiter '%ld']; %#ok<AGROW>
end
formatstr = [formatstr '\n'];

% read and parse numeric data chunk by chunk
D = [];
while ~feof(fid)
    a = fscanf(fid,formatstr,[ncol chunksize])';
    D = cat(1,D,a);
end
D = double(D);
out.data = D;
out.textdata = headercell;
out.colheaders = headercell;

fclose(fid);

end