function setupDirectories(username)
%SETUPDIRECTORIES add relevant paths to the MATLAB path
% if username is specified, add user specific/username too; assuming MATLAB
% is launched in TrackingMicroscopeAnalysis/

addpath(genpath(pwd));
rmpath(genpath('user specific'));
if (exist('username','var') && ~isempty(username))
    usersubdir = fullfile(pwd, 'user specific', username, '');
    if ~exist(usersubdir,'dir')
        mkdir(usersubdir);
    end
    addpath(genpath(usersubdir));
end

end