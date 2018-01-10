% This file should be run before running the HOWTO scripts
function startup

current_directory = fileparts(mfilename('fullpath'));
entries = strsplit(genpath(current_directory), pathsep);

% helper function to find matching entries
entry_matcher = @(entries, pattern) ~cellfun(@isempty, strfind(entries, pattern));

% Remove all .git and .svn entries from the path
entries(entry_matcher(entries, '.git') | entry_matcher(entries, '.svn')) = [];

% as well as the documentation folder
entries(entry_matcher(entries, strcat(filesep, 'documentation'))) = [];

if ~exist('OCTAVE_VERSION', 'builtin')
    % Matlab is running
    % Thus remove also the octave helper functions
    octave_path = strcat(filesep, 'external', filesep, 'octave_only');
    entries(entry_matcher(entries, octave_path)) = [];
else
    % Octave is running
    more off;
    pkg load tsa;
    pkg load signal;
end

% remove backward compatible functions
MATLABVERSION=version('-release');
MATLABVERSION=MATLABVERSION(1:end-1);
if str2double(MATLABVERSION)>=2015
    entries(entry_matcher(entries, 'backcompatibility_2015')) = [];
end

% addpaths
addpath(strjoin(entries, pathsep));
end