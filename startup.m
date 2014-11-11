% This file should be run before running the HOWTO scripts

addpath(genpath(pwd));

% Remove all .git entries from the path ...
entries = regexp(path, ['[^',pathsep,']+'],'match');
idx = find(~cellfun(@isempty,regexp(entries,'\.git')) | ~cellfun(@isempty,regexp(entries,'\.svn')));
if ~isempty(idx); rmpath(entries{idx}); end

% ... s well as the documentation folder
rmpath(genpath([pwd '/documentation']));

if ~exist('OCTAVE_VERSION')
    % Matlab is running
    % Thus remove also the octave helper functions
    rmpath(genpath([pwd '/external/octave_only']));
else
    % Octave is running
    more off;
    pkg load tsa;
    pkg load signal;
end
