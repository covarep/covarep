addpath(genpath(pwd));

rmpath(genpath([pwd '/.git']));

rmpath(genpath([pwd '/documentation']));

if ~exist('OCTAVE_VERSION')
    % Matlab is running
    % Thus remove the ocatve helper functions
    rmpath(genpath([pwd '/external/octave_only']));
end
