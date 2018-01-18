% COVAREP parallel feature formant extraction function
%
% Description
%  This function performs COVAREP_feature_formant_extraction
%  on all .wav files in the specified folder in parallel
%
% Input
%  in_dir [directory path] : Path to directory containing wav files to be analysed
%  options:     Struct with optional arguments.
%
% Output
%  None
%
% Author
%   Florian Helmhold <Florian@Helmhold.de> based on
%   COVAREP_feature_formant_extraction Torsten Wörtwein <twoertwe@cs.cmu.edu>
%   COVAREP_feature_extraction.m (John Kane <kanejo@.tcd.ie>)
%   COVAREP_formant_extraction.m (Stefan Scherer <scherer@ict.usc.edu>)

function COVAREP_parallel_feature_formant_extraction(in_dir, options)

%% Initial settings
if nargin < 2
    sample_rate=0.01; % Default to 10 ms sampling
end

% Analysis settings
fileList=dir([in_dir filesep '*.wav']);
N=length(fileList);

if N==0
    disp('No wav files in inputted directory!!!')
end

parfor n=1:N
    basename=regexp(fileList(n).name,'\.wav','split');
    basename=char(basename(1));
    disp(['Analysing file: ' basename])
    try
        COVAREP_feature_formant_extraction([in_dir filesep basename '.wav'], options)
        disp([basename ' successfully analysed'])
    catch err
        warning(['An error occurred while analysing ' basename ': ' getReport(err)])
    end
end
