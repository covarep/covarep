% COVAREP parallel feature formant extraction function
%
% Description
%  This function performs COVAREP_feature_formant_extraction
%  on all .wav files in the specified folder in parallel
%
% Input
%  in_dir [directory path] : Path to directory containing wav files to be analysed
%  options:     Struct with optional arguments (applied to all files in folder).
%
% Output
%  None
%
% 'options' can have the following fields:
%  feature_fs:  Sampling rate of the returned features (default: 0.01ms)
%  gender:      Used for vowelSpace (0: female (default); 1: male; 2: child)
%  features:    The list of features to be extracted (default: all).
%               Possible features are: 'f0', 'F', 'VUV', 'NAQ', 'QOQ',
%               'H1H2', 'PSP', 'MDQ', 'HRF', 'peakSlope', 'Rd', 'Rd_conf',
%               'MCEP', 'HMPDM', 'HMPDD', 'vowelSpace', 'VAD'
%  mcep_order:  Order of the MCEPs (default: 24)
%  hmpdm_order:	Order of the HMPDMs (default: 24)
%  hmpdd_order:	Order of the HMPDDs (default: 12)
%  channel:     If the audio file has multiple channels (default: 1)
%  save_mat:    Save the results as a MAT file (default: false)
%  save_csv:    Save the results as a CSV file (default: true)
%                   If 'save_mat' or 'save_csv' are strings, they are
%                   interpreted as paths to the file to be created.
% Author
%   Florian Helmhold <Florian@Helmhold.de> based on
%   COVAREP_feature_formant_extraction_perfile (Torsten Wörtwein <twoertwe@cs.cmu.edu>)
%   COVAREP_feature_extraction.m (John Kane <kanejo@.tcd.ie>)
%   COVAREP_formant_extraction.m (Stefan Scherer <scherer@ict.usc.edu>)

function COVAREP_feature_formant_extraction(in_dir, options)

%% Default settings
if nargin < 2
    options = struct; % Default empty options struct
end

%% Sanity check options
if isfield(options, 'save_mat') && ~islogical(options.save_mat)
    warning('Non boolean save_mat options will overwrite the output file.')
end

if isfield(options, 'save_csv') && ~islogical(options.save_csv)
    warning('Non boolean save_csv options will overwrite the output file.')
end

% Analysis settings
fileList=dir([in_dir filesep '*.wav']);
N=length(fileList);

if N==0
    disp('No wav files in inputted directory!!!')
end

if(hasParallelComputing)
    parfor n=1:N
        basename=regexp(fileList(n).name,'\.wav','split');
        basename=char(basename(1));
        disp(['Analysing file: ' basename])
        try
            COVAREP_feature_formant_extraction_perfile([in_dir filesep basename '.wav'], options)
            disp([basename ' successfully analysed'])
        catch err
            warning(['An error occurred while analysing ' basename ': ' getReport(err)])
        end
    end
else
    for n=1:N
        basename=regexp(fileList(n).name,'\.wav','split');
        basename=char(basename(1));
        disp(['Analysing file: ' basename])
        try
            COVAREP_feature_formant_extraction_perfile([in_dir filesep basename '.wav'], options)
            disp([basename ' successfully analysed'])
        catch err
            warning(['An error occurred while analysing ' basename ': ' getReport(err)])
        end
    end
end

%% Helper functions

function has_parallel = hasParallelComputing
% get all installed toolbox names
v = ver;
% collect the names in a cell array
[installedToolboxes{1:length(v)}] = deal(v.Name);

% check
has_parallel = all(ismember({'Parallel Computing Toolbox'},installedToolboxes));
