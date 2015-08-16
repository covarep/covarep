% Specialized COVAREP feature extraction script for Formants and Vowel
% space
%
% Description
%  This script extracts formant features available with the COVAREP repository.
%  For each .wav file in the input-directory, .mat files are produced containing
%  extracted features. Function is parallelizable using 'parfor' instead of
%  simple 'for'. Function will automatically run in parallel mode if the
%  Parallel Computing Toolbox is installed.
%
% Input
%  in_dir [directory path] : Path to directory containing wav files to be analysed
%  sample_rate [seconds] : feature sampling rate in seconds (optional)
%
% Output
%         : No arguments are outputted with this script, though .mat files
%         are saved corresponding to each .wav file. .mat files contain the
%         feature matrices: formants [number of frames x 5] and vowelSpace [1x1]
%
% Example
%   in_dir='//home/test_dir'; % Specify directory of wavs
%   sample_rate=0.01; % State feature sampling rate
%   COVAREP_formant_extraction(in_dir,sample_rate); % Launch feature extraction
%
% Copyright (c) 2014 University of Southern California, Insitute for
% Creative Technologies
%
% License
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General
%  Public License as published by the Free Software Foundation, either version 3
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% Note
%  This function has been developed and tested on speech signals sampled
%  at 16 kHz. Though the analysis should be sampling frequency independent
%  we cannot guarantee optimal performance on non-16 kHz signals
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Author
%  Stefan Scherer scherer@ict.usc.edu


function COVAREP_formant_extraction(in_dir,sample_rate)

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

%% Do processing
if(hasParallelComputing)
    parfor n=1:N
        basename=regexp(fileList(n).name,'\.wav','split');
        basename=char(basename(1));
        disp(['Analysing file: ' basename])
        try
            % Load file and set sample locations
            [x,fs]=audioread([in_dir filesep basename '.wav']);
            feature_sampling=round((sample_rate/2)*fs):round(sample_rate*fs):length(x);
            
            [formantPeaks_int, vowelSpace] = do_computation(x,fs,sample_rate,feature_sampling);
            
            wrap_up(basename, in_dir, formantPeaks_int, vowelSpace);
            
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
            % Load file and set sample locations
            [x,fs]=audioread([in_dir filesep basename '.wav']);
            feature_sampling=round((sample_rate/2)*fs):round(sample_rate*fs):length(x);
            
            [formantPeaks_int, vowelSpace] = do_computation(x,fs,sample_rate,feature_sampling);
            
            wrap_up(basename, in_dir, formantPeaks_int, vowelSpace);
            
            disp([basename ' successfully analysed'])
        catch err
            warning(['An error occurred while analysing ' basename ': ' getReport(err)])
        end
    end
end

function [formantPeaks_int, vowelSpace] = do_computation(x,fs,sample_rate,feature_sampling)
% Polarity detection
polarity = polarity_reskew(x,fs);
x=polarity*x; % Correct polarity if necessary

fprintf('extracting formants\n');
tic
[formantPeaks,~]=formant_CGDZP(x,fs,30,sample_rate*1000);
formantPeaks_int=zeros(length(feature_sampling),size(formantPeaks,2));
for m=1:size(formantPeaks,2)
    formantPeaks_int(:,m) = interp1(round(linspace(1,length(x),size(formantPeaks,1))),formantPeaks(:,m),feature_sampling);
end
toc

fprintf('calculating vowel space\n');
tic
[vowelSpace,~]=getVowelSpace(formantPeaks_int(:,1:2));
toc

function wrap_up(basename, in_dir, formantPeaks_int, vowelSpace)

% Create feature matrix and save
formants = formantPeaks_int;
formants(isnan(formants))=0;
save([in_dir filesep basename '.mat'],'formants','vowelSpace')
clear formants

function has_parallel = hasParallelComputing
% get all installed toolbox names
v = ver;
% collect the names in a cell array
[installedToolboxes{1:length(v)}] = deal(v.Name);

% check
has_parallel = all(ismember({'Parallel Computing Toolbox'},installedToolboxes));

