% COVAREP feature formant extraction function
%
% Description
%  This function can extract most of COVAREP's features. It
%  processes an audio file and returns the results as a MATLAB table.
%  Options to save the result as MAT file or CSV exist. You can specify a
%  subset of the features to be extracted to shorten the possibly long
%  processing time (only the necessary computations are done).
% 
% Input
%  audio_file:  Path to the audio file.
%  options:     Struct with optional arguments.
%
% Output
%  MATLAB table with a column for 'time' and each feature. 'vowelSpace' is
%  duplicated to fit in the table.
%
% 'options' can have the following fields:
%  feature_fs:  Sampling rate of the returned features (default: 0.01ms)
%  gender:      Used for vowelSpace (0: female (default); 1: male; 2: child)
%  features:    The list of features to be extracted (default: all). 
%               Possible features are: 'f0', 'F', 'VUV', 'NAQ', 'QOQ', 
%               'H1H2', 'PSP', 'MDQ', 'HRF', 'peakSlope', 'Rd', 'Rd_conf',
%               'MCEP', 'HMPDM', 'HMPDD', 'vowelSpace', 'VAD'
%  start:       Extract features after this point in time (in seconds;
%               default: 0).
%  stop:        Extract features before this point in time (in seconds). If
%               negative, no stop limit is enforced (default: -1)
%  mcep_order:  Order of the MCEPs (default: 24)
%  hmpdm_order:	Order of the HMPDMs (default: 24)
%  hmpdd_order:	Order of the HMPDDs (default: 12)
%  channel:     If the audio file has multiple channels (default: 1)
%  save_mat:    Save the results as a MAT file (default: false)
%  save_csv:    Save the results as a CSV file (default: true)
%                   If 'save_mat' or 'save_csv' are strings, they are
%                   interpreted as paths to the file to be created.
%
% Author
%  Torsten Wörtwein <twoertwe@cs.cmu.edu> based on
%   COVAREP_feature_extraction.m (John Kane <kanejo@.tcd.ie>) and
%   COVAREP_formant_extraction.m (Stefan Scherer <scherer@ict.usc.edu>)

function results = COVAREP_feature_formant_extraction_perfile(audio_file, options)

all_features = {'f0', 'VUV', 'NAQ', 'QOQ', 'H1H2', 'PSP', 'MDQ', 'HRF', ...
    'peakSlope','Rd', 'Rd_conf','MCEP', 'HMPDM', 'HMPDD', ...
    'vowelSpace', 'F', 'VAD'};

%% input arguments
if ~exist('options', 'var'), options = struct(); end
if ~isfield(options, 'feature_fs'), options.feature_fs = 0.01; end
if ~isfield(options, 'gender'), options.gender = 0; end
if ~isfield(options, 'features'), options.features = all_features; end
if ~isfield(options, 'start'), options.start = 0; end
if ~isfield(options, 'stop'), options.stop = -1; end
if ~isfield(options, 'mcep_order'), options.mcep_order = 24; end
if ~isfield(options, 'hmpdm_order'), options.hmpdm_order = 24; end
if ~isfield(options, 'hmpdd_order'), options.hmpdd_order = 12; end
if ~isfield(options, 'channel'), options.channel = 1; end
if ~isfield(options, 'save_mat'), options.save_mat = false; end
if ~isfield(options, 'save_csv'), options.save_csv = true; end

% F0 settings
if ~isfield(options, 'F0min'), options.F0min = 50; end
if ~isfield(options, 'F0max'), options.F0max = 500;
% IAIF settings
if ~isfield(options, 'hpfilt'), options.hpfilt = 1; end
if ~isfield(options, 'd'), options.d = 0.99; end
% LP settings
if ~isfield(options, 'LP_winLen'),  options.LP_winLen = 0.025; end
if ~isfield(options, 'LP_winShift'),  options.LP_winShift = 0.005; end
% Envelope/Rd MSP settings
opt = sin_analysis();
opt.fharmonic = true;
opt.use_ls = false; % Use Peak Picking
opt.dftlen = 4096;  % Force the DFT length
opt.frames_keepspec = true; % Keep the computed spectra in the frames structure
opt.debug = 0;
if isfield(options, 'fharmonic'),  opt.fharmonic = options.fharmonic; end
if isfield(options, 'use_ls'),  opt.use_ls = options.use_ls; end
if isfield(options, 'dftlen'),  opt.dftlen = options.dftlen; end
if isfield(options, 'frames_keepspec'),  opt.dftlen = options.frames_keepspec; end

assert(any(options.gender == [0, 1, 2]), 'Not a valid gender argument');
assert(all(ismember(options.features, all_features)), ...
    ['Possible features are: ' strjoin(all_features, ' ')])

%% Load file and select channel
[x, fs] = audioread(audio_file);
if ismatrix(x), x = x(:, options.channel); end

%% Trim audio to start/stop
start_index = round(fs * options.start + 1);
options.start = (start_index - 1) / fs;
stop_index = round(fs * options.stop);
if options.stop < 0
    stop_index = numel(x);
end
stop_index = min(stop_index, numel(x));
x = x(start_index:stop_index);

%% Polarity detection
x = polarity_reskew(x, fs) * x; % Correct polarity if necessary

%% pitch_srh: required for f0, VUV, and all frame-based features
if any(ismember(options.features, {'f0', 'VUV', 'vowelSpace', 'Rd', ...
        'MCEP', 'HMPDM', 'HMPDD', 'NAQ', 'QOQ', 'H1H2', 'PSP', 'MDQ', 'HRF'}))
    [srh_f0, srh_vuv, ~, srh_time] = pitch_srh(x, fs, options.F0min, ...
        options.F0max, options.feature_fs*1000);
    srh_f0(srh_f0 <= options.F0min) = options.F0min;
end

%% gci_sedreams and inverse filtering: required for all inverse-filtered features
if any(ismember(options.features, {'NAQ', 'QOQ', 'H1H2', 'PSP', 'MDQ', 'HRF'}))
    % gci_sedreams
    F0med = median(srh_f0(srh_f0>options.F0min & srh_f0<options.F0max & srh_vuv==1));
    GCI = gci_sedreams(x, fs, F0med, 1); % SEDREAMS GCI detection
    GCI = round(GCI*fs);
    GCI(GCI<1 | isnan(GCI)==1 | isinf(GCI)==1) = [];
    VUV_int = interp1(round(srh_time*fs), srh_vuv, 1:length(x));
    VUV_int(isnan(VUV_int))=0;
    GCI(VUV_int(GCI)<.5) = []; % Remove GCIs in detected unvoiced regions
    GCI = unique(GCI); % Remove possible duplications
    
    % Iterative and adaptive inverse filtering (IAIF) & LP inverse filtering
    p_gl = 2 * round(fs/4000);
    p_vt = 2 * round(fs/2000) + 4;
    [g_iaif, gd_iaif] = iaif_gci(x, fs, GCI/fs, p_vt, p_gl, options.d, options.hpfilt);
    res = lpcresidual(x, options.LP_winLen*fs, options.LP_winShift*fs, fs/1000+2); % LP residual
end

%% Glottal source parameterisation
if any(ismember(options.features, {'NAQ', 'QOQ', 'H1H2', 'PSP', 'HRF'}))
    % Estimate conventional glottal parameters
    [NAQ, QOQ, H1H2, HRF, PSP] = get_vq_params(g_iaif, gd_iaif, fs, GCI/fs);
end

%% Maxima dispersion quotient measurement
if any(strcmp(options.features, 'MDQ'))
    MDQ = mdq(res, fs, GCI/fs); 
end

%% Peak Slope
if any(strcmp(options.features, 'peakSlope'))
    PS = peakslope(x, fs);
end

%% sin_analysis: required for all frame-based features
if any(ismember(options.features, {'Rd', 'MCEP', 'HMPDM', 'HMPDD'}))
    frames = sin_analysis(x, fs, [srh_time(:), srh_f0(:)], opt);
end

%% Rd parameter estimation of the LF glottal model using Mean Squared Phase (MSP)
if any(strcmp(options.features, 'Rd'))
    rds = rd_msp(frames, fs);
end

%% MCEP
if any(strcmp(options.features, 'MCEP'))
    % Spectral envelope parameterisation
    M = numel(frames);
    MCEP = zeros(M, options.mcep_order+1);
    TE_orders = round(0.5*fs./[frames.f0]); % optimal cepstral order
    spec = hspec2spec(vertcat(frames.S));
    TE_orders_unique = unique(TE_orders);
    for m=1:numel(TE_orders_unique)
        idx = TE_orders_unique(m)==TE_orders;
        MCEP(idx,:) = hspec2fwcep(env_te(spec(idx,:), TE_orders_unique(m))', ...
            fs, options.mcep_order)';
    end
end

%% HMPDM and HMPDD
if any(ismember(options.features, {'HMPDM', 'HMPDD'}))
    hmpdopt = hmpd_analysis();
    hmpdopt.debug = 0;
    hmpdopt.usemex = false;
    hmpdopt.amp_enc_method = 2;
    hmpdopt.amp_log = true;
    hmpdopt.pdd_log = true; 
    hmpdopt.pdm_log = true; 
    hmpdopt.amp_order = 39;
    hmpdopt.pdd_order = options.hmpdd_order;
    hmpdopt.pdm_order = options.hmpdm_order;
    [hmpdf0s, ~, HMPDM, HMPDD] = hmpd_analysis_features(frames, fs, hmpdopt);
end

%% Voice activation detection
if any(strcmp(options.features, 'VAD'))
    [VAD, ~, ~, ~, VAD_time] = VAD_Drugman(x, fs, false);
end

%% Formants
if any(ismember(options.features, {'F', 'vowelSpace'}))
    formantPeaks = formant_CGDZP(x, fs, 30, options.feature_fs*1000);
end

%% vowel space
if any(ismember(options.features, {'vowelSpace'}))
    formant_time = linspace(1, length(x), size(formantPeaks,1));
    VUV_int = interp1(round(srh_time*fs), srh_vuv, formant_time) > 0.5;
    vowelSpace = getVowelSpace(formantPeaks(VUV_int, 1:2), options.gender);
end

%% Create Table and sample features
feature_sampling = (options.feature_fs/2*fs):(options.feature_fs*fs):length(x);
results = array2table(feature_sampling'/fs, 'VariableNames', {'time'});
for ifeature = 1:numel(options.features)
    feature = options.features{ifeature};
    if any(startsWith(results.Properties.VariableNames, feature)), continue; end
    
    % sampling
    names = {};
    signal = [];
    special_sampling = false;
    
    switch (feature)
        case {'f0', 'VUV'}
            time = round(srh_time*fs);
            if any(strcmp('f0', options.features))
                signal = srh_f0';
                names = {'f0'};
            end
            if any(strcmp('VUV', options.features))
                signal = [signal srh_vuv'];
                names = [names 'VUV'];
            end
        case {'HMPDM', 'HMPDD'}
            special_sampling = true;
            time = hmpdf0s(:,1);
            feature_sampling_ = (feature_sampling-1) / fs;
            if any(strcmp('HMPDM', options.features))
                signal = irregsampling2uniformsampling(time, HMPDM, feature_sampling_, @unwrap, @wrap, 'linear', 0, hmpdopt.usemex);
                names = arrayfun(@(x) strcat('HMPDM_', num2str(x)), 0:options.hmpdm_order, 'UniformOutput', false);
            end
            if any(strcmp('HMPDD', options.features))
                signal = [signal irregsampling2uniformsampling(time, HMPDD, feature_sampling_, [], [], 'linear', 0, hmpdopt.usemex)];
                names = [names arrayfun(@(x) strcat('HMPDD_', num2str(x)), 0:options.hmpdd_order, 'UniformOutput', false)];
            end
        case {'NAQ', 'QOQ', 'H1H2', 'PSP', 'HRF'}
            time = NAQ(:, 1)*fs;
            if any(strcmp('NAQ', options.features))
                signal = NAQ(:, 2);
                names = {'NAQ'};
            end
            if any(strcmp('QOQ', options.features))
                signal = [signal QOQ(:, 2)];
                names = [names 'QOQ'];
            end
            if any(strcmp('H1H2', options.features))
                signal = [signal H1H2(:, 2)];
                names = [names 'H1H2'];
            end
            if any(strcmp('PSP', options.features))
                signal = [signal PSP(:, 2)];
                names = [names 'PSP'];
            end
            if any(strcmp('HRF', options.features))
                signal = [signal HRF(:, 2)];
                names = [names 'HRF'];
            end
        case {'Rd'}
            names = {'Rd', 'Rd_conf'};
            time = rds(:, 1)*fs;
            signal = [rds(:, 2) rds(:, 3)];
        case {'MDQ'}
            names = {feature};
            time = MDQ(:, 1)*fs;
            signal = MDQ(:, 2);
        case {'peakSlope'}
            names = {feature};
            time = PS(:, 1)*fs;
            signal = PS(:, 2);
        case {'MCEP'}
            names = arrayfun(@(x) strcat('MCEP_', num2str(x)), 0:options.mcep_order, 'UniformOutput', false);
            time = linspace(1, length(x), size(MCEP,1));
            signal = MCEP;
        case {'VAD'}
            names = {feature};
            signal = VAD';
            time = VAD_time*fs;
        case {'F'}
            names = arrayfun(@(x) strcat('F', num2str(x)), 1:size(formantPeaks,2), 'UniformOutput', false);
            time = linspace(1, length(x), size(formantPeaks,1));
            signal = formantPeaks;
        case {'vowelSpace'}
            special_sampling = true;
            names = {feature};
            signal = repmat(vowelSpace, [size(feature_sampling, 2) 1]);
    end
    
    assert(size(signal, 2) == numel(names))
    if ~special_sampling
        signal = interp1(time, signal, feature_sampling);
    end
    if size(signal, 1) == numel(names)
        signal = signal';
    end
    signal(isnan(signal)) = 0;
    
    % binarize VUV after interpolation
    vuv_index = strcmp('VUV', names);
    if any(vuv_index)
        signal(:, vuv_index) = signal(:, vuv_index) > 0.5;
    end
    
    % add to table
    results = [results array2table(signal, 'VariableNames', names)];
end

%% adjust time to account for different start point
results.time = results.time + options.start;

%% save results
[dirs, name, ~] = fileparts(audio_file);
if isempty(dirs), dirs = '.'; end
mat_file = strcat(dirs, filesep, name, '.mat');
csv_file = strcat(dirs, filesep, name, '.csv');
if ischar(options.save_mat) || options.save_mat
    if ischar(options.save_mat)
        mat_file = options.save_mat;
    end
    save(mat_file, 'results');
end
if ischar(options.save_csv) || options.save_csv
    if ischar(options.save_csv)
        csv_file = options.save_csv;
    end
    writetable(results, csv_file);
end
end
