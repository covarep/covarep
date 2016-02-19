% General COVAREP feature extraction script
%
% Description
%  This script extracts features to do with glottal source and spectral
%  envelope available with the COVAREP repository. For each .wav file in
%  the inputted directory path, .mat files are produced containing 
% 
% Input
%  in_dir [directory path] : Path to directory containing wav files to be analysed
%  sample_rate [seconds] : feature sampling rate in seconds (optional)
%
% Output
%         : No arguments are outputted with this script, though .mat files
%         are saved corresponding to each .wav file. .mat files contain the
%         feature matrix: features [number of frames X 35] and names:
%         containing the feature name correspond to each column of the
%         feature matrix
% Example
%   in_dir='//home/john/Desktop/test_dir'; % Specify directory of wavs
%   sample_rate=0.01; % State feature sampling rate
%   COVAREP_feature_extraction(in_dir,sample_rate); % Launch feature extraction
%
% Copyright (c) 2013 Trinity College Dublin - Phonetics & Speech Lab
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
%  John Kane <kanejo@.tcd.ie>


function COVAREP_feature_extraction(in_dir,sample_rate)

%% Initial settings
if nargin < 2
    sample_rate=0.01; % Default to 10 ms sampling
end

% F0 settings
F0min = 50; % Minimum F0 set to 50 Hz
F0max = 500; % Maximum F0 set to 500 Hz

% IAIF settings
hpfilt = 1;
d = 0.99;

% LP settings
LP_winLen=0.025;
LP_winShift=0.005;

% Rd MSP settings
opt = sin_analysis();
opt.fharmonic  = true;
opt.use_ls     = false;
opt.debug = 0;

% Envelope settings
opt.use_ls     = false; % Use Peak Picking
opt.dftlen     = 4096;  % Force the DFT length
opt.frames_keepspec = true; % Keep the computed spectra in the frames structure
MCEP_ord=24;

% Analysis settings
fileList=dir([in_dir filesep '*.wav']);
N=length(fileList);
names={'F0','VUV','NAQ','QOQ','H1H2','PSP','MDQ','peakSlope','Rd', ...
    'Rd_conf','creak','MCEP_0','MCEP_1','MCEP_2','MCEP_3','MCEP_4','MCEP_5', ...
    'MCEP_6','MCEP_7','MCEP_8','MCEP_9','MCEP_10','MCEP_11','MCEP_12', ...
    'MCEP_13','MCEP_14','MCEP_15','MCEP_16','MCEP_17','MCEP_18', ...
    'MCEP_19','MCEP_20','MCEP_21','MCEP_22','MCEP_23','MCEP_24',...
    'HMPDM_0','HMPDM_1','HMPDM_2','HMPDM_3','HMPDM_4','HMPDM_5', ...
    'HMPDM_6','HMPDM_7','HMPDM_8','HMPDM_9','HMPDM_10','HMPDM_11','HMPDM_12', ...
    'HMPDM_13','HMPDM_14','HMPDM_15','HMPDM_16','HMPDM_17','HMPDM_18', ...
    'HMPDM_19','HMPDM_20','HMPDM_21','HMPDM_22','HMPDM_23','HMPDM_24',...
    'HMPDD_0','HMPDD_1','HMPDD_2','HMPDD_3','HMPDD_4','HMPDD_5', ...
    'HMPDD_6','HMPDD_7','HMPDD_8','HMPDD_9','HMPDD_10','HMPDD_11','HMPDD_12'};

if N==0
    disp('No wav files in inputted directory!!!')
end

%% Do processing
for n=1:N
   
    basename=regexp(fileList(n).name,'\.wav','split');
    basename=char(basename(1));
    disp(['Analysing file: ' basename])
    try
        % Load file and set sample locations
        [x,fs]=audioread([in_dir filesep basename '.wav']);
        feature_sampling=round((sample_rate/2)*fs):round(sample_rate*fs):length(x);
        
        % Check if signal is mono or stereo
        if(size(x, 2) ~= 1)
            warning(['file: ' basename ' is not a mono signal. processing only first channel.']);
            x = x(:,1);
        end

        % Polarity detection
        polarity = polarity_reskew(x,fs);
        x=polarity*x; % Correct polarity if necessary

        % F0/GCI detection 
        [srh_f0,srh_vuv,~,srh_time] = pitch_srh(x,fs,F0min,F0max, ...
            sample_rate*1000);
        F0med=median(srh_f0(srh_f0>F0min&srh_f0<F0max&srh_vuv==1));
        F0 = interp1(round(srh_time*fs),srh_f0,feature_sampling);
        VUV = interp1(round(srh_time*fs),srh_vuv,feature_sampling);
        VUV_int = interp1(round(srh_time*fs),srh_vuv,1:length(x));
        VUV(isnan(VUV)==1)=0; VUV_int(isnan(VUV_int)==1)=0; 
        VUV(VUV>=.5)=1; VUV(VUV<.5)=0;

        GCI = gci_sedreams(x,fs,F0med,1); % SEDREAMS GCI detection
        GCI=round(GCI*fs); GCI(GCI<1|isnan(GCI)==1|isinf(GCI)==1)=[];
        GCI(VUV_int(GCI)<.5)=[]; % Remove GCIs in detected unvoiced regions
        GCI=unique(GCI); % Remove possible duplications

        % Iterative and adaptive inverse filtering (IAIF) & LP inverse
        % filtering
        p_gl = 2*round(fs/4000);
        p_vt = 2*round(fs/2000)+4;
        [g_iaif,gd_iaif] = iaif_gci(x,fs,GCI/fs,p_vt,p_gl,d,hpfilt);
        res = lpcresidual(x,LP_winLen*fs,LP_winShift*fs,fs/1000+2); % LP residual

        % Glottal source parameterisation
        [NAQ,QOQ,H1H2,HRF,PSP] = get_vq_params(g_iaif,gd_iaif,fs,GCI/fs); % Estimate conventional glottal parameters

        % Wavelet-based parameters
        MDQ = mdq(res,fs,GCI/fs); % Maxima dispersion quotient measurement
        PS = peakslope(x,fs);   % peakSlope extraction
        MDQ=interp1(MDQ(:,1)*fs,MDQ(:,2),feature_sampling);
        PS=interp1(PS(:,1)*fs,PS(:,2),feature_sampling);

        % Rd parameter estimation of the LF glottal model using Mean Squared
        % Phase (MSP)
        srh_f0(srh_f0==0) = 100;
        frames = sin_analysis(x, fs, [srh_time(:),srh_f0(:)], opt);
        rds = rd_msp(frames, fs);

        % Creaky voice detection
        warning off
        try
            creak_pp = detect_creaky_voice(x,fs); % Detect creaky voice
            creak_pp=interp1(creak_pp(:,2),creak_pp(:,1),feature_sampling);
        catch
            creak_pp=zeros(length(feature_sampling),1);
        end
        warning on

        % Spectral envelope parameterisation
        M=numel(frames);
        MCEP=zeros(M,MCEP_ord+1);
        TE_orders = round(0.5*fs./[frames.f0]); % optimal cepstral order
        spec = hspec2spec(vertcat(frames.S));
        TE_orders_unique = unique(TE_orders);
        for m=1:numel(TE_orders_unique)
            idx = TE_orders_unique(m)==TE_orders;
            MCEP(idx,:) = hspec2fwcep(env_te(spec(idx,:), TE_orders_unique(m))',...
                fs, MCEP_ord)';
        end

        % Interpolate features to feature sampling rate
        NAQ=interp1(NAQ(:,1)*fs,NAQ(:,2),feature_sampling);
        QOQ=interp1(QOQ(:,1)*fs,QOQ(:,2),feature_sampling);
        H1H2=interp1(H1H2(:,1)*fs,H1H2(:,2),feature_sampling);
        PSP=interp1(PSP(:,1)*fs,PSP(:,2),feature_sampling);
        Rd=interp1(rds(:,1)*fs,rds(:,2),feature_sampling);
        Rd_conf=interp1(rds(:,1)*fs,rds(:,3),feature_sampling);

        MCEP_int=zeros(length(feature_sampling),MCEP_ord+1);
        for m=1:MCEP_ord+1
            MCEP_int(:,m) = interp1(round(linspace(1,length(x),size(MCEP,1))),MCEP(:,m),feature_sampling);
        end

        % Add PDM and PDD
        hmpdopt = hmpd_analysis();
        hmpdopt.debug = 0;
        hmpdopt.usemex = false;
        hmpdopt.amp_enc_method=2; hmpdopt.amp_log=true; hmpdopt.amp_order=39;
        hmpdopt.pdd_log=true; hmpdopt.pdd_order=12;% MFCC-like phase variance
        hmpdopt.pdm_log=true; hmpdopt.pdm_order=24;% Number of log-Harmonic coefs
        [hmpdf0s, ~, HMPDM, HMPDD] = hmpd_analysis_features(frames, fs, hmpdopt);
        HMPDM = irregsampling2uniformsampling(hmpdf0s(:,1), HMPDM, (feature_sampling-1)/fs, @unwrap, @wrap, 'linear', 0, hmpdopt.usemex);
        HMPDD = irregsampling2uniformsampling(hmpdf0s(:,1), HMPDD, (feature_sampling-1)/fs, [], [], 'linear', 0, hmpdopt.usemex);

        % Create feature matrix and save
        features=[F0(:) VUV(:) NAQ(:) QOQ(:) H1H2(:) PSP(:) MDQ(:) PS(:) ...
            Rd(:) Rd_conf(:) creak_pp(:) MCEP_int HMPDM HMPDD];
        features(isnan(features))=0;
        save([in_dir filesep basename '.mat'],'features','names')
        clear features

        disp([basename ' successfully analysed'])
        
    catch err
        warning(['An error occurred while analysing ' basename ': ' getReport(err)])
    end
    
end
