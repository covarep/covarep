% Features computation for the Harmonic Model + Phase Distortion (HMPD)
%
% This the main entry point for estimating the parameters of the HMPD vocoder.
%
% Please read the README.txt file for general information about HMPD before
% using it.
%
% Inputs
%  frames : [Nxstruct] N structures containing the sinusoidal parameters,
%           as provided by hmpd_analysis_harmonic.m or another sinusoidal
%           analysis.
%  fs     : [Hz] The sampling rate of the analyzed waveform
%  [opt]  : Additional options (see code below)
%           Even though this argument is optional, it is highly recommended to
%           adapt the f0min and f0max options if no f0s argument is provided.
% 
% Outputs
%  f0s   : [s, Hz] [Nx2] A time/data column vector, as above, containing the
%          analysis instants and the f0 curve. 
%  AE    : [NxD] A matrix containing the amplitude envelope.
%          D is either opt.dftlen/2+1 (from hmpd_features_compute.m), or
%          opt.amp_order+1, depending if compression is disabled or enabled.
%  PDM   : [NxD] A matrix containing the Phase Distortion Mean.
%          D is either opt.dftlen/2+1 or opt.pdm_order+1, depending if
%          compression is disabled or enabled.
%  PDD   : [NxD] A matrix containing the Phase Distortion Deviation.
%          D is either opt.dftlen/2+1 or opt.pdd_order+1, depending if
%          compression is disabled or enabled.
%  opt   : The options used for the analysis.
%
% Copyright (c) 2013 University of Crete - Computer Science Department (UOC-CSD)
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
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function [f0s, AE, PDM, PDD, opt] = hmpd_analysis_features(frames, fs, opt)
    if nargin<3
        opt.sin_nbat = 4; % TODO

        opt.dftlen    = 1024;     % Used DFT length (e.g. uncompressed features)

        % Amplitude
        opt.amp_enc_method = 1; % 1:envelope; % 2:cepstral
        opt.amp_f0norm = true;
        opt.amp_log   = false;
        opt.amp_order = opt.dftlen/2; % 24mfcc smells for 32kHz, 32 good
        opt.amp_logfn = @frq2mel;

        opt.pd_vtf_rm = true; % Remove the VTF phase from the phase measurement

        % Phase
        opt.dc_phase   = 0;  % keep ori. phase if empty; otherwise, set to value
        opt.polarity_inv = false; % (applied after dc_phase is set)

        opt.pdm_nbper  = 6;  % TODO
        opt.pdm_log    = false; % If true, compress the phase coefficients using
                                % a log scale on a harmonic scale
        opt.pdm_log_hb = 8; % Below this harmonic limit, the scale is linear
                            % Then it is logarithmic (similar to the mel scale)
        opt.pdm_order  = opt.dftlen/2;

        opt.pdd_nbper = 2; % TODO
        opt.pdd_log    = false; % lin&log scale (kind of a mel scale)
        opt.pdd_order  = opt.dftlen/2; %32 for log TODO


        opt.regularstepsize = 0.005; % [s] step size for features resampling

        % Misc
        opt.usemex    = false; % Use interp1ordered TODO to false
        opt.debug     = false;
    end
    if nargin==0; f0s=opt; return; end
    if opt.amp_enc_method==1; opt.amp_order=opt.dftlen/2; end

    if opt.debug>0; disp('HMPD Vocoder: Features computation ...'); end

    irregf0s = [[frames.t]', [frames.f0]']; % Retrieve f0 from the frames

    if opt.debug>0; disp('    Estimate amplitude envelope ...'); end
    [AE, frames] = hmpd_amplitude_envelope_estimate(frames, fs, opt);

    if opt.debug>0; disp('    Estimate phase envelope ...'); end
    rpspdopt = opt;
    rpspdopt.pd_method = 1;
    rpspdopt.pd_vtf_rm = false; % Already done in hmpd_amplitude_envelope_estimate
    rpspdopt.harm2freq = true;
    rpspdopt.usemex = opt.usemex;
    PE = phase_rpspd(frames, fs, rpspdopt);

    if opt.debug>0; disp('    Compute phase statistics ...'); end

    % Phase's mean
    PDM = hmpd_phase_mean(PE, opt.pdm_nbper*opt.sin_nbat); % Below 3, doubles the voice, why?

    % Phase's standard-deviation
    % Remove first the trend from PE, otherwise the std measure is overestimated
    PEtrend = hmpd_phase_smooth(PE, opt.pdd_nbper*opt.sin_nbat); % Compute the trend
    PE = PE - PEtrend;

    PDD = hmpd_phase_deviation(PE, opt.pdd_nbper*opt.sin_nbat);


    if opt.debug>0; disp('    Compress features ...'); end
    [AE, PDM, PDD] = hmpd_features_compress(irregf0s, AE, PDM, PDD, fs, opt);


    if opt.debug>0; disp('    Resample features with uniform sampling ...'); end
    uT = (0:opt.regularstepsize:irregf0s(end,1))'; % New uniform time instants

    % First, resample the f0 curve using uniform step size
    uf0s = exp(interp1(irregf0s(:,1), log(irregf0s(:,2)), uT, 'linear', NaN));

    % Drop the time instants outside of the features
    idx = find(~isnan(uf0s));
    uT = uT(idx);
    f0s = [uT, uf0s(idx)];

    % Then, resample amplitude, PD's mean and PD's deviation
    AE = irregsampling2uniformsampling(irregf0s(:,1), AE, uT, [], [], 'linear', NaN, 1);
    PDM = irregsampling2uniformsampling(irregf0s(:,1), PDM, uT, @unwrap, @wrap, 'linear', 0, 1);
    PDD = irregsampling2uniformsampling(irregf0s(:,1), PDD, uT, [], [], 'linear', 0, 1);

return
