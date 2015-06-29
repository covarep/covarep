% Parameter analysis for the Harmonic Model + Phase Distortion (HMPD)
%
% Octave compatible
%
% This the main entry point for estimating the parameters of the HMPD vocoder.
%
% Please read the README.txt file for general information about HMPD before
% using it.
%
% Inputs
%  wav   : The samples of the waveform to analyse
%  fs    : [Hz] The sampling rate of the waveform
% [f0s]  : [s, Hz] [Mx2] A time/data column vector with time instants t[s] and
%          fundamental frequency f0[Hz] estimated at each time t.
%          It can be empty. In this case, it will be estimated by the SRH method
%          provided in COVAREP.
%          WARNING: Even if you do not provide your own f0 curve, it is still
%                   important to adapt the f0min and f0max values in the options
%                   according to the analysized voice!
%          It is very likely that the output features will not be estimated
%          at the time instants t. They will be given at regular time instants
%          (see option opt.regularstepsize), regardless the f0s argument.
% [opt]  : Additional options (see code below)
%          Even though this argument is optional, it is highly recommended to
%          adapt the f0min and f0max options if no f0s argument is provided.
% 
% Outputs
%  f0s   : [s, Hz] [Nx2] A time/data column vector, as above, containing the
%          used fundamental frequency curve. 
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
% Copyright (c) 2012 University of Crete - Computer Science Department (UOC-CSD)
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

function [f0s, AE, PDM, PDD, opt] = hmpd_analysis(wav, fs, f0s, opt)
    if nargin<4
        % Options for the f0 estimation
        % (used only if f0s not provided in argument)
        opt.f0min   = 60;
        opt.f0max   = 440;
        opt.f0refine = false;

        % Options for the harmonic analysis
        opt.highpass = []; % High-pass the signal before sinusoidal analysis

        opt.sin = sin_analysis();
        opt.sin.use_f0time = true; % Use the times adapted to f0 curve.
        opt.sin.fharmonic  = true; % Force strict harmonicity
        opt.sin.use_ls     = true; % By default, use Least Square (LS) solution
        opt.sin.fadapted   = true; % and adapt the frequency basis to f0 curve
        opt.sin.debug      = 0;    % Do not print anything while running

        % features computation
        opt.sin_nbat  = 4;    % Number of analysis instant per period
        opt.dftlen    = 1024; % Used DFT length (e.g. for uncompressed features)
        opt.regularstepsize = 0.005;% [s] step size for features resampling

        % Amplitude
        opt.amp_enc_method = 1; % 1:envelope; % 2:cepstral
        opt.amp_f0norm = true;
        opt.amp_log   = false;
        opt.amp_order = opt.dftlen/2; % 24mfcc smells for 32kHz, 32 good
        opt.amp_logfn = @frq2mel;
        opt.amp_def   = -300; % [dB] def value that replaces zero values

        opt.pd_vtf_rm = true; % Remove the VTF phase from the phase measurement

        % Phase
        opt.dc_phase   = 0;  % keep ori. phase if empty; otherwise, set to value
        opt.polarity_inv = false; % (applied after dc_phase is set)

        opt.pdm_nbper  = 6;  % Number of period to consider for PDM
        opt.pdm_log    = false; % If true, compress the phase coefficients using
                             % a log scale on a harmonic scale
        opt.pdm_log_hb = 8;  % Below this harmonic limit, the scale is linear
                             % Then it is logarithmic (similar to the mel scale)
        opt.pdm_order  = opt.dftlen/2;% PDM's order for compression

        opt.pdd_nbper  = 2;            % Number of period to consider for PDD
        opt.pdd_log    = false;        % lin&log scale (kind of a mel scale)
        opt.pdd_order  = opt.dftlen/2; % PDD's order for compression

        % Misc
        opt.usemex    = false; % Use mex fn, faster but use linear interpolation
        opt.debug     = 1;
    end
    if nargin==0; f0s=opt; return; end
    if nargin<3; f0s=[]; end
    if opt.amp_enc_method==1; opt.amp_order=opt.dftlen/2; end
    if exist('OCTAVE_VERSION'); opt.usemex=false; end

    % Estimate sinusoidal harmonic parameters
    frames = hmpd_analysis_harmonic(wav, fs, f0s, opt);

    % Compute amplitude envelope and phase statistics
    [f0s, AE, PDM, PDD] = hmpd_analysis_features(frames, fs, opt);

return
