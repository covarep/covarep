% Waveform synthesis for the Harmonic Model + Phase Distortion (HMPD)
%
% Octave compatible
%
% This the main entry point for synthesizing a waveform using HMPD parameters.
%
% Please read the README.txt file for general information about HMPD before
% using it.
%
% Inputs
%  f0s    : [s, Hz] [Nx2] A time/data column vector, containing the
%           fundamental frequency f0 values. 
%  AE     : [NxD] A matrix containing the amplitude envelope.
%           D is either opt.dftlen/2+1 (from hmpd_features_compute.m), or
%           opt.amp_order+1, depending if compression is disabled or enabled.
%  PDM    : [NxD] A matrix containing the Phase Distortion Mean.
%           D is either opt.dftlen/2+1 or opt.pdm_order+1, depending if
%           compression is disabled or enabled.
%  PDD    : [NxD] A matrix containing the Phase Distortion Deviation.
%           D is either opt.dftlen/2+1 or opt.pdd_order+1, depending if
%           compression is disabled or enabled.
%  fs     : [Hz] The sampling rate of the waveform.
%           (It has to be the same as the one used during analysis)
% [wavlen]: Length of the synthesized waveform [in number of samples].
%           If this option is omited, it is predicted from the last time instant
%           of the f0s argument.
% [opt]   : Additional options for synthesis (see code below).
%           The option opt.enc has to contain absolutely the options used
%           during analysis.
%
% Outputs
%  syn   : The samples of the synthesized waveform
%  opt   : The options used during synthesis.
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

function [syn, opt] = hmpd_synthesis(f0s, AE, PDM, PDD, fs, wavlen, opt)

    if nargin<7
        % Options
        opt.enc = hmpd_analysis();
        opt.amp_minphasegain = 1; % 1:minimum-phase, 0:zero-phase, -1:max-phase
        opt.pd_gain = 1;
        opt.rps_randunvoiced = false; % Randomize the RPS in unvoiced segments
                                      % Unvoiced segments are indicated by
                                      % f0=0 values.
                                      % (not used by default with HMPD)
        opt.pdv_corrthresh = 0.75;    % PDD correction threshold

        opt.defdbamp = -300;  % [dB]
        opt.usemex   = false; % Use mex fn, faster but use linear interpolation
        opt.debug    = 1;
    end
    if nargin==0; syn=opt; return; end
    if nargin==1;
        syn = hmpd_synthesis();
        syn.enc = f0s;
        return;
    end

    if opt.debug>0; disp('HMPD Vocoder: Synthesis ...'); end
    if opt.debug>1; disp(opt); end

    % ==========================================================================

    if nargin<5 || isempty(fs);     fs=44100; end
    if nargin<6 || isempty(wavlen); wavlen=round(fs*f0s(end,1)); end

    unvoiced = f0s(:,2)==0; % Remember the unvoiced indices ...
    f0s = fillf0(f0s);      % ... and replace them with interpolated f0

    [AH, APH] = hmpd_amplitude_decompress(AE, f0s, fs, opt); % Amp decoding
    Hmax = size(AH,2)-1; % Max number of harmonics (without DC)

    % Reconstruct a Relative Phase Shift (RPS) of the voice source
    RPS = hmpd_phase_decompress(PDM, PDD, f0s, fs, opt); % Phase decoding

    % If asked, force full RPS randomization in unvoiced segments
    if opt.rps_randunvoiced
        if opt.debug>0; disp('    Randomize unvoiced segments (where f0=0)');end
        idx = find(unvoiced);
        RPS(idx,:) = wrap(2.*pi.*rand(length(idx),size(RPS,2)));
    end

    % Add the VTF phase
    RPS = wrap(opt.pd_gain*RPS + opt.amp_minphasegain*APH);

    if opt.debug>0; disp(['    HM synthesis (' num2str(Hmax) ' harmonics)']);end
    deflogamp = log(db2mag(opt.defdbamp));

    % Compute instantaneous fundamental frequency along the whole recording
    f1 = interp1(f0s(:,1), f0s(:,2), (0:wavlen-1)'/fs, 'linear');
    f1 = interp1_fixnan(f1)'; % Fix bounds
    p1 = filter(1, [1 -1], (2*pi/fs)*f1);

    syn = zeros(wavlen,1); % Allocate the synthetic signal
    % For each harmonic
    for h=1:Hmax % here h is the harmonic number

        % Synthesize only indices where the frequency doesn't go above Nyquist
        idxsyn = find(f1*h<fs/2);
        if length(idxsyn)==0; break; end
        synts = (idxsyn-1)/fs;

        % Amplitude ------------------------------------------------------------
        % Interpolate on a log scale
        Y = AH(:,1+h);
        if opt.usemex
            am = interp1ordered(f0s(:,1), Y, synts, deflogamp).';
        else
            am = interp1(f0s(:,1), Y, synts, 'linear', deflogamp);
        end

        % Phase ----------------------------------------------------------------
        % Interpolate on RPS domain
        Y = unwrap(RPS(:,1+h));
        if opt.usemex
            rpst = interp1ordered(f0s(:,1), Y, synts, 0).'; % Lose the splines !
        else
            rpst = interp1(f0s(:,1), Y, synts, 'spline');
        end
        % Instantaneous phase = linear phase + relative phase shift
        ph = h.*p1(idxsyn) + rpst;

        % Harmonic synthesis ---------------------------------------------------
        syn(idxsyn) = syn(idxsyn) + 2*exp(am).*cos(ph);
    end

return
