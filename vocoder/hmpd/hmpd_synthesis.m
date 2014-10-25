% Harmonic Model + Phase Distortion (HMPD)
%
% Copyright (c) 2012 University of Crete - Computer Science Department
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

function [syn, opt] = hmpd_synthesis(f0s, AE, PEM, PED, fs, wavlen, opt)

    if nargin<7
        % Options
        opt.enc = hmpd_features_compute();
        opt.amp_minphasegain = 1; % 1:min, 0:zero-phase, -1:max
        opt.pd_gain = 1;
        opt.rps_randunvoiced = false; % Randomize the RPS in unvoiced segments
                                      % Unvoiced segments are indicated by
                                      % f0=0 values.
                                      % (not used by default with HMPD)
        opt.pdv_corrthresh = 0.75;

        opt.defdbamp = -300; % [dB]
        opt.usemex   = false; % Use interp1ordered TODO to false
        opt.debug    = false;
    end
    if nargin==0; syn=opt; return; end
    if nargin==1;
        syn = hmpd_synthesis();
        syn.enc = f0s;
        return;
    end

    disp('HMPD Vocoder: Synthesis ...');
    if opt.debug; disp(opt); end

    % ==========================================================================

    if nargin<5 || isempty(fs);     fs=44100; end
    if nargin<6 || isempty(wavlen); wavlen=round(fs*f0s(end,1)); end

    unvoiced = f0s(:,2)==0; % Remember the unvoiced indices ...
    f0s = fillf0(f0s);      % ... and replace them with interpolated f0

    [AH, APH] = hmpd_amplitude_decompress(AE, f0s, fs, opt); % Amp decoding
    Hmax = size(AH,2)-1; % Max number of harmonics (without DC)

    % Reconstruct a Relative Phase Shift (RPS) of the voice source
    RPS = hmpd_phase_decompress(PEM, PED, f0s, fs, opt); % Phase decoding

    % If asked, force full RPS randomization in unvoiced segments
    if opt.rps_randunvoiced
        disp('    Randomize unvoiced segments (where f0=0)');
        idx = find(unvoiced);
        RPS(idx,:) = wrap(2.*pi.*rand(length(idx),size(RPS,2)));
    end

    % Add the VTF phase
    RPS = wrap(opt.pd_gain*RPS + opt.amp_minphasegain*APH); % TODO the wrap is useless

    % TODO with APH only, the waveform seems weird, not just delayed diracs

    disp(['    HM synthesis (' num2str(Hmax) ' harmonics)']);
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
