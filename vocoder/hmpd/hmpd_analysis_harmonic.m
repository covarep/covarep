% Harmonic model analysis for the Harmonic Model + Phase Distortion (HMPD)
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
%  TODO
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

function [frames, opt] = hmpd_analysis_harmonic(wav, fs, f0s, opt)
    if nargin<4
        % Options for the f0 estimation
        % (used only if f0 not provided as argument)
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
        opt.sin.debug      = false;% Do not print anything while running

        opt.sin_nbat       = 4; % 4 analysis per period for statistics estimates
    end
    if nargin==0; frames=opt; return; end
    if nargin<3; f0s=[]; end

    if opt.debug>0; disp('HMPD Vocoder: Harmonic analysis ...'); end
    
    % Get fundamental frequency (f0)
    if ~isempty(f0s)
        % If an f0 curve is provided ...
        f0s = f0s(:,1:2); % Keep only the necessary columns
        f0s = fillf0(f0s); % Ensure there is no zero-values
    else
        % If no f0 curve is provided ...
        if opt.debug>0; disp('    No f0 provided in arguments, estimate it using SRH ...'); end
        times = (0:length(wav)-1)/fs;

        [f0s,VUVDecisions,SRHVal,f0times] = pitch_srh(wav, fs, opt.f0min, opt.f0max, 5);
        f0s = [f0times(:), f0s(:)];
        f0s(find(~VUVDecisions),2) = 0;
        f0s = fillf0(f0s);

        if opt.f0refine
            airopt = ahm_air_analysis2();
            airopt.do_final_ahm_step = false;
            f0s = ahm_air_analysis2(wav, fs, f0s, airopt);
        end

        if 0
            times = (0:length(wav)-1)'/fs;
            plot(times, wav, 'k');
            hold on;
            plot(f0s(:,1), log2(f0s(:,2))-log2(440), 'r');
            title('Waveform');
            xlabel('Time [s]');
            keyboard
        end
    end


    if opt.debug>0; disp('    Estimation of sinusoidal parameters ...'); end

    if ~isempty(opt.highpass)
        [b, a] = ellip(6, .5, 60, 2*opt.highpass/fs, 'high');
        wav = filtfilt(b, a, wav);
    end

    % Generate analysis time instants for the harmonic analysis
    tmargin = 0.5*opt.sin.win_durnbper/max(opt.f0min,0.66*min(f0s(1,2),f0s(end,2)));
    hrts = gen_analysis_times(wav, fs, tmargin, true, 2, f0s, opt.sin_nbat);
    f0shr = interp1td(f0s, hrts);

    frames = sin_analysis(wav, fs, f0shr, opt.sin);

return
