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

function [f0s, AE, PDM, PDD, opt] = hmpd_analysis(wav, fs, f0s, opt)
    if nargin<4
        opt = hmpd_features_compute();

        opt.uniform_step = 0.005; % [s] step size for uniform features resampling      

        % Options for the f0 estimation (if not provided in argument)
        opt.f0min   = 60;
        opt.f0max   = 440;
        opt.f0refine = false;

        % Options for the harmonic analysis
        opt.highpass = []; % High-pass the signal before sinusoidal analysis

        opt.sin = sin_analysis();
        opt.sin.use_f0time = true;
        opt.sin.harmonic   = true;
        opt.sin.use_ls     = true;
        opt.sin.fadapted   = true;
        opt.sin.debug      = false;

        opt.sin_nbat       = 4; % 4 analysis per period for statistics estimates

    end
    if nargin==0; f0s=opt; return; end
    if nargin<3; f0s=[]; end


    disp('HMPD Vocoder: Analysis ...');
    
    % Get fundamental frequency (f0)
    if ~isempty(f0s)
        % If an f0 curve is provided ...
        f0s = f0s(:,1:2); % Keep only the necessary columns
        f0s = fillf0(f0s); % Ensure there is no zero-values
    else
        % If no f0 curve is provided ...
        disp('    No f0 provided in arguments, estimate it using SRH ...');
        times = (0:length(wav)-1)/fs;

        [f0s,VUVDecisions,SRHVal,f0times] = pitch_srh(wav, fs, opt.f0min, opt.f0max, opt.uniform_step*1000);
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


    disp('    Harmonic analysis ...');

    if ~isempty(opt.highpass)
        [b, a] = ellip(6, .5, 60, 2*opt.highpass/fs, 'high');
        wav = filtfilt(b, a, wav);
    end

    % Generate analysis time instants for the harmonic analysis
    tmargin = 0.5*opt.sin.win_durnbper/max(opt.f0min,0.66*min(f0s(1,2),f0s(end,2)));
    hrts = gen_analysis_times(wav, fs, tmargin, true, 2, f0s, opt.sin_nbat);
    f0shr = interp1td(f0s, hrts);
    frames = sin_analysis(wav, fs, f0shr, opt.sin);


    disp('    Compute features ...');
    [irregf0s, AE, PDM, PDD] = hmpd_features_compute(frames, fs, opt);

    disp('    Compress features ...');
    [AE, PDM, PDD] = hmpd_features_compress(irregf0s, AE, PDM, PDD, fs, opt);


    disp('    Resample features with uniform sampling times ...');
    uT = (0:opt.uniform_step:irregf0s(end,1))'; % New uniform time instants

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
