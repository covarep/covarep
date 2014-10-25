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
%
% Output
%  All the outputs are uniformely sampled (see option uniform_step)
%  uf0s : [Nx2] [times, f0] 
%  uAE  : [N x amp_order] Amplitude (see option amp_enc_method)
%         (as harm amplitudes, envelope bins)
%  uAEV : Amplitude variance
%  uPE  : [N x pdm_order] Phase Distortion (PD) or Relative Phase Shift (RPS)
%         (see pd_method)
%         There is one value per harmonic (harmonic scale)
%  uPEV : Phase variance based on the phase measurement given by opt pd_method.
%         With frequency scale

function [f0s, AE, PEM, PEV, opt] = hmpd_features_compute(frames, fs, opt)

    % Parameters
    if nargin<3
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

        % Misc
        opt.usemex    = false; % Use interp1ordered TODO to false
        opt.debug     = false;
    end
    if nargin==0; f0s=opt; return; end

    if opt.amp_enc_method==1; opt.amp_order=opt.dftlen/2; end

    % Take from the frames
    f0s = [[frames.t]', [frames.f0]'];

    disp('    Estimate amplitude envelope ...');
    [AE, frames] = hmpd_amplitude_envelope_estimate(frames, fs, opt);

    disp('    Estimate phase envelope ...');
    rpspdopt = opt;
    rpspdopt.pd_method = 1;
    rpspdopt.pd_vtf_rm = false; % Already done in hmpd_amplitude_envelope_estimate
    rpspdopt.harm2freq = true;
    rpspdopt.usemex = opt.usemex;
    PE = phase_rpspd(frames, fs, rpspdopt);

    disp('    Compute phase statistics ...');
    % Phase's mean
    PEM = hmpd_phase_mean(PE, opt.pdm_nbper*opt.sin_nbat); % Below 3, doubles the voice, why?

    % Phase's standard-deviation
    % Remove first the trend from PE, otherwise the std measure is overestimated
    PEtrend = hmpd_phase_smooth(PE, opt.pdd_nbper*opt.sin_nbat); % Compute the trend
    PE = PE - PEtrend;

    PEV = hmpd_phase_deviation(PE, opt.pdd_nbper*opt.sin_nbat);

return
