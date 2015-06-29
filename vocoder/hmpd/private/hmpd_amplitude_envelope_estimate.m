% HMPD: Estimate amplitude envelopes given sinusoidal parameters for all frames
% 
% Inputs
%  frames : [Nxstruct] N structures containing the sinusoidal parameters,
%           as provided by sin_analysis.m
%  fs     : [Hz] The sampling rate of the analyzed waveform
%  opt    : Additional options (see hmpd_features_compute.m)
% 
% Outputs
%  AE    : [NxD] A matrix containing the amplitude envelope.
%          D is either opt.dftlen/2+1 (from hmpd_features_compute.m), or
%          opt.amp_order+1, depending if compression is disabled or enabled.
%  frames: The same frames as in the input arguments with some possible
%          modification.
%          For example, if opt.pd_vtf_rm is True, the minimum-phase response of
%          the estimated amplitude envelope is removed from the phase parameters
%  opt   : The used options.
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

function [AE, frames, opt] = hmpd_amplitude_envelope_estimate(frames, fs, opt)

    if opt.usemex
        phiinterp1fn = @(x, y, xi, yid) interp1ordered(x.', y.', xi.', yid);
        ampinterp1fn = @(x, y, xi, varargin) interp1ordered(x.', y.', xi.', NaN);
    else
        phiinterp1fn = @(x, y, xi, yid) interp1(x, y, xi, 'linear', yid);
        ampinterp1fn = @(x, y, xi) interp1(x, y, xi, 'linear', NaN);
    end

    AE = zeros(numel(frames),opt.dftlen/2+1); % Pre-allocate the necessary space
    F = fs*(0:opt.dftlen/2)/opt.dftlen; % The bins' frequency of the envelope

    if opt.debug>0; pg=progressbar(numel(frames)); end
    for n=1:numel(frames)

        % Estimate an amplitude envelope
        E = env_interp(frames(n).sins(1:2,2:end), fs, opt.dftlen, true, 'interp1fn', ampinterp1fn);
        E(isnan(E)) = db2mag(opt.amp_def);

        if opt.amp_f0norm; E = E./sqrt(frames(n).f0); end

        AE(n,:) = abs(E);

        % If asked, remove the min-phase of the envelope
        if opt.pd_vtf_rm
            lE = hspec2minphaseloghspec(E);
            vtfp = phiinterp1fn(F, imag(lE), frames(n).sins(1,:), 0);
            frames(n).sins(3,:) = wrap(frames(n).sins(3,:) - vtfp(:)');
            frames(n).sins = [frames(n).sins; vtfp(:)']; % Keep the "vtf" phase
        end

        if opt.debug>0; pg=progressbar(pg, n); end
    end 
    
return
