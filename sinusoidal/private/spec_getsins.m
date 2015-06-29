% Estimate all sinusoids of a spectrum by peak peaking
%  
% Octave compatible
% 
% Inputs
%  S         : Spectrum
%  fs        : [Hz] Sampling frequency
%  
%  The remaining arguments are forwarded to spec_fit_freq_amp.m
%
% Outputs
%  sins      : vector of partials structures [Hz; linear; rad]
%                 + the is always included !
%
% Copyright (c) 2011 University of Crete - Computer Science Department
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

function sins = spec_getsins(S, fs, varargin)

    S = S(:).';

    dftlen = length(S);

    idx = v_findpeaks(abs(S(1:floor(end/2)))); % TODO Use 'q'
    
    sins = zeros(5,length(idx));
    for k=1:length(idx)
        [sins(1,k) sins(2,k)] = spec_fit_freq_amp(S, fs, idx(k), varargin{:});
    end
    sins(5,:) = true; % they are all from spectral peaks

    % Estimate the phase from linear interpolation of the phase spectrum
    sins(3,:) = wrap(interp1q(fs*(0:dftlen/2)'/dftlen, unwrap(angle(S(1:end/2+1))).', sins(1,:)'));

    % Add the DC
    sins = [[0; abs(S(1)); angle(S(1)); 0; abs(S(1))>abs(S(2))], sins];

return

