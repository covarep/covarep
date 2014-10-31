% Create a spectrum with a linear phase, a given delay in samples
%
% Description
%  If the spectrum length is even, the Nyquist's phase is managed as follow:
%   sign(real(exp((delay*1i*pi))))
%
% Inputs
%  delay    : [samples]
%  fftlen   : length of the fft
%
% Outputs
%  shift    : the delay-spectrum
%
% Example
%  S = fft(...);
%  D = delay2spec(12, 1024);
%  S = S.*D;  % shift the time signal by 12 samples to the left
%
% Copyright (c) 2008 Ircam/CNRS-UMR9912-STMS
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
%  Gilles Degottex <gilles.degottex@ircam.fr>
%

function shift = delay2spec(delay, dftlen)
    if delay==0
        shift = ones(1,dftlen);
    else
        % odd length
        if mod(dftlen,2)==1
            shift = exp((delay*2i*pi./dftlen).*(1:(dftlen-1)/2));
            shift = [1, shift(1:end), conj(shift(end:-1:1))];

        % even length
        else
            shift = exp((delay*2i*pi./dftlen).*(1:dftlen/2));
            shift = [1, shift(1:end-1), sign(real(shift(end))), conj(shift(end-1:-1:1))];
        end
    end
return

