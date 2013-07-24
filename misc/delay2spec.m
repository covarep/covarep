% Create a spectrum with a linear phase, a given delay in samples
%
% DESCRIPTION
%
%  If the spectrum length is even, the Nyquist's phase is managed as follow:
%   sign(real(exp((delay*1i*pi))))
%
% USAGE
%  shift = delay2spec(delay, fftlen)
%
% INPUT
%  delay    : [samples]
%  fftlen   : length of the fft
%
% OUPUT
%  shift    : the delay-spectrum
%
% REQUIRED
%
% EXAMPLE
%  S = fft(...);
%  D = delay2spec(12, 1024);
%  S = S.*D;  % shift the time signal by 12 samples to the left
%
% AUTHOR
%  degottex@ircam.fr
%
% $Id$

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

