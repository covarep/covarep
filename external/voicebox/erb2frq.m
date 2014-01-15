function [frq,bnd] = erb2frq(erb)
%ERB2FRQ  Convert ERB frequency scale to Hertz FRQ=(ERB)
%	frq = erb2frq(erb) converts a vector of ERB-rate values
%	to the corresponding frequencies in Hz.
%   [frq,bnd] =  erb2frq(erb) also calculates the ERB bandwidths
%
%       Note that erb values must not exceed 42.79
%
% See also: frq2erb

%	We have df/de = 6.23*f^2 + 93.39*f + 28.52
%	where the above expression gives the Equivalent Rectangular
%	Bandwidth (ERB)in Hz  of a human auditory filter with a centre
%	frequency of f kHz.
%
%	By integrating the reciprocal of the above expression, we
%	get:
%		e = a ln((f/p-1)/(f/q-1))/c
%
%	where p and q are the roots of the equation: -0.312 and -14.7
%  	and c = (6.23*(p-q))/1000 = 0.08950404
%
%	from this we can derive:
%
%	f = a/(b-exp(c*e)) - d
%
%	where a = 1000 q (1 - q/p) = 676170.4
%	      b = q/p = 47.06538
%	      d = -1000q = 14678.49
%	and f is in Hz
%
%	References:
%
%	  [1] B.C.J.Moore & B.R.Glasberg "Suggested formula for
%		calculating auditory-filter bandwidth and excitation
%		patterns", J Acoust Soc America V74, pp 750-753, 1983
%
%	  [2] O. Ghitza, "Auditory Models & Human Performance in Tasks
%		related to Speech Coding & Speech Recognition",
%		IEEE Trans on Speech & Audio Processing, Vol 2,
%		pp 115-132, Jan 1994
%	

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: erb2frq.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frq = sign(erb).*(676170.4*(47.06538-exp(0.08950404*abs(erb))).^(-1) - 14678.49);
bnd=6.23e-6*frq.^2 + 93.39e-3*abs(frq) + 28.52;
if ~nargout
    plot(erb,frq,'-x');
    xlabel(['Frequency (' xticksi 'Erb-rate)']);
    ylabel(['Frequency (' yticksi 'Hz)']);
end