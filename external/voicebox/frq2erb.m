function [erb,bnd] = frq2erb(frq)
%FRQ2ERB  Convert Hertz to ERB frequency scale ERB=(FRQ)
%	erb = frq2erb(frq) converts a vector of frequencies (in Hz)
%	to the corresponding values on the ERB-rate scale on which
%  	the human ear has roughly constant resolution as judged by
%  	psychophysical measurements of the cochlear filters.

%	We have df/de = 6.23*f^2 + 93.39*f + 28.52
%	where the above expression gives the Equivalent Rectangular
%	Bandwidth (ERB)in Hz  of a human auditory filter with a centre
%	frequency of f kHz.
%
%	By integrating the reciprocal of the above expression, we
%	get:
%		e = a ln((f/p-1)/(f/q-1))
%
%	where p and q are the roots of the equation: -0.312 and -14.7
%  	and a = 1000/(6.23*(p-q)) = 11.17268
%
%	We actually implement e as
%
%		e = a ln (1 + b*f/(f+c))
%
%	where b = q/p - 1 = 46.06538
%	      c = -1000q = 14678.49
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
%      Version: $Id: frq2erb.m 713 2011-10-16 14:45:43Z dmb $
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
g=abs(frq);
erb=11.17268*sign(frq).*log(1+46.06538*g./(g+14678.49));
bnd=6.23e-6*g.^2 + 93.39e-3*g + 28.52;
if ~nargout
    plot(frq,erb,'-x');
    xlabel(['Frequency (' xticksi 'Hz)']);
    ylabel(['Frequency (' yticksi 'Erb-rate)']);
end
