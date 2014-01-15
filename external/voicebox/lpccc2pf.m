function pf=lpccc2pf(cc,np,nc)
%LPCCC2PF Convert complex cepstrum to power spectrum PF=(CC,NP,NC)
% cc = complex cepstrum coefs
% np+2 = number of frequency values 0 to nyquist (default: np=size(cc,2))
% nc = number of cepstral coefs to use (default: nc=np);
% for high speed make np one less than a power of 2


%      Copyright (C) Mike Brookes 1998
%      Version: $Id: lpccc2pf.m 713 2011-10-16 14:45:43Z dmb $
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

[nf,mc]=size(cc);
if nargin<2 np=mc+2; end
if nargin<3 nc=np; end
if nc==mc
   pf=exp(rsfft([zeros(nf,1) cc].',2*np+2).');
else
   pf=exp(rsfft([zeros(nf,1) lpccc2cc(cc,nc)].',2*np+2).');
end




