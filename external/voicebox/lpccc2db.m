function db=lpccc2db(cc,np,nc)
%LPCCC2DBf Convert complex cepstrum to dB power spectrum DB=(CC,NP,NC)
% cc = complex cepstrum coefs
% np+2 = number of frequency values 0 to nyquist (default: np=size(cc,2))
% nc = number of cepstral coefs to use (default: nc=np);
% for high speed make np one less than a power of 2


%      Copyright (C) Mike Brookes 1998
%      Version: $Id: lpccc2db.m,v 1.4 2007/05/04 07:01:38 dmb Exp $
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
k=10/log(10);
if nc==mc
   db=k*rsfft([zeros(nf,1) cc].',2*np+2).';
else
   db=k*rsfft([zeros(nf,1) lpccc2cc(cc,nc)].',2*np+2).';
end
if nargout==0
   plot((0:np+1)/(2*np+2),db.');
   ylabel('dB');
end





