function cc=lpcar2cc(ar,np)
%LPCAR2CC LPC: Convert ar filter to complex cepstrum CC=(AR,NP)
% the "real" cepstrum is half the complex cepstrum
% cc() does not include c0 whose value can be calculated
% from the prediction residual energy, e, as ln(e)/2
% for both real and complex cepstrum.



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: lpcar2cc.m,v 1.4 2007/05/04 07:01:38 dmb Exp $
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

[nf,p1]=size(ar);
p=p1-1;
if (nargin<2) np=p; end
cc=zeros(nf,np);
cm=(1:np).^(-1);
if np>p
  xm=-(1:p);
  nz=np-p;
  for k=1:nf
    cc(k,:)=filter(1,ar(k,:),[ar(k,2:p1).*xm zeros(1,nz)]).*cm;
  end
else
  p1=np+1;
  xm=-(1:np);
  for k=1:nf
    cc(k,:)=filter(1,ar(k,:),ar(k,2:p1).*xm).*cm;
  end
end
