function [t,s]=zerocros(x,m)
%ZEROCROS finds the zeros crossings in a signal [T,S]=(X,M)% find zero crossings in a signal
% Inputs:  x = input waveform
%          m = mode string containing:
%              'p' - positive crossings only
%              'n' - negative crossings only
%              'b' - both (default)
%              'r' - round to integer values
%
% Outputs: t = sample positions of zero crossings (not necessarily integers)
%          s = estimated slope of x at the zero crossing
%
% This routine uses linear interpolation to estimate the position of a zero crossing
% A zero crossing occurs between x(n) and x(n+1) iff (x(n)>=0) ~= (x(n+1)>=0)

% Example: x=sin(2*pi*(0:1000)/200); x(1:100:1001)=0; zerocros(x);
% Note that we get a zero crossing at the end but not at the start.

%	   Copyright (C) Mike Brookes 2003
%      Version: $Id: zerocros.m 713 2011-10-16 14:45:43Z dmb $
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

if nargin<2
    m='b';
end
s=x>=0;
k=s(2:end)-s(1:end-1);
if any(m=='p')
    f=find(k>0);
elseif any(m=='n')
    f=find(k<0);
else
    f=find(k~=0);
end
s=x(f+1)-x(f);
t=f-x(f)./s;
if ~nargout
    n=length(x);
    plot(1:n,x,'-',t,zeros(length(t),1),'o');
end
