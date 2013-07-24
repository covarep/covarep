function q=qrmult(q1,q2)
%QRMULT multiplies together two real quaternions q=[q1,q2]
%
% Inputs:
%
%     q1(4,1), q2(4,1)  Two real quaternions in the form [r, i, j, k]' where i^2=j^2=k^2=ijk=-1
%
% Outputs: 
%
%     q(4,1)   Product of q1 and q2

%      Copyright (C) Mike Brookes 2000-2008
%      Version: $Id: qrmult.m,v 1.1 2008/12/03 10:07:53 dmb Exp $
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
persistent a b c d
if isempty(a)
    a=[5 8 9 10 15 13];
    b=[6 7 11 12 14 16];
    c=[1 2 3 4 6 7 11 12 16 14];
    d=[1 2 3 4 5 8 9 10 13 15];
end
t=q1*q2.';
s=zeros(4,4);
s(a)=-t(b);
s(c)=t(d);
q=sum(s,2);
