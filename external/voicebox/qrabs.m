function [m,q]=qrabs(q1)
%QRABS absoloute value and normalization of a real quaternions [m,q]=[q1]
%
% Inputs:
%
%     q1(4,1)  A real quaternion in the form [r, i, j, k]' where i^2=j^2=k^2=ijk=-1
%
% Outputs: 
%
%     m   Magnitude of the quaternion, q1.
%     q   Normalized version of q1 with unit magnitude.

%      Copyright (C) Mike Brookes 2000-2008
%      Version: $Id: qrabs.m,v 1.1 2008/12/03 10:07:53 dmb Exp $
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
m=sqrt(q1'*q1);
if m>0
    q=q1/m;
else
    q=[1 0 0 0]';
end
