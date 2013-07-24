function mc=rotqc2mc(qc)
%ROTQC2MC converts a matrix of complex quaternion vectors to quaternion matrices
% Inputs: 
%
%     QC(2m,n)   mxn matrix of real quaternion vectors (each 2x1)
%
% Outputs: 
%
%     MC(2m,2n)   mxn matrix of real quaternion matrices (each 2x2)
%
% In matrix form, quaternions can be multiplied and added using normal matrix 
% arithmetic. Each element of an mxn matrix of quaternions is itself a 2x2 block
% so the total dimension of MC is 2m x 2n.

% 
%      Copyright (C) Mike Brookes 2000-2006
%      Version: $Id: rotqc2mc.m,v 1.3 2007/11/23 18:47:45 dmb Exp $
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
[m,n]=size(qc);
n2=2*n;
mc=zeros(m,2*n);
mc(:,1:2:n2)=qc;
mc(1:2:m,2:2:n2)=-conj(qc(2:2:m,:));
mc(2:2:m,2:2:n2)=conj(qc(1:2:m,:));
