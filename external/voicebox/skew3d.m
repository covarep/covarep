function x=skew3d(y)
%SKEW3D Convert between a 3-dimensional vector and the corresponding skew-symmetric matrix
% If A and B are 3-dimensional vectors, the vector cross product is
%                   A x B = skew3d(A)*B = -skew3d(B)*A
% skew3d() is its own inverse, i.e. skew3d(skew3d(A))=A
% see also http://www.ee.ic.ac.uk/hp/staff/dmb/matrix/special.html#skew_symmetric

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: skew3d.m,v 1.4 2007/05/04 07:01:39 dmb Exp $
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

k=prod(size(y));
if k==3
   x=zeros(3,3);
   x([6 7 2])=y(:)';
   x([8 3 4])=-y(:)';
elseif k==9
   x=y([6 7 2]');
else
   error('input must be a vector or matrix of order 3');
end