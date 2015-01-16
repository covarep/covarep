% Replace zeros f0 values by linear interpolation
%
% Copyright (c) 2012 University of Crete - Computer Science Department (UOC-CSD)
%
% License
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function f0sout = fillf0(f0sin)

    f0sout = f0sin;

    idx = find(f0sin(:,2)~=0);
    if length(idx)==1
        f0sout(:,2) = f0sin(idx,2);
    elseif length(idx)>1
        f0sout(:,2) = 2.^(interp1(f0sin(idx,1), log2(f0sin(idx,2)), f0sin(:,1), 'linear', NaN));
    end

    f0sout(:,2) = interp1_fixnan(f0sout(:,2));

return
