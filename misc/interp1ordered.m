% Fast linear interpolation of ordered values (one dimendionnal)
%
% This matlab version is made to illustrate the behavior of the C version.
% Do not expect speed from this Matlab version  !
% Also, if this function is updated, its C version has to be updated
% accordingly.
%
% Note tha the C version makes absolutely NO check at all.
% Thus, it can crash very easily. 
%
% Build the C version using
%  $ mex interp1ordered.c
%
% Input
%  x   ; [Nx1] Ordered reference abscissa
%  y   ; [Nx1] Ordered reference ordinates
%  xi  : [Mx1] Ordered abscissa where to interpolate
%  yid : [Mx1] Default out of bounds value
%
% Output
%  yi  : Interpolated output values
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

function yi = interp1ordered(x, y, xi, yid)

    yi = yid*ones(size(xi));

    i=1;
    j=1;

    while xi(i)<x(1)
        yi(i) = yid;
        i=i+1;
    end

    while i<=length(xi)
        while j<length(x) && x(j+1)<xi(i);
            j=j+1;
        end
        if j==length(x); break; end

        g = (xi(i)-x(j))./(x(j+1)-x(j));

        yi(i) = y(j,:)*(1-g) + y(j+1,:)*g;

        i=i+1;
    end

    if i<=length(xi)
        yi(i:end) = yid;
    end

return
