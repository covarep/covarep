% Replace NaN values at boundaries of a vector by constant extrapolation
% 
% Octave compatible
% 
% Copyright (c) 2012 University of Crete - Computer Science Department
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

function x = interp1_extrapbounds(x)

    % Check if bounds are defined, and replace nan values
    if isnan(x(1))
        firstnonnan = min(find(~isnan(x)));
        x(1:firstnonnan-1) = x(firstnonnan);
    end
    if isnan(x(end))
        lastnonnan = max(find(~isnan(x)));
        x(lastnonnan+1:end) = x(lastnonnan);
    end

return

