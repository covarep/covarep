% Replace NaN values by extrapolation using the nearest value
%
% Copyright (c) 2011 University of Crete - Computer Science Department (UOC-CSD)
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

function x = interp1_fixnan(x, method)

    if nargin<2; method='nearest'; end

    idx = find(~isnan(x));
    if length(idx)<length(x)
        x = interp1(idx, x(idx), (1:length(x)), method, 'extrap');
    end

return
