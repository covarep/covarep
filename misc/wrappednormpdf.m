% 
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

function p = wrappednormpdf(x, m, s, N)

    % N=10 gives a log error smaller than 1e-12 for s=8
    % with s=8 distribution is "almost" uniform (std < 1e-12)
    if nargin<4; N=10; end

    p = 0;
    for k=-N:N
        p = p + normpdf(x+2.*pi.*k, m, s);
    end

return
