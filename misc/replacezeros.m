% Replace zero values in a vector with arbitrary small values
%
% Description
%  This is a convenient function to avoid methods crashing because
%  of zero values in a signal.
% 
% Input
%  x    : A vector containing possibly zero values.
%
% Output
%  y    : The same vector as x with zero values replaced.
%
% References
%  https://github.com/covarep/covarep/issues/75
%
% Copyright (c) 2011 University of Crete - Computer Science Department
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
%  Gilles Degottex <gad27@cam.ac.uk>


function y = replacezeros(x)

    y = x;
    idx = find(x==0);
    y(idx) = eps*rand(1,length(idx));
