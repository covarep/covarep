% Convert dB to linear value
%
% Input
%  x      : [dB] decibel value to convert
%
% Output
%  y      : decibel value on a linear scale
%
% Example
%  y = db2lin(-6)
%
% See also
%  mag2db
%
% Copyright (c) 2011 University of Crete - Computer Science Department
% 
% License
%  This file is part of libphoni. libphoni is free software: you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. libphoni is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details (see the COPYING* files in the base directory).
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%
% $Id: db2lin.m 134 2012-07-31 15:12:33Z degottex $

function y = db2mag(x)

    if strcmp(class(x),'single');   e = realmin('single');
    else                            e = realmin;            end

	y = (10.^(x/20))-e;

return
