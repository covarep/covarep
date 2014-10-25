% Convert a bandwidth to a root radius
%
% Input
%  bw     : [Hz] bandwidth to convert
%  fs     : [Hz] sampling frequency
%
% Output
%  radius : corresponding decay rate
%
% Example
%  radius = bandwidth2rootradius(50, 44100)
%
% Copyright (c) 2007 Ircam-CNRS UMR9912-STMS
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
%  Gilles Degottex <degottex@ircam.fr>
%

function radius = bandwidth2rootradius(bw, fs)

	radius = exp((-pi/fs).*bw);

return
