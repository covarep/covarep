% Create a colormap which is consistent for circular values
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

function map = circmap()

    % Current "most pleasant"
    reds   = interp1([-pi -2*pi/3 0 2*pi/3 pi], [1      0     0    1     1], (-pi:2*pi/64:pi-2*pi/64));
    greens = interp1([-pi -2*pi/3 0 2*pi/3 pi], [1      1     0    0     1], (-pi:2*pi/64:pi-2*pi/64));
    blues  = interp1([-pi -2*pi/3 0 2*pi/3 pi], [0      0     1    0     0], (-pi:2*pi/64:pi-2*pi/64));

    %  reds   = interp1([-pi -2*pi/3 0 2*pi/3 pi], [0      1     0    1     0], (-pi:2*pi/64:pi-2*pi/64));
    %  greens = interp1([-pi -2*pi/3 0 2*pi/3 pi], [1      1     0    0     1], (-pi:2*pi/64:pi-2*pi/64));
    %  blues  = interp1([-pi -2*pi/3 0 2*pi/3 pi], [0      0     1    0     0], (-pi:2*pi/64:pi-2*pi/64));

    % Awful, put the purpule at pi
    %  reds   = interp1([-pi -2*pi/3 0 2*pi/3 pi], [1      0     0    1     1], (-pi:2*pi/64:pi-2*pi/64));
    %  greens = interp1([-pi -2*pi/3 0 2*pi/3 pi], [0      0     1    0     0], (-pi:2*pi/64:pi-2*pi/64));
    %  blues  = interp1([-pi -2*pi/3 0 2*pi/3 pi], [1      1     0    0     1], (-pi:2*pi/64:pi-2*pi/64));

    map = [reds', greens', blues'];

    %  cyan = 0 0.75 1

return
