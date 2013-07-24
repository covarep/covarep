% Wrap a phase value in [-pi,pi]
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
% $Id$

function phase = wrap(phase)

    phase = phase - round(phase/2/pi)*2*pi;

    if phase>pi;        phase=phase-2*pi;
    elseif phase<-pi;   phase=phase+2*pi; end

return
