% Helper function for oq_egg.m, estimate Oq at regular intervals using oq_egg
%
% Octave compatible
% 
% (please see oq_egg.m for documentation and HOWTO_egg.m for a usage example)
%
%
% Copyright (c) 2013 University of Crete - Computer Science Department
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

function [Oqeggmess,f0mess, atimes] = oq_egg_sig(s, fs, f0, method)

    if nargin<3
        f0 = [];
    end

    winLen=2*round(20/1000*fs /2)+1;
    winShift=round(5/1000*fs);

    % Do processing
    start=1;
    stop=start+winLen-1;
    cnt=1;

    while stop <= length(s)

        % Get windowed frame
        x_frame = s(start:stop);
        atimes(cnt) = (0.5*(stop+start))/fs;

        [Oqeggmess(cnt), f0mess(cnt)] = oq_egg(x_frame, fs, f0, method);

        % Increment
        start=start+winShift-1;
        stop=start+winLen-1;
        cnt=cnt+1;
    end

return
