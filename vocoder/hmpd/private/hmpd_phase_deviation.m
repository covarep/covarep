% Harmonic Model + Phase Distortion (HMPD)
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

function PDD = hmpd_phase_deviation(PD, nbat)

    winlen = round(nbat/2)*2+1;
    % win = hann(winlen);
    win = ones(winlen,1); % Rectangular window is better than smooth window
                          % It better discriminates voiced/unvoiced segments
    win = win./sum(win);

    % Compute the std in polar coordinates
    % Compute first the center of gravity of the exp(i*angles)
    PDc = filtfilt(win, 1, cos(PD));
    PDs = filtfilt(win, 1, sin(PD));

    z = abs(PDc + 1j*PDs); % For the variance, we need only the magnitude

    z = min(1,z); % BUG FIX: without this: abs(z)=1 can imply PDD=0+0i

    PDD = sqrt(-2*log(z)); % Fisher's standard-deviation

return
