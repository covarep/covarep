% Harmonic Model + Phase Distortion (HMPD)
%
% Copyright (c) 2013 University of Crete - Computer Science Department(UOC-CSD)/ 
%                    Foundation for Research and Technology-Hellas - Institute
%                    of Computer Science (FORTH-ICS)
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

% Compute mean phase
function PDM = hmpd_phase_mean(PD, nbat)

    winlen = round(nbat/2)*2+1;
    % win = hann(winlen);
    win = ones(winlen,1); % Rectangular window seems better than smooth window
                          % To study ... TODO
    win = win./sum(win);

    % Compute the mean in polar coordinates
    % Compute the center of gravity of the exp(i*angles)
    PDc = filtfilt(win, 1, cos(PD));
    PDs = filtfilt(win, 1, sin(PD));

    % Return the angle of the center of gravity, i.e. the mean value
    PDM = atan2(PDs, PDc);

return
