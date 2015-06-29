% HMPD: Estimate the short-term deviation of Phase Distortion
%
% Inputs
%  PD   : [NxM rad] A matrix of Phase Distortion to measure the deviation from.
%         N is the number of frames, M is the order of the PD (either the
%         maximum number of harmonics or the number of bins).
%  nbat : The number of frames to consider in the smoothing window, i.e. the
%         window size.
%  
% Outputs
%  PDD  : [NxM rad] The Phase Distortion Deviation (PDD)
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

    PDD = zeros(size(PDc));
    idx = find(z<1); % To avoid neg sqrt or sqrt(-log(0))
    PDD(idx) = sqrt(-2*log(z(idx))); % Fisher's standard-deviation

return
