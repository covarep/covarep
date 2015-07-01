% HMPD: Compute a smooth estimate of the Phase Distortion (a local trend)
%  
%  Very similar to the PDM computation.
%  However, the window size is different and this one makes use of median and
%  hanning filter in order to avoid outliers and Gibbs phenomenon.
%  Unfortunately, this particular filtering was forgot in the publications [1,2]
%  
% Inputs
%  PD   : [NxM rad] A matrix of Phase Distortion to be smoothed
%         N is the number of frames, M is the order of the PD (either the
%         maximum number of harmonics or the number of bins).
%  nbat : The number of frames to consider in the smoothing window, i.e. the
%         window size.
%  
% Outputs
%  PD   : [NxM rad] The smoothed Phase Distortion
%
% Copyright (c) 2013 University of Crete - Computer Science Department(UOC-CSD)/ 
%                    Foundation for Research and Technology-Hellas - Institute
%                    of Computer Science (FORTH-ICS)
%
% References
%  [1] G. Degottex and D. Erro, "A uniform phase representation for the harmonic
%      model in speech synthesis applications", EURASIP, Journal on Audio, Speech,
%      and Music Processing - Special Issue: Models of Speech - In Search of Better
%      Representations, 2014.
%  [2] G. Degottex and D. Erro, "A Measure of Phase Randomness for the Harmonic
%      Model in Speech Synthesis", In Proc. Interspeech, Singapore. International
%      Speech Communication Association (ISCA), September 2014.
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

% Smooth the PD, which is similar to computing its local trend
function PD = hmpd_phase_smooth(PD, nbat)

    winlen = round(nbat/2)*2+1;
    win = hann(winlen); % ERRATUM: Was not mentioned in the papers [1,2]
    win = win./sum(win);

    % Smooth the values in polar coordinates
    PDc = cos(PD);
    PDs = sin(PD);

    % Apply first a median filter:
    % Because: i) The glottal pulse is not supposed to drastically change from
    %             one frame to the next.
    %         ii) It keeps the outliers in the residual phase, which optimize
    %             the variance measurement.
    % ERRATUM: Was not mentioned in the original papers [1,2]
    PDc = medfilt1(PDc, winlen);
    PDs = medfilt1(PDs, winlen);
%      PDc = medfilt1(PDc, winlen, [], 1);
%      PDs = medfilt1(PDs, winlen, [], 1);

    % Then smooth the steps of the median filter
    PDc = filtfilt(win, 1, PDc);
    PDs = filtfilt(win, 1, PDs);

    % Move back to radians
    PD = atan2(PDs, PDc);

return
