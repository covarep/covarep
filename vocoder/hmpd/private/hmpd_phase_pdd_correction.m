% HMPD: Correct Phase Distortion Deviation
%
% The PDD is always under-estimated, for reason which need to be investigated.
% This function does a simple engeenering solution: It amplify the PDD's 
% distribution above a given threshold.
%  
% Inputs
%  PDD  : [NxM rad] A matrix of Phase Distortion Deviation (PDD) to correct.
%         N is the number of frames, M is the order of the PDD (either the
%         maximum number of harmonics or the number of bins).
%  thr  : The PDD's threshold. Below this value, PDD is kept, above, PDD is
%         amplified.
%  
% Outputs
%  PDD  : [NxM rad] The corrected Phase Distortion Deviation (PDD)
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

function PDD = hmpd_phase_pdd_correction(PDD, thr)

    if thr==Inf; return; end

    if 1
        % Use cst above thr
        idx = find(PDD>thr);
        PDD(idx) = 2;
    elseif 0
        % Use a linear factor above thr (same effect as above)
        idx = find(PDD>thr);
        PDD(idx) = 10*(PDD(idx)-thr)+thr; % multiply pdd by 10
    end

return
