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

% Compress phase coefficents using hlin2hlog.m
function philog = philin2philog(philin, Hb, Hmax, order)

    hsl = hlin2hlog(Hb, Hmax, order);
    hsl = hsl(1:min(Hmax,length(philin)));

    hsl = round(hsl);
    philog = nan(1,order);
    for h=1:order
        idx = find(hsl==h);
        philog(h) = angle(mean(exp(1i*philin(idx))));
    end

    % There is maybe not enough values to reach the Hmax on the linear scale
    % Thus, put random phase where values are missing
    idx = find(isnan(philog));
    philog(idx) = wrap((2*pi)*rand(length(idx),1));
    
return
