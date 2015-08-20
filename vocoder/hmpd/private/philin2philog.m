% Compress harmonic phase on a log-harmonic scale
%
% Without the vocal tract filter phase, the harmonic phases are only dependent
% on the glottal pulse shape. Because the pulse's shape is scaled with
% respect to f0, a harmonic phase value keep its harmonic number whatever the f0
% changes (strictly from a signal processing point of view, regardless the
% impact of the physiology wrt f0).
% Therefore, a log-harmonic scale (a scale which is compressed as the harmonic
% number grows) seems more appropriate for the compression of the Phase
% Distortion Mean (PDM), than a log-hertz scale (e.g. as in MFCC).
% 
% This function makes such a compression.
%
% Inputs
%  philin : The phase values to compress on a log-harmonic scale.
%  Hb     : Below this harmonic limit, the scale is linear
%           (similar to the mel scale which is linear below 1000Hz)
%           (e.g. 12)
%           Based on observation of the LF model, the asymptotic behavior of the
%           spectrum starts around the 12th harmonic (for the most tense voice)
%  Hmax   : The maxmimum number of harmonic considered during synthesis
%           (e.g. 256)
%  order  : The reduced number of phase coefficients (e.g. 24)
%
% Outputs
%  philog : The compressed coefficients
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

function philog = philin2philog(philin, Hb, Hmax, order)
    if nargin<3
        % Hb contains directly the whole log-harmonnic scale (pre-computed)
        hsl = Hb;
        Hmax = length(Hb);
        order = hsl(end);
    else
        hsl = hlin2hlog(Hb, Hmax, order);
    end
    hsl = hsl(1:min(Hmax,length(philin)));

    hsl = round(hsl);
    philog = nan(1,order);
    for h=1:order
        idx = hsl==h;
        philog(h) = angle(mean(exp(1i*philin(idx))));
    end

    % There is maybe not enough values to reach the Hmax on the linear scale
    % Thus, put random phase where values are missing
    idx = find(isnan(philog));
    philog(idx) = wrap((2*pi)*rand(length(idx),1));
    
return
