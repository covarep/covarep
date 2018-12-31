% Psychoacoustic scaling function from [1]
% 
% Description
% Based of the experiment in [1]_ wherein participants adjusted a
% second tone until it was half the pitch of the first. The functional
% approximation to the scale is implemented with the formula from
% [2]_ (being honest, from Wikipedia):
% 
% ..math:: s = 1127 \ln \left(1 + \frac{f}{700} \right)
% 
% Where `s` is the scale and `f` is the frequency in Hertz.
% 
% References
% [1] S. S. Stevens, J. Volkmann & E. B. Newman (1937). A Scale for
%    the Measurement of the Psychological Magnitude Pitch. The Journal
%    of the Acoustical Society of America, 8, 185-190.
% [2] O'Shaughnessy, D. (1987). Speech communication: human and
%    machine. Addison-Wesley Pub. Co.
%
%
% Copyright (c) 2018 Department of Computer Science,
%                    University of Toronto, Canada,
%                    Vector Institute, Canada
%
% License
% This file is under the LGPL license,  you can
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
%  Yingxue Wang <yingxue@cs.toronto.edu>
%  Sean Robertson <sdrobert@cs.toronto.edu>
%

classdef MelScaling < ScalingFunction
    methods (Static)
        function hertz = scale_to_hertz(scale)
            hertz = 700 * (exp(scale / 1127) - 1);
        end
        
        function scale = hertz_to_scale(hertz)
            scale = 1127 * log(1 + hertz / 700);
        end
    end
end
