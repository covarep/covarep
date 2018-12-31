% A unit-normalized Hann window
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

classdef HannWindow < WindowFunction

    methods (Static)
        function window = get_impulse_response(obj, width)
            % Impulse response of a unit-normalized Hann window
            % 
            % Inputs
            % width : [int] width of the window
            % 
            % Outputs
            % window : returns the N-point symmetric Hann window in a 
            %           column vector.
            %
            % Example
            % >> width = 64;
            % >> window = HannWindow.get_impulse_response(width);
            % 
            
            window = hann(width);
            window = window / 0.5 * max(1, width - 1);
        end
    end
end