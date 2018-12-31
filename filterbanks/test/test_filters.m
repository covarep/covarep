% Test file for Gamma Window Filters
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

classdef test_filters <  matlab.unittest.TestCase
    properties (TestParameter)
        window_size = {10, 100, 1000}
        peak_ratio = {0.5, 0.75, 0.9} 
        order = {2,4}
    end
    
    methods (Test)
        function test_gamma_window_peak_matches(testCase, window_size, peak_ratio, order)
            expected_max_idx = floor(window_size * peak_ratio);
            window = GammaWindow(order, peak_ratio).get_impulse_response(window_size);
            [max_window, max_idx] = max(window);
            testCase.verifyTrue(expected_max_idx == max_idx || expected_max_idx == max_idx+1);
        end
    end

end