% Test file for Scaling Functions
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

classdef test_scales < matlab.unittest.TestCase
    
    properties
        scaling_function
    end 
    
    methods(TestMethodSetup)
        function setScalingFunction(testCase)
            testCase.scaling_function = MelScaling();
        end
    end

    methods (Test)
        function test_scales_invertible(testCase)
             for hertz = 200:100:2000
                scale = testCase.scaling_function.hertz_to_scale(hertz);
                testCase.verifyEqual(hertz, testCase.scaling_function.scale_to_hertz(scale),'relTol',sqrt(eps));
             end
        end
    end
end
