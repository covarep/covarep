% A lowpass filter based on the Gamma function
% 
% Description
% A Gamma function is defined as:
% .. math:: p(t; \alpha, n) \defeq t^{n - 1} e^{-\alpha t} u(t)
% 
% Where :math:`n` is the order of the function, :math:`\alpha`
% controls the bandwidth of the filter, and :math:`u` is the step
% function.
% 
% This function returns a window based off a reflected Gamma function.
% :math:`\alpha` is chosen such that the maximum value of the window
% aligns with `peak`. The window is clipped to the width. For
% reasonable values of `peak` (i.e. in the last quarter of samples),
% the majority of the support should lie in this interval anyways.
% 
% Properties
% order : int
% peak : int
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
    
classdef GammaWindow < WindowFunction
    
    properties
        order
        peak
    end
    
    methods

        function obj = GammaWindow(order, peak)
            % Constructor of a GammaWindow object
            % 
            % Inputs
            % order : [int]
            % peak : [double]
            %     ``peak * width``, where ``width`` is the length of the window
            %     in samples, is where the approximate maximal value of the window
            %     lies
            % 
            % Outputs
            % obj : [GammaWindow object]
            % 
            % Example
            % >> window_size = 100;
            % >> peak_ratio = 0.75
            % >> order = 2;
            % >> window = GammaWindow(order, peak_ratio).get_impulse_response(window_size);
            % 
    
            obj = obj@WindowFunction();
            
            % Assign default parameters
            if nargin == 0
                obj.order = 4;
                obj.peak = 0.75;
            elseif nargin == 1
                obj.order = order;
                obj.peak = 0.75;
            else
                obj.order = order;
                obj.peak = peak;
            end
        end
        
        function window = get_impulse_response(obj, width)
            % Validate parameters
            if nargin < 1
                error('Please specify width of the window!');
            end
            if width <= 0
                window = [];
            elseif width == 1
                window = [1];
            else
                peak = obj.peak * width;
                ret = width - 1: -1:0; 
                if obj.order > 1
                    alpha = (obj.order - 1) / (width - peak);
                    offs = width - 1;
                else
                    % align alpha roughly with a support in M
                    alpha = 5 / width;
                    offs = width;
                end
                ln_c = obj.order * log(alpha);
                ln_c = ln_c - log(factorial(obj.order - 1));
                ret(1:offs) = ret(1:offs) .^ (obj.order - 1) .* exp(-alpha * ret(1:offs) + ln_c);
                window = ret;
            end
        end
    end
end