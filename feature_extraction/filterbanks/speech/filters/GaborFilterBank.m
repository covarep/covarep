% Gabor filters with ERBs between points from a scale
% 
% Description
% Gabor filters are complex, mostly analytic filters that have a
% Gaussian envelope in both the time and frequency domains. They are
% defined as
% 
% .. math::
% 
%      f(t) = C \sigma^{-1/2} \pi^{-1/4}
%             e^{\frac{-t^2}{2\sigma^2} + i\xi t}
% 
% in the time domain and
% 
% .. math::
% 
%      \widehat{f}(\omega) = C \sqrt{2\sigma} \pi^{1/4}
%                            e^{\frac{-\sigma^2(\xi - \omega)^2}{2}}
% 
% in the frequency domain. Though Gaussians never truly reach 0, in
% either domain, they are effectively compactly supported. Gabor
% filters are optimal with respect to their time-bandwidth product.
% 
% `scaling_function` is used to split up the frequencies between
% `high_hz` and `low_hz` into a series of filters. Every subsequent
% filter's width is scaled such that, if the filters are all of the
% same height, the intersection with the precedent filter's response
% matches the filter's Equivalent Rectangular Bandwidth (erb == True)
% or its 3dB bandwidths (erb == False). The ERB is the width of a
% rectangular filter with the same height as the filter's maximum
% frequency response that has the same L_2 norm.
% 
% Properties
% centers_hz : [2-D Array]
% is_real : [bool]
% is_analytic : [bool]
% num_filts : [int]
% sampling_rate : [double]
% supports_hz : [double] 2-D Array
% supports : [double] 2-D Array
% supports_ms : [double] 2-D Array
% scaled_l2_norm : [bool]
% erb : [bool]
% 
% See Also
% Config.EFFECTIVE_SUPPORT_THRESHOLD : the absolute
%     value below which counts as zero
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
%  Sean Robertson <sdrobert@cs.toronto.ca>
%
    
classdef GaborFilterBank < LinearFilterBank
      
    properties
        centers_hz
        is_real
        is_analytic = false
        num_filts = 40
        sampling_rate = 16000
        supports_hz
        supports
        supports_ms
    end
    
    properties (Access=private)
        scaled_l2_norm = false
        erb = false
        wrap_below
        centers_ang
        stds
        supports_ang
        wrap_supports_ang
        scale_l2_norm
    end
    
    methods
        
        function obj = GaborFilterBank(scaling_function,...
                num_filts,...
                low_hz,...
                sampling_rate,...
                scale_l2_norm,...
                erb,...
                high_hz)
            
            % Constructor for GaborFilterBank           
            %            
            % Input
            % scaling_function : [filterbanks/scales/ScalingFunction]
            %     Dictates the layout of filters in the Fourier domain. Can be
            %     a ScalingFunction
            % num_filts : [int]
            %     The number of filters in the bank
            % high_hz, low_hz : [double], optional
            %     The topmost and bottommost edge of the filters, respectively.
            %     The default for high_hz is the Nyquist
            % sampling_rate : [double], optional
            %     The sampling rate (cycles/sec) of the target recordings
            % scale_l2_norm : [bool]
            %     Whether to scale the l2 norm of each filter to 1. Otherwise the
            %     frequency response of each filter will max out at an absolute
            %     value of 1.
            % erb : [bool]
            %
            % Output
            % obj : [GaborFilterBank object]
            %
            % Example
            % >>scaling_function = MelScaling();
            % >>num_filts = 11;
            % >>low_hz = 0;
            % >>sampling_rate = 8000;
            % >>scale_l2_norm = false;
            % >>erb = false;
            % >>bank = GaborFilterBank(scaling_function,num_filts,... 
            %                low_hz, sampling_rate, scale_l2_norm, erb);
            
            if nargin == 1
                num_filts = 40;
                low_hz = double(20);
                sampling_rate = double(16000);
                scale_l2_norm = false;
                erb = false;
                high_hz = double.empty; 
            elseif nargin == 2
                low_hz = double(20);
                sampling_rate = double(16000);
                scale_l2_norm = false;
                erb = false;
                high_hz = double.empty; 
            elseif nargin == 3
                sampling_rate = double(16000);
                scale_l2_norm = false;
                erb = false;
                high_hz = double.empty; 
            elseif nargin == 4
                scale_l2_norm = false;
                erb = false;
                high_hz = double.empty; 
            elseif nargin == 5
                erb = false;
                high_hz = double.empty; 
            elseif nargin == 6
                high_hz = double.empty; 
            end

            obj = obj@LinearFilterBank();
            obj.scaled_l2_norm = scale_l2_norm;
            obj.erb = erb;
            obj.sampling_rate = sampling_rate;
            if isempty(high_hz)
               high_hz = floor(sampling_rate/2);
            end
            
            % Input validation
            if low_hz < 0 || high_hz <= low_hz || high_hz > floor(sampling_rate/2)
                error('Invalid frequency range: (%d, %d)', low_hz, high_hz);
            end
            scale_low = scaling_function.hertz_to_scale(low_hz);
            scale_high = scaling_function.hertz_to_scale(high_hz);
            scale_delta = (scale_high - scale_low) / (num_filts + 1);
            % edges dictate the points where filters should intersect. We
            % make a pretend intersection halfway between low_hz and
            % the first filter center in the scaled domain. Likewise with
            % high_hz and the last filter center. Intersections are spaced
            % uniformly in the scaled domain
            edges = [];
            for idx= 1:num_filts+1
                edges = [edges ; scaling_function.scale_to_hertz(scale_low + scale_delta * ((idx-1) + .5))];
            end
            centers_hz = [];
            centers_ang = [];
            stds = [];
            supports_ang = [];
            supports = [];
            wrap_supports_ang = [];
            obj.wrap_below = false;
            log_2 = log(2);
            log_pi = log(pi);
            t_support_const = -2 * log(Config.EFFECTIVE_SUPPORT_THRESHOLD);
            f_support_const = t_support_const;
            
            if scale_l2_norm
               f_support_const = f_support_const + log_2 + .5 * log_pi;
               t_support_const = t_support_const - .5 * log_pi;
            else
               t_support_const = t_support_const - log_2 - log_pi;
            end
            if erb
               bandwidth_const = sqrt(pi) / 2;
            else
               bandwidth_const = sqrt(3 / 10 * log(10));
            end
          
            zip_edges = [edges(1:end-1) edges(2:end)];
            
            for ind = 1:size(zip_edges,1)  % 1st dimension of zip_edges 
                left_intersect = zip_edges(ind);
                right_intersect = zip_edges(ind + size(zip_edges,1));
                center_hz = (left_intersect + right_intersect)/2;
                center_ang = Util.hertz_to_angular(center_hz, obj.sampling_rate);
                std = bandwidth_const / Util.hertz_to_angular(center_hz - left_intersect, obj.sampling_rate);
                log_std = log(std);
                
                if scale_l2_norm
                    diff_ang = sqrt(log_std + f_support_const) / std;
                    wrap_diff_ang = sqrt(log_std + f_support_const + log_2) / std;
                    diff_samps = ceil(std * sqrt(t_support_const - log_std));
                else
                    diff_ang = sqrt(f_support_const) / std;
                    wrap_diff_ang = sqrt(f_support_const + log_2) / std;
                    diff_samps = ceil(std * sqrt(t_support_const - 2 * log_std));
                end
                supp_ang_low = center_ang - diff_ang;
                if supp_ang_low < 0
                    obj.wrap_below = true;
                end
                centers_hz = [centers_hz; center_hz];
                centers_ang = [centers_ang; center_ang];
                supports_ang = [supports_ang; [center_ang - diff_ang center_ang + diff_ang]];
                wrap_supports_ang = [wrap_supports_ang; 2 * wrap_diff_ang];
                supports = [supports; [-diff_samps diff_samps]];
                stds = [stds; std];
            end
            
            obj.centers_ang = centers_ang;
            obj.centers_hz = centers_hz;
            obj.stds = stds;
            obj.supports_ang = supports_ang;
            obj.wrap_supports_ang = wrap_supports_ang;
            obj.num_filts = num_filts;
            
            supports_hz =[];
            for idx = 1:size(supports_ang,1) % 1st dimension of supports_ang 
                supports_hz = [supports_hz; ...
                    [Util.angular_to_hertz(supports_ang(idx), obj.sampling_rate) ...
                    Util.angular_to_hertz(supports_ang(idx+size(supports_ang,1)), obj.sampling_rate)]];
            end
            obj.supports_hz = supports_hz;
            obj.supports = supports;
        end
 
    end
   
    methods
        function res = get_impulse_response (obj, filt_idx, width)
            % Helper function to get impulse response of a GaborFilterBank
            %
            % Inputs
            % filt_idx  : [int] filter index
            % width     : [int] width of the signal     
            % 
            % Outputs
            % res : [double]
            % 
            % Example
            % >>scaling_function = MelScaling();
            % >>num_filts = 11;
            % >>low_hz = 0;
            % >>sampling_rate = 8000;
            % >>scale_l2_norm = false;
            % >>erb = false;
            % >>bank = GaborFilterBank(scaling_function,num_filts,... 
            %                low_hz, sampling_rate, scale_l2_norm, erb);
            % >>x = bank.get_impulse_response(filt_idx, dft_size);
            %
            center_ang = obj.centers_ang(filt_idx);
            std = obj.stds(filt_idx);
            res = zeros(width,1); % we want it to store complex values
            if obj.scale_l2_norm
                const_term = -0.5 * log(std) - 0.25 * log(pi);
            else
                const_term = -0.5 * log(2 * pi) - log(std);
            end
            denom_term = 2 * std ^ 2;
            
            for t = 1: width+1
                val = complex(-(t-1) ^ 2 / denom_term + const_term, center_ang * (t-1));
                val = exp(val);
                if t ~= width+1
                    res(t) = res(t) + val;
                end
                if width-t+2 > 0 && t>1
                    res(end-t+2) = res(end-t+2) + conj(val);
                end
            end
        end
        

        function res = get_frequency_response(obj, filt_idx, width, half)
            % Helper function to get frequency response of a GaborFilterBank
            %
            % Inputs
            % filt_idx  : [int] filter index
            % width     : [int] width of the signal     
            % half      : [bool] (def. false)
            %
            % Outputs
            % res : [hz] double
            % 
            % Example
            % >>scaling_function = MelScaling();
            % >>num_filts = 11;
            % >>low_hz = 0;
            % >>sampling_rate = 8000;
            % >>scale_l2_norm = false;
            % >>erb = false;
            % >>bank = GaborFilterBank(scaling_function,num_filts,... 
            %                low_hz, sampling_rate, scale_l2_norm, erb);
            % >>X = bank.get_frequency_response(filt_idx, dft_size, false);
            %
            if nargin == 3
                half = false;
            end
            center_ang = obj.centers_ang(filt_idx);
            lowest_ang = obj.supports_ang(filt_idx,1);
            highest_ang = obj.supports_ang(filt_idx,2);
            lowest_ang = filt_idx;
            std = obj.stds(filt_idx);
            dft_size = width;
            
            if half == true
                if mod(width,2) ~= 0
                    dft_size = floor((width + 1) / 2);
                else
                    dft_size = floor(width / 2) + 1;
                end
            end
            res = zeros(dft_size, 1);
            if obj.scale_l2_norm
                const_term = .5 * log(2 * std) + .25 * log(pi);
            else
                const_term = 0;
            end
            num_term = -(std ^ 2) / 2;
            
            for idx = 1: dft_size
                for period =  -1 - floor(max(-lowest_ang, 0) / (2 * pi)): 2 + floor(highest_ang / (2 * pi)) - 1
                    omega = ((double(idx)-1) / double(width) + double(period)) * 2 * pi;
                    val = num_term * (center_ang - omega) ^ 2 + const_term;
                    val = exp(val);
                    res(idx) = res(idx) + val;
                end
            end
        end
        
        
        function  [start_idx, truncated_response] = get_truncated_response(obj, filt_idx, width)
            % Helper function to get truncated response of a GaborFilterBank
            % 
            % Desription
            % wrap_supports_ang contains the angular supports of each filter
            % if the effective support threshold were halved. If this
            % support exceeds the 2pi period, overlap from aliasing in the
            % periphery will exceed the effective support, meaning the
            % entire period lies in the support
            %
            % Inputs
            % filt_idx  : [int] filter index
            % width     : [int] width of the signal     
            % 
            % Outputs
            % start_idx : [int]
            % truncated_response : [1-D column vector]
            % 
            % Example
            % >>scaling_function = MelScaling();
            % >>num_filts = 11;
            % >>low_hz = 0;
            % >>sampling_rate = 8000;
            % >>scale_l2_norm = false;
            % >>erb = false;
            % >>bank = GaborFilterBank(scaling_function,num_filts,... 
            %                low_hz, sampling_rate, scale_l2_norm, erb);
            % [start_idx, truncated_response] = bank.get_truncated_response(filt_idx, dft_size);
            %
            % Please refer to test_fbank.test_truncated_matches_full
            % for more information
            %
            
            if obj.wrap_supports_ang(filt_idx) >= 2 * pi
                start_idx = 1;
                truncated_response = obj.get_frequency_response(filt_idx, width);            
            else
                center_ang = obj.centers_ang(filt_idx);
                std = obj.stds(filt_idx);
                lowest_ang = obj.supports_ang(filt_idx, 1);
                highest_ang = obj.supports_ang(filt_idx, 2);
                left_idx = ceil(width * lowest_ang / (2 * pi));
                right_idx = round(width * highest_ang / (2 * pi));
                res = zeros(1 + right_idx - left_idx, 1);
                
                if obj.scale_l2_norm
                    const_term = .5 * log(2 * std) + .25 * log(pi);
                else
                    const_term = 0;
                end
                num_term = -(std ^2) / 2;
                
                for idx = left_idx: right_idx
                    for period = - round(max(-lowest_ang, 0) / (2 * pi)): 1 + round(highest_ang / (2 * pi)) - 1 % -1 because the syntax is different in matlab and python
                        omega = ((idx-1) / width + period) * 2 * pi;
                        val = num_term * (center_ang - omega) ^ 2 + const_term;
                        val = exp(val);
                        res(idx - left_idx + 1) = res(idx - left_idx + 1) + val; % Need to check if this is right
                    end
                end
                start_idx = mod(left_idx, width);
                truncated_response = res;
            end
        end
    end
    
end



