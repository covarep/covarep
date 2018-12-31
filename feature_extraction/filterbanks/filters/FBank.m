% A mel-triangular filter bank that is square-rooted
% 
% Description
% An ``Fbank`` instance is intended to replicate the filters from Kaldi and
% HTK. Its scale is fixed to Mel-scale. Like a
% ``TriangularOverlappingFilterBank``, ``Fbank`` places the vertices of
% triangular filters uniformly along the target scale. However, an ``Fbank``
% is triangular in the Mel-scale, whereas the triangular bank is triangular
% in frequency.
% 
% .. note:: In a standard mel-filterbank spectrogram, the power spectrum is
% calculated before filtering. This module's spectrogram takes the
% power spectrum after filtering. To recreate the frequency response
% of the alternate order, we can take the pointwise square root of the
% frequency response.
% 
% Properties
% centers_hz : [Hz]
%     The point of maximum gain in each filter's frequency response
%     This property gives the so-called "center frequencies" - the
%     point of maximum gain - of each filter.
% is_real : [bool]
% is_analytic : [bool]
% num_filts : [int]
% sampling_rate : [2-D Array]
% supports_hz : [tuple]
% supports : [2-D Array]
% supports_ms : [2-D Array]
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

classdef FBank < LinearFilterBank
   
    properties
        centers_hz
        is_real
        is_analytic 
        num_filts
        sampling_rate 
        supports_hz
        supports
        supports_ms
    end
    
    properties (Access=private)
        vertices
        is_zero_phase
    end
    
    methods
        
        function obj = FBank(num_filts, high_hz, low_hz, sampling_rate, ...
                analytic)
            % Constructor Function for FBank
            %
            % Inputs
            % num_filts : int, optional
            %     The number of filters in the bank
            % high_hz, low_hz : float, optional
            %     The topmost and bottommost edge of the filters, respectively.
            %     The default for high_hz is the Nyquist
            %     (def. high_hz = floor(sampling_rate/2)
            %     (def. low_hz = 20)
            % sampling_rate : float, optional
            %     The sampling rate (cycles/sec) of the target recordings
            %     (def. 16000)
            % analytic : bool, optional
            %     Whether to use an analytic form of the bank. The analytic form
            %     is easily derived from the real form in [1]_ and [2]_. Since
            %     the filter is compactly supported in frequency, the analytic
            %     form is simply the suppression of the ``[-pi, 0)`` frequencies
            %     (def. false)
            % 
            % Outputs
            % obj : An FBank Object
            % 
            % Example
            % >>num_filts = 30;
            % >>fbank = FBank(num_filts);
            %     
            
            scaling_function = MelScaling();
            
            % Provide values for defult argument if not provided
            if nargin == 1
                sampling_rate = 16000;
                high_hz = floor(sampling_rate/2);
                low_hz = 20;
                analytic = false;
            elseif nargin == 2
                low_hz = 20;
                sampling_rate = 16000;
                analytic = false;
            elseif nargin == 3
                sampling_rate = 16000;
                analytic = false;
            elseif nargin == 4
                analytic = false;
            end
     
            if (low_hz < 0) || (high_hz <= low_hz) || (high_hz > floor(sampling_rate/2))
                error('Invalid frequency range: (%d, %d)', low_hz, high_hz);
            end
            
            % Call superclass constructor before accessing object
            obj = obj@LinearFilterBank();

            % Compute vertices
            scale_low = scaling_function.hertz_to_scale(low_hz);
            scale_high = scaling_function.hertz_to_scale(high_hz);
            scale_delta = (scale_high - scale_low) / (num_filts + 1);
            vertices = [];
            for idx = 1:num_filts + 2
                vertices = [vertices; scaling_function.scale_to_hertz(scale_low + scale_delta * (idx-1))];
            end
            obj.vertices = vertices;   

            obj.is_analytic = analytic;
            obj.is_real = ~obj.is_analytic;
            obj.is_zero_phase = true;
            obj.num_filts = length(obj.vertices)-2;
            obj.centers_hz = obj.vertices(2:length(obj.vertices)-1);
            obj.sampling_rate = sampling_rate;
            if length(vertices) > 2
                obj.centers_hz = vertices(2:end-1);
            else
                obj.centers_hz = [];
            end
            
            obj.supports_hz = [vertices(1:end-2) vertices(3:end)];
            
            supports = [];

            for idx = 1:length(obj.vertices) - 2
                left = Util.hertz_to_angular(obj.vertices(idx), obj.sampling_rate);
                mid = Util.hertz_to_angular(obj.vertices(idx + 1), obj.sampling_rate);
                right = Util.hertz_to_angular(obj.vertices(idx + 2), obj.sampling_rate);
                K = right - left + 2 * ((right - mid) * (mid - left))^2;
                K = K /(Config.EFFECTIVE_SUPPORT_THRESHOLD)^2 * pi;
                K = K /(right - mid) * (mid - left);
                K = K / sqrt(Config.EFFECTIVE_SUPPORT_THRESHOLD);
                K = K / sqrt(mid - left) * sqrt(right - mid);
                K = K^.3333;
                K = ceil(K);
                supports = [supports; -floor(K/2)-1 floor(K/2) + 1];
            end
            obj.supports = supports;
        end % Function
    end % Methods
    
    methods
        function frequency_response = get_frequency_response(obj, filt_idx, width, half)
            % Helper function to get frequency response of a FBank
            %
            % Inputs
            % filt_idx  : [int] filter index
            % width     : [int] width of the signal     
            % half      : [bool] 
            % 
            % Outputs
            % frequency_response : [Hz] double
            % 
            % Example
            % Initialize FBank object
            % >> scaling_function = MelScaling();
            % >> num_filts = 11;
            % >> high_hz = 8000;
            % >> low_hz = 0;
            % >>  sampling_rate = 16000;
            % >> analytic = true;
            % >> bank = FBank(num_filts, high_hz, low_hz, sampling_rate, analytic);
            % >> Xh = bank.get_frequency_response(filt_idx, dft_size, true); 
            %
            
            if nargin < 3
                error('Not enough input argument, requires: filt_idx, width')
            elseif nargin == 3
                half = false;
            end
            
            scaling_function = MelScaling();
            left_hz = obj.vertices(filt_idx);
            mid_hz = obj.vertices(filt_idx + 1);
            right_hz = obj.vertices(filt_idx + 2);
            left_mel = scaling_function.hertz_to_scale(left_hz);
            mid_mel = scaling_function.hertz_to_scale(mid_hz);
            right_mel = scaling_function.hertz_to_scale(right_hz);
            left_idx = max(ceil(width * left_hz / obj.sampling_rate), 1);
            right_idx = floor(width * right_hz / obj.sampling_rate);
            assert(obj.sampling_rate * (left_idx - 1) / width <= left_hz);
            assert(obj.sampling_rate * (right_idx + 1) / width >= right_hz);
            assert(width>0);
            
            dft_size = width;
            
            if half == true
                if mod(width,2)~=0
                    dft_size = floor((width + 1) / 2);
                else
                    dft_size = floor(width / 2) + 1;
                end
            end
            
            frequency_response = zeros(dft_size, 1); % size of dft_size*1 matrix
            
            for idx = left_idx : min(dft_size, right_idx+1)
                hz = obj.sampling_rate * (idx-1) / width;
                mel = scaling_function.hertz_to_scale(hz);
                if mel <= mid_mel
                    val = (mel - left_mel) / (mid_mel - left_mel);
                else
                    val = (right_mel - mel) / (right_mel - mid_mel);
                end
                frequency_response(idx) = val ^ .5;

                if half==false && obj.is_analytic==false
                    frequency_response(end-idx+1) = val ^.5;
                end
            end
        end
        
        function [left_idx, truncated_response] = get_truncated_response(obj, filt_idx, width)
            % Helper function to get truncated response of a FBank
            % 
            % Inputs
            % filt_idx  : [int] filter index
            % width     : [int] width of the signal     
            % 
            % Outputs
            % left_idx : [int]
            % truncated_response : [1-D column vector]
            % 
            % Example
            % Initialize FBank object
            % >> scaling_function = MelScaling();
            % >> num_filts = 11;
            % >> high_hz = 8000;
            % >> low_hz = 0;
            % >>  sampling_rate = 16000;
            % >> analytic = true;
            % >> bank = FBank(num_filts, high_hz, low_hz, sampling_rate, analytic);
            % [bin_idx, truncated] = bank.get_truncated_response(filt_idx, dft_size);
            %
            % Please refer to test_fbank.test_truncated_matches_full
            % for more information
            %
            
            scaling_function = MelScaling();
            left_hz = obj.vertices(filt_idx);        
            mid_hz = obj.vertices(filt_idx + 1);
            right_hz = obj.vertices(filt_idx + 2);
            left_mel = scaling_function.hertz_to_scale(left_hz);
            mid_mel = scaling_function.hertz_to_scale(mid_hz);
            right_mel = scaling_function.hertz_to_scale(right_hz);
            left_idx = max(ceil(width * left_hz / obj.sampling_rate), 1);
            right_idx = ceil(width * right_hz / obj.sampling_rate);
           
            assert(obj.sampling_rate * (left_idx - 1) / width <= left_hz);
            assert(obj.sampling_rate * (right_idx + 1) / width >= right_hz);
            assert(width > 0);
            
            truncated_response = zeros(min(width, right_idx + 1) - left_idx, 1);
            
            for idx = left_idx:min(width-1, right_idx)
                hz = obj.sampling_rate * (idx-1) / width;
                mel = scaling_function.hertz_to_scale(hz);
                if mel <= mid_mel
                    truncated_response((idx - left_idx)+1) = (mel - left_mel) / (mid_mel - left_mel);
                else
                    truncated_response((idx - left_idx)+1) = (right_mel - mel) / (right_mel - mid_mel);
                end
            end
            truncated_response = truncated_response.^0.5;
        end
        
        
        function impulse_response = get_impulse_response(obj, filt_idx, width)
            % Helper function to get impulse response of a FBank
            %
            % Inputs
            % filt_idx  : [int] filter index
            % width     : [int] width of the signal     
            % 
            % Outputs
            % impulse_response : [double]
            % 
            % Example
            % Initialize FBank object
            % >> scaling_function = MelScaling();
            % >> num_filts = 11;
            % >> high_hz = 8000;
            % >> low_hz = 0;
            % >>  sampling_rate = 16000;
            % >> analytic = true;
            % >> bank = FBank(num_filts, high_hz, low_hz, sampling_rate, analytic);
            % >> X = bank.get_impulse_response(filt_idx, dft_size, true); 
            % 
            
            % For the ease of computation, we simply invert the frequency response
            if obj.is_analytic
                freq_response = obj.get_frequency_response(...
                    filt_idx, width, false);
                impulse_response = ifft(freq_response);
            else
                freq_response = obj.get_frequency_response(...
                    filt_idx, width, true);
                % v_irfft is from the Voicebox toolkit
                impulse_response = v_irfft(freq_response, width); 
            end
        end
    end
    
end



