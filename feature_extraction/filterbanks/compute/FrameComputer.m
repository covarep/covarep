% Construct features from a signal from fixed-length segments
% 
% Description
% A signal is treated as a (possibly overlapping) time series of
% frames. Each frame is transformed into a fixed-length vector of
% coefficients.
% 
% Properties
% frame_style       : {'causal', 'centered'}
%         Dictates how the signal is split into frames
% 
%         If ``'causal'``, the ``k``th frame is computed over the indices
%         ``signal[k * frame_shift:k * frame_shift + frame_length]`` (at
%         most). If ``'centered'``, the ``k``th frame is computed over
%         the indices
%         ``signal[
%             k * frame_shift - (frame_length + 1) // 2 + 1:
%             k * frame_shift + frame_length // 2 + 1]``. Any range
%         beyond the bounds of the signal is generated in an
%         implementation-specific way.
%
% sampling_rate     : float
%         Number of samples in a second of a target recording
% 
% frame_length      : int
%         Number of samples which dictate a feature vector
% 
% frame_length_ms   : float
%         Number of milliseconds of audio which dictate a feature vector
% 
% frame_shift       : int
%         Number of samples absorbed between successive frame computations
%
% frame_shift_ms    : float
%         Number of milliseconds between succecssive frame computations
%
% num_coeffs        : int
%         Number of coefficients returned per frame
%
% started           : bool
%         Boolean variable indicating if the computing process has been started
% 
% Example
% Features can be computed one at a time, for example:
% >>> chunk_size = 2 ** 10
% >>> while len(signal):Z
% >>>     segment = signal[:chunk_size]
% >>>     feats = computer.compute_chunk(segment)
% >>>     # do something with feats
% >>>     signal = signal[chunk_size:]
% >>> feats = computer.finalize()
% 
% Or all at once (which can be much faster, depending on how the
% computer is optimized):
% >>> feats = computer.compute_full(signal)
% 
% The ``k``th frame can be roughly localized to the signal offset
% to about ``signal[k * computer.frame_shift]``. The signal's exact
% region of influence is dictated by the `frame_style` property.
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

classdef FrameComputer < handle 
    properties (Abstract)
%         frame_style
%         sampling_rate
%         frame_length
%         frame_length_ms
%         frame_shift
%         frame_shift_ms
%         num_coeffs
%         started
    end
    
    methods (Abstract)
        compute_chunk(obj,chunk)
        finalize(obj)        
    end
        
    methods
        function coeffs = compute_full(obj, signal)
            coeffs = Util.frame_by_frame_calculation(obj, signal);
        end
    end
    
end