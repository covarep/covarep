% A collection of linear, time invariant filters
% 
% Description
% A ``LinearFilterBank`` instance is expected to provide factory
% methods for instantiating a fixed number of LTI filters in either
% the time or frequency domain. Filters should be organized lowest
% frequency first.
% 
% Properties
% is_real         : [bool]    Whether the filters are real or complex
% is_analytic     : [bool]    Whether the filters are (approximately) analytic
% is_zero_phase   : [bool]    Whether the filters are zero phase or not
%         Zero phase filters are even functions with no imaginary part
%         in the fourier domain. Their impulse responses center around 0.
% num_filts       : [int]     Number of filters in the bank
% sampling_rate   : [double]  Number of samples in a second of target recording
% centers_hz      : [2-D Array]
% supports_hz     : [Hz] Boundaries of effective support of filter freq responses
%         A tuple of length `num_filts` containing pairs of floats
%         of the low and high frequencies. Frequencies outside the span
%         have a response of approximately (with magnitude up to
%         `speech.EFFECTIVE_SUPPORT_SIGNAL`) zero.
% 
%         The boundaries need not be tight, i.e. the region inside the
%         boundaries could be zero. It is more important to guarantee that
%         the region outside the boundaries is approximately zero.
% 
%         The boundaries ignore the Hermitian symmetry of the filter if it
%         is real. Bounds of ``(10, 20)`` for a real filter imply that the
%         region ``(-20, -10)`` could also be nonzero.
% 
%         The user is responsible for adjusting the for the periodicity
%         induced by sampling. For example, if the boundaries are
%         ``(-5, 10)`` and the filter is sampled at 15Hz, then all bins
%         of an associated DFT could be nonzero.
% 
% supports : [samples] Boundaries of effective support of filter 
%            impulse resps, in samples
%         A tuple of length `num_filts` containing pairs of
%         integers of the first and last (effectively) nonzero samples.
% 
%         The boundaries need not be tight, i.e. the region inside the
%         boundaries could be zero. It is more important to guarantee that
%         the region outside the boundaries is approximately zero.
% 
%         If a filter is instantiated using a buffer that is unable to
%         fully contain the supported region, samples will wrap around the
%         boundaries of the buffer.
% 
%         Noncausal filters will have start indices less than 0. These
%         samples will wrap to the end of the filter buffer when the
%         filter is instantiated.
%
% supports_ms : [ms] Boundaries of effective support of filter 
%               impulse resps, in ms
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

classdef (Abstract) LinearFilterBank < handle

    properties (Abstract)
%         centers_hz
%         is_real
%         is_analytic
%         num_filts
%         sampling_rate
%         supports_hz
%         supports
%         supports_ms
    end
    
    methods (Abstract)
%         get_impulse_response(obj, filt_idx, width)
%         get_frequency_response(obj, filt_idx, width, half) %half = False
%         get_truncated_response(obj, filt_idx, width) 
    end
        
end