% Compute features of a signal by integrating STFTs
% 
% Description
% Computations are per frame and as follows:
% 1. The current frame is multiplied with some window (rectangular,
%    Hamming, Hanning, etc)
% 2. An DFT is performed on the result
% 3. For each filter in the provided input bank:
%    a. Multiply the result of 2. with the frequency response of the
%       filter
%    b. Sum either the pointwise square or absolute value of elements
%       in the buffer from 3a.
%    c. Optionally take the log of the sum
% 
% .. warning:: This behaviour differs from that of [1]_ or [2]_ in
%              three ways. First, the sum (3b) comes after the
%              filtering (3a), which changes the result in the squared
%              case. Second, the sum is over the full power spectrum,
%              rather than just between 0 and the Nyquist. This
%              doubles the value at the end of 3c. if a real filter is
%              used. Third, frame boundaries are calculated
%              diffferently.
% 
% Properties
% bank : LinearFilterBank
% frame_length_ms : float, optional
%     The length of a frame, in milliseconds. Defaults to the length
%     of the largest filter in the bank
% frame_shift_ms : float, optional
%     The offset between successive frames, in milliseconds
% frame_style : {'causal', 'centered'}, optional
%     Defaults to ``'centered'`` if `bank.is_zero_phase`, ``'causal'``
%     otherwise.
% include_energy : bool, optional
% pad_to_nearest_power_of_two : bool, optional
%     Whether the DFT should be a padded to a power of two for
%     computational efficiency
% window_function : pydrobert.speech.filters.WindowFunction, dict, or str
%     The window used in step 1. Can be a WindowFunction or something
%     compatible with
%     `pydrobert.speech.alias_factory_subclass_from_arg`. Defaults to
%     `pydrobert.speech.filters.GammaWindow` when ``frame_style`` is
%     ``'causal'``, otherwise `pydrobert.speech.filters.HannWindow`.
% use_log : bool, optional
%     Whether to take the log of the sum from 3b.
% use_power : bool, optional
%     Whether to sum the power spectrum or the magnitude spectrum
% kaldi_shift : bool, optional
%     If ``True``, the ``k``th frame will be computed using the signal
%     between ``signal[
%         k - frame_length // 2 + frame_shift // 2:
%         k + (frame_length + 1) // 2 + frame_shift // 2]``.
%     These are the frame bounds for Kaldi [1]_.
% 
% 
% References
% [1] Povey, D., et al (2011). The Kaldi Speech Recognition
%        Toolkit. ASRU
% [2] Young, S. J., et al (2006). The HTK Book, version 3.4.
%        Cambridge University Engineering Department
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


classdef ShortTimeFourierTransformFrameComputer <  LinearFilterBankFrameComputer
    
    properties
        bank
        frame_length_ms = double.empty
        frame_shift_ms = 10
        frame_style = []
        include_energy
        pad_to_nearest_power_of_two = true
        use_log = true
        use_power = false
        window_function = [] % None
        kaldi_shift = false
        sampling_rate
        frame_length
        frame_shift
        num_coeffs
        started
    end
    
    properties
        log
        power
        real
        first_frame
        buf_len
        chunk_dtype
        buf
        window
        dft_size
        nonlin_op
        truncated_filts
        filt_start_idxs
    end
    
    methods
        
        function obj = ShortTimeFourierTransformFrameComputer(...
            bank, ...
            frame_length_ms, ...
            frame_shift_ms, ...
            frame_style, ...
            include_energy,...
            pad_to_nearest_power_of_two, ...
            use_log, ...
            use_power, ...
            window_function,...
            kaldi_shift)
        
            % Constructor of ShortTimeFourierTransformFrameComputer
            % 
            % Inputs
            % bank : FilterBank objects
            % frame_length_ms : float
            % frame_shift_ms : float
            % frame_style : {'causal', 'centered'}
            % include_energy : bool
            % pad_to_nearest_power_of_two : bool
            % use_log : bool
            % use_power : bool
            % window_function : float
            % frame_shift : window funcion object
            % kaldi_shift : bool
            %
            % Outputs
            % obj : ShortTimeFourierTransformFrameComputer object
            %
            % Exmaple
            % >> bank = GaborFilterBank(MelScaling());
            % >> frame_length_ms = 25;
            % >> frame_shift_ms = 10;
            % >> frame_style = 'causal';
            % >> include_energy = true;
            % >> pad_to_nearest_power_of_two = true;
            % >> use_log = true;
            % >> use_power = true;
            % computer = ShortTimeFourierTransformFrameComputer(...
            %                                         bank, ...
            %                                         frame_length_ms, ...
            %                                         frame_shift_ms, ...
            %                                         frame_style, ...
            %                                         include_energy, ...
            %                                         pad_to_nearest_power_of_two, ...
            %                                         use_log, ...
            %                                         use_power);
            %
                                    
            % Assign Default Values
            if nargin == 9
                kaldi_shift = false;
            elseif nargin == 8
                kaldi_shift = false;
                window_function = [];
            elseif nargin == 7
                kaldi_shift = false;
                window_function = [];
                use_power = false;
            elseif nargin == 6
                kaldi_shift = false;
                window_function = [];
                use_power = false;
                use_log = true;
            elseif nargin == 5
                kaldi_shift = false;
                window_function = [];
                use_power = false;
                use_log = true;
                pad_to_nearest_power_of_two = true;
            elseif nargin == 4
                kaldi_shift = false;
                window_function = [];
                use_power = false;
                use_log = true;
                pad_to_nearest_power_of_two = true;
                include_energy = false;           
            end

            % Call superclass for construction
            obj = obj@LinearFilterBankFrameComputer(bank, include_energy);
            
            % Assign Properties
            obj.bank = bank; %%% How to pass object as an argument
            obj.sampling_rate = bank.sampling_rate;
            obj.frame_shift = round(0.001 * frame_shift_ms * obj.sampling_rate);
            obj.log = use_log;
            obj.power = use_power;
            obj.real = bank.is_real;
            obj.started = false;
            obj.first_frame = true;
            obj.buf_len = 0;
            obj.chunk_dtype = 'double';
            obj.kaldi_shift = kaldi_shift;
            obj.include_energy = include_energy;
            
            if isempty(frame_style)
                if bank.is_zero_phase == true
                    frame_style = 'centered';
                else
                    frame_style = 'causal';
                end
            elseif ~ismember(frame_style, {'centered', 'causal'})
                error('Invalid frame_style type! ');
            end
            obj.frame_style = frame_style;
            
            if isempty(frame_length_ms)
                res = zeros(length(bank.supports));
                for x = bank.supports % size(a) = (2*n) % Matlab has a weried shape so pay attention
                    append(res, x(2)-x(1));
                end
                max_length = max(res);
                
                res = zeros(length(bank.supports_hz));
                for x = bank.supports_hz % size(a) = (2*n) % Matlab has a weried shape so pay attention
                    append(res, x(2)-x(1));
                end
                min_length_hz = min(res);
                calc_length = max(max_length, fix(2*obj.sampling_rate/min_length_hz));
                
                obj.frame_length = max(max_length, calc_length);
            else
                obj.frame_length = fix(0.001 * frame_length_ms * bank.sampling_rate);
            end
            obj.buf = zeros(obj.frame_length, 1); % double.empty
          
            if isempty(window_function)
                if strcmp(frame_style, 'causal')
                    window_function = GammaWindow(); 
                else
                    window_function = HannWindow(); 
                end
            else
                obj.window_function = window_function;
            end
            obj.window = window_function.get_impulse_response(obj.frame_length); 
            
            if pad_to_nearest_power_of_two
                obj.dft_size = round(2^ceil(log2(obj.frame_length)));
            else
                obj.dft_size = obj.frame_length;
            end
              
            if obj.power
                obj.nonlin_op = @(x) norm(x) ^ 2;
            else
                obj.nonlin_op = @(x) sum(abs(x));
            end
            
            obj.truncated_filts = [];
            obj.filt_start_idxs = [];
            
            obj.num_coeffs = bank.num_filts + obj.include_energy;
            for filt_idx = 1:bank.num_filts
                [start_idx, truncated_filt] = bank.get_truncated_response(filt_idx, obj.dft_size); % Here!
                obj.filt_start_idxs = [obj.filt_start_idxs;start_idx];
                obj.truncated_filts = [obj.truncated_filts;truncated_filt];
            end
        end  
    end
    
    
    methods
        function coeffs = compute_full(obj, signal)
            % Produces features using full signal. 
            %
            % Inputs
            % signal : 1D column array of double 
            % 
            % Outputs
            % coeffs : Computed result, 
            %       feature vecotr with size(num_frames, obj.num_coeffs)
            %
            % Example:
            % Create computer object
            % >> computer = ShortTimeFourierTransformFrameComputer(...
            %                                         bank, ...
            %                                         frame_length_ms, ...
            %                                         frame_shift_ms, ...
            %                                         frame_style, ...
            %                                         include_energy, ...
            %                                         pad_to_nearest_power_of_two, ...
            %                                         use_log, ...
            %                                         use_power);
            % >> feats_full = computer.compute_full(buff);
            %
            
            if obj.started
            	error('Already started computing frames');
            end
            frame_length = obj.frame_length;
            frame_shift = obj.frame_shift;
            if length(signal) < floor(frame_length/2) + 1
                coeffs = double.empty(0,obj.num_coeffs);
            else
                if strcmp(obj.frame_style, 'causal')
                    pad_left = 0;
                elseif obj.kaldi_shift == true
                    pad_left = floor(frame_length/2) - floor(frame_shift/2);
                else
                    pad_left = floor((obj.frame_length + 1)/2) - 1;
                end
                num_frames = max(0, floor((length(signal) + floor(frame_shift/2))/frame_shift));
                total_len = (num_frames - 1) * frame_shift - pad_left + frame_length;
                pad_right = max(0, total_len - length(signal));
                if pad_left ~= 0 || pad_right ~= 0
                    signal = padarray(signal,pad_left,'symmetric','pre');
                    signal = padarray(signal,pad_right,'symmetric','post');
                end
                coeffs = zeros(num_frames, obj.num_coeffs);

                for frame_idx = 1:num_frames
                    frame_left = (frame_idx-1) * frame_shift + 1;
                    temp = obj.compute_frame(...
                        signal(frame_left:frame_left + frame_length - 1),...
                        coeffs(frame_idx,:)...
                    );
                    coeffs(frame_idx,:) = temp(:);
                end
            end
        end % function
        
        function coeffs = compute_frame(obj, frame, coeffs)
            % Produces features using a piece of signal, i.e. a frame. 
            %
            % Inputs
            % frame : 1D column array of double 
            % coeffs : Uncalculated feature vecotr
            %
            % Outputs
            % coeffs : Computed result, 
            %       Calculated feature vecotr
            %
            % Example:
            % Please refer to
            % ShortTimeFourierTransformFrameComputer.compute_full for
            % detailed usage.
            %
            % Create computer object
            % >> computer = ShortTimeFourierTransformFrameComputer(...
            %                                         bank, ...
            %                                         frame_length_ms, ...
            %                                         frame_shift_ms, ...
            %                                         frame_style, ...
            %                                         include_energy, ...
            %                                         pad_to_nearest_power_of_two, ...
            %                                         use_log, ...
            %                                         use_power);
            % >> feats_frame = computer.compute_frame(frame, coeffs);
            %
            
            assert(length(frame) == obj.frame_length, "Length of frame must equal to object's frame_length attribute");
            assert(length(coeffs) == obj.num_coeffs, "Length of coeffs must equal to object's num_coeffs attribute");

            energy_coeff = [];
            if obj.include_energy
                coeffs(1) = Util.inner(frame, frame) / obj.frame_length;
                energy_coeff = coeffs(1);
                if obj.power == false
                    coeffs(1) = coeffs(1) ^.5;
                    energy_coeff = coeffs(1); 
                end
                if obj.log == true
                    coeffs(1) = log(max(coeffs(1), Config.LOG_FLOOR_VALUE));
                    energy_coeff = coeffs(1);
                end
                coeffs = coeffs(2:end);
            end
            
            % v_rfft() is from the Voicebox toolkit
            half_spect = transpose(v_rfft(frame .* obj.window, obj.dft_size));  
            half_len = length(half_spect);
            for filt_idx = 1:length(obj.filt_start_idxs)
                start_idx = obj.filt_start_idxs(filt_idx)+1;
                truncated_filt = obj.truncated_filts(filt_idx);
                trunc_len = length(truncated_filt);
                consumed = 0;  % index
                conjugate = false;
                val = 0;
                while consumed < trunc_len
                    if conjugate
                        seg_len = min(...
                            start_idx + trunc_len - consumed,...
                            half_len - 2 + mod(half_len, 2)...
                            )-start_idx;
                        seg_len = max(0, seg_len);
                        if seg_len
                            val = val + ...
                                obj.nonlin_op(conj(half_spect(...
                                (-2 + mod(half_len, 2) - start_idx):...
                                (-2 + mod(half_len, 2) - start_idx - seg_len):...
                                end))) * truncated_filt(consumed:...
                                consumed + seg_len);
                        end
                        start_idx = start_idx - half_len - 2 + ...
                            mod(half_len, 2);
                    else
                        seg_len = min(start_idx + trunc_len - consumed, half_len);
                        seg_len = seg_len - start_idx;
                        seg_len = max(0, seg_len);                        
                        if seg_len ~= 0
                            val = vpa(obj.nonlin_op(...
                                vpa(half_spect(start_idx:...
                                (start_idx + seg_len))) ...
                                * vpa(truncated_filt(...
                                consumed+1:min(consumed+1 + seg_len, end)))...
                                ));
                            start_idx = start_idx - half_len;
                        end
                    end
                    conjugate = ~conjugate;
                    consumed = consumed + seg_len;
                    start_idx = max(0, start_idx);
                end
                
                if obj.real
                    val = val*2;
                end
                if obj.log
                    val = log(max(val, Config.LOG_FLOOR_VALUE));
                end
                coeffs(filt_idx) = val;
            end
            coeffs = [energy_coeff coeffs];
        end % function
        
        function coeffs = compute_chunk(obj, chunk)
            % Produces features using a chunk of signal.
            % 
            % Description
            % Another method to computer feature matrix.
            % Algorithm should work when frame shift is greater than frame
            % length - buf_len may be negative, which will skip samples
            %
            % Inputs
            % chunk : A 1D double array of a chunk of the target signal
            % 
            % Outputs
            % coeffs : computed feature matrix of the chunks of signal,
            %           in a dimension of (num_frames, obj.num_coeffs)
            %
            % Example
            % Please refer to Util.frame_by_frame_calculation() for example
            % of usage
            %
                        
            buf_len = obj.buf_len;
            chunk_len = length(chunk);
            total_len = chunk_len + buf_len;
            noncausal_first = strcmp(obj.frame_style, 'centered');
            noncausal_first = noncausal_first && obj.first_frame;
            
            if noncausal_first
                if obj.kaldi_shift
                    frame_length = floor((obj.frame_length + 1) / 2);
                    frame_length = frame_length + floor(obj.frame_shift / 2);
                else
                    frame_length = floor(obj.frame_length/2) + 1;
                end
            else
                frame_length = obj.frame_length;
            end
            
            frame_shift = obj.frame_shift;
            num_frames = max(0, floor((total_len - frame_length)/frame_shift)+1);
            
            coeffs = zeros(num_frames, obj.num_coeffs); 
    
            for frame_idx = 1:num_frames
                frame_start_idx = (frame_idx-1) * frame_shift+1;
                if frame_start_idx < buf_len
                    frame = [... % Concatenate
                        obj.buf(end-(buf_len - frame_start_idx):end);...
                        chunk(1:frame_length - buf_len + frame_start_idx-1)]; % Added end to accomodate the index
                else
                    frame = chunk(...
                        frame_start_idx - buf_len: ...
                        frame_start_idx - buf_len + frame_length - 1);
                end
                           
                if noncausal_first
                    % the first frame's l.h.s is a reflection of its r.h.s.
                    % shove the reflection into the buf - later frames
                    % might need it
                    chunk = chunk((frame_length - buf_len):end);
                    chunk_len = chunk_len - (frame_length - buf_len);
                    frame_length = obj.frame_length;
                    if obj.kaldi_shift
                        obj.buf(:) = padarray(frame, ...
                            floor(obj.frame_length /2) - floor(obj.frame_shift / 2),...
                            'symmetric', 'pre');
                    else                        
                        obj.buf(:) = padarray(frame, ...
                            floor((frame_length + 1) / 2) - 1,...
                            'symmetric', 'pre');
                    end
                    frame = obj.buf;
                    total_len = chunk_len + frame_length;
                    buf_len = frame_length;
                    noncausal_first = false;          
                end
                coeffs(frame_idx,:) = obj.compute_frame(frame, coeffs(frame_idx,:));
                obj.first_frame = false;
            end
            rem_len = total_len - num_frames * frame_shift;
            assert(rem_len < frame_length);
            if rem_len > 0
                throw_away = total_len - rem_len;

                if throw_away < buf_len
                    rem_ring_len = buf_len - throw_away;
                    assert(rem_ring_len < rem_len);
                    assert(~(rem_ring_len <= rem_len && isempty(chunk)));
                    assert(rem_ring_len <= rem_len);
                    obj.buf(...
                        obj.frame_length - rem_len:...
                        obj.frame_length - rem_len + rem_ring_len - 1 ...
                        ) = obj.buf(obj.frame_length - rem_ring_len + 1 :end);
                    
                    obj.buf(...
                        obj.frame_length - (rem_len - rem_ring_len):end)...
                        = chunk;
                else
                    obj.buf(end-(rem_len-1):end) = chunk(end-(rem_len-1):end);
                end
            end
            obj.buf_len = rem_len;
            obj.started = true;
        end % function
        
        function coeffs = finalize(obj)
            % Finializing feature matrix after calls of compute_chunk
            %
            % Outputs
            % coeffs : computed feature matrix of the chunks of signal,
            %           in a dimension of (num_frames, obj.num_coeffs)
            %
            % Example
            % Please refer to Util.frame_by_frame_calculation() for example
            % of usage
            %
            buf_len = obj.buf_len;
            frame_length = obj.frame_length;
            frame_shift = obj.frame_shift;
            
            if strcmp(obj.frame_style, 'causal')
                pad_left = 0;
            elseif obj.kaldi_shift
                pad_left = floor(frame_length/ 2) - floor(frame_shift/ 2);
            else
                pad_left = floor((frame_length + 1)/ 2) - 1;
            end
            num_frames = buf_len + floor(frame_shift/ 2);
            
            if obj.first_frame == false
                num_frames = num_frames - pad_left;
                pad_left = 0;
            end
            num_frames = floor(num_frames/ frame_shift);
            
            if num_frames >= 1
                pad_right = (num_frames - 1) * frame_shift + ...
                    frame_length - buf_len;
                pad_right = pad_right - pad_left;
                coeffs = zeros(...
                    num_frames, obj.num_coeffs);
                
                frames = padarray(obj.buf(end-buf_len+1:end), ...
                    pad_left, 'symmetric','pre');
                frames = padarray(obj.buf(end-buf_len+1:end), ...
                    pad_right, 'symmetric','post');
                
                for frame_idx = 1: num_frames
                    frame = frames(...
                        (frame_idx-1) * frame_shift + 1:...
                        (frame_idx-1) * frame_shift + frame_length ...
                    );
                    coeffs(frame_idx, :) = obj.compute_frame(frame, coeffs(frame_idx, :));
                end
            else
                coeffs = double.empty(0, obj.num_coeffs); % dtype=self._chunk_dtype)
            end
            obj.buf_len = 0;
            obj.started = false;
            obj.first_frame = true;
        end %function

    end % methods
end % class



 
            
            
            
            
            
            
            
            
             
            
            
            
            
            
            
            
            
             
            
            
            
            
            
            
            
            
            



