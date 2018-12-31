% Some Utility Functions used in the Library
% 
% Description
% Including hertz to angular coversion, 
% frame_by_frame_calculation used for computer conputation,
% and inner() to calculate inner product of two matrix.
% 
% Methods
% angular = hertz_to_angular(hertz, samp_rate)
% hertz = angular_to_hertz(angle, samp_rate)
% coeffs = frame_by_frame_calculation(computer, signal, chunk_size)
% y = inner(a,b)
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

classdef Util < handle

    methods (Static)
        function angular = hertz_to_angular(hertz, samp_rate)
            % Convert cycles/sec to radians/sec
            %
            % Input
            % samp_rate : [Hz] sampling rate in Hz
            %
            % Output
            % hertz : [radian] sampling rate in radian
            %
            % Example
            % >> samp_rate = 16
            % >> Angular = hertz_to_angular(samp_rate);
            %
            angular = hertz * 2 * pi / samp_rate;
        end
        
        function hertz = angular_to_hertz(angle, samp_rate)
            % Convert radians/sec to cycles/sec
            %
            % Input
            % samp_rate : [radian] sampling rate in radian
            %
            % Output
            % hertz : [Hz] sampling rate in Hz
            %
            % Example
            % >> samp_rate = 16
            % >> Hertz = angular_to_hertz(samp_rate);
            %
            hertz = angle * samp_rate / (2 * pi);
        end

        function coeffs = frame_by_frame_calculation(computer, signal, chunk_size) 
            % Compute feature representation of entire signal iteratively
            % 
            % Description
            % This function constructs a feature matrix of a signal through
            % successive calls to `computer.compute_chunk`. Its return value
            % should be identical to that of calling
            % `computer.compute_full(signal)`, but is possibly much slower.
            % `computer.compute_full` should be favoured.
            % 
            % Inputs
            % computer   : Computer Object
            % signal     : array-like,  A 1D float array of the entire signal
            % chunk_size : int
            %              Length of the signal buffer to process at a given time
            % 
            % Outputs
            % coeffs :  A 2D float array of shape ``(num_frames, num_coeffs)``.
            % ``num_frames`` is nonnegative (possibly 0). Contains some number
            % of feature vectors, ordered in time over axis 0.
            % 
            % Example
            % >> bank = GaborFilterBank(MelScaling());
            % >> frame_length_ms = 25;
            % >> frame_shift_ms = 10;
            % >> frame_style = 'causal';
            % >> include_energy = true;
            % >> pad_to_nearest_power_of_two = true;
            % >> use_log = true;
            % >> use_power = true;
            % >> computer = ShortTimeFourierTransformFrameComputer(...
            %                                         bank, ...
            %                                         frame_length_ms, ...
            %                                         frame_shift_ms, ...
            %                                         frame_style, ...
            %                                         include_energy, ...
            %                                         pad_to_nearest_power_of_two, ...
            %                                         use_log, ...
            %                                         use_power);
            % feats_framewise = Util.frame_by_frame_calculation(computer, buff);
            %
            
            if nargin == 2
                chunk_size = 2 ^ 10;
            end
            if (computer.started == true)
                error('Already started computing frames');
            end
            coeffs = [];
            
            while length(signal) > 0
                coeff = computer.compute_chunk(signal(1:min(chunk_size, length(signal))));
                coeffs = [coeffs; ...
                            coeff];
                if chunk_size <= length(signal)
                    signal = signal(chunk_size+1:end);
                else
                    signal = [];
                end
            end
            
            chunk_finalize = computer.finalize();
            
            coeffs = [coeffs; chunk_finalize];
            coeffs = horzcat(coeffs);
        end
       
        function res = inner(a,b)     
            % Compute the inner product of two vectors a and b. 
            %
            % Inputs
            % a : vector
            % b : vector
            %   If a and b are nonscalar, their last dimensions must match.
            %
            % Outputs: 
            % res : The value of the inner product of a and b.
            %       res.shape = a.shape[:-1] + b.shape[:-1]
            %
            % Example
            % >> y = inner(a,b) or inner(a,b)
            %
            
            c = 0;
            n = length(a);		
            for k=1:n		
                c=c+a(k)*b(k);	
            end			
            res = c;
        end
        
    end
end