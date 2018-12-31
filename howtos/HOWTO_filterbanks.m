% Some example for running filterbanks and computers such as stft computer
%
% This filterbank library is the implementation of paper:
% "Learning Filter Banks from Raw Speech for Phone Recognition": 
% https://arxiv.org/abs/1711.01161
% 
% Please find the path management in the startup.m script in the root directory
% of this repository. Note that by starting matlab in the root directory, this
% script should automatically run. If it is not the case, you can also move to the
% root directory and run this script manually. 
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


clear all;

%% Using Scaling Functions
% Mel Scaling Function
% Please use the Voicebox 
% (http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html)
% to convert between Hz, Mel, Erb and MIDI frequency scales
disp('Example - Using Mel Scaling Function')
scaling_function = MelScaling();
hertz = 200;
scale = scaling_function.hertz_to_scale(hertz);
disp('done')

%% Using Window Functions
% Hann Window, Hamming Window, Gamma Window
disp('Example - Using Window Function: Impulse Response of a Gamma Window')
window_size = 1000;
peak_ratio = 0.75; 
order = 2;
window = GammaWindow(order, peak_ratio).get_impulse_response(window_size);
disp('done')

%% Using Compute Method to Compute FilterBanks
% Example: Compute a Standard GaborFilterBank
disp('Example - Compute a Standard GaborFilterBank')

bank = GaborFilterBank(MelScaling());
frame_length_ms = 25;
frame_shift_ms = 10;
frame_style = 'causal'; %{'causal', 'centered'}
include_energy = true;  %{true, false} 
pad_to_nearest_power_of_two = true; % {true, false}
use_log = true; %{true, false}
use_power = true; %{true, false}
                       
% Create a random signal, for example
signal_len = 2 ^ 10;  
signal = rand(signal_len,1);    % Example Lengths: 'empty buffer', 
                                % 'length 1 buffer', 
                                % 'medium buffer'(2^8), 
                                % 'large buffer'(2^10)
% disp(signal)

% Create Computer
computer = ShortTimeFourierTransformFrameComputer(...
                            bank, ...
                            frame_length_ms, ...
                            frame_shift_ms, ...
                            frame_style, ...
                            include_energy, ...
                            pad_to_nearest_power_of_two, ...
                            use_log, ...
                            use_power);

% Below are different methods to compute the filterbank
    % Use computer.compute_full to obtain the features
    feats_full = computer.compute_full(signal);
    % disp(feats_full) 

    % Use frame_by_frame_calculation to obtain the features
    feats_framewise = Util.frame_by_frame_calculation(computer, signal);
    % disp(feats_framewise) 

    % Compute a piece of signal
    feats_chunk = computer.compute_chunk(signal(1:ceil(length(signal)/4)));
    feats_chunk = computer.finalize();
    % disp(feats_chunk)

disp('done')

        
        
        
        
        
        
        
        
        
        
        