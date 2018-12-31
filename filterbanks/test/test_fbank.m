% Test file for FBank
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

classdef test_fbank < matlab.unittest.TestCase
    properties
        bank
    end 
    
    properties (TestParameter)
        scaling_function = {MelScaling()}
        num_filts = {1, 11}
        high_hz = {8000}
        low_hz = {0}
        sampling_rate = {16000}
        analytic = {true}
    end

    methods (Test)
        function test_validate(testCase, ...
            num_filts, high_hz, low_hz, sampling_rate, analytic)
        
            testCase.bank = FBank(num_filts, high_hz, low_hz, ...
                sampling_rate, analytic);
            testCase.verifyEqual(testCase.bank.isvalid, true);
            
        end
        
        function test_half_response_matches_full(testCase, ...
            num_filts, high_hz, low_hz, sampling_rate, analytic)
        
            testCase.bank = FBank(num_filts, high_hz, low_hz, ...
                sampling_rate, analytic);
            
            for filt_idx = 1:testCase.bank.num_filts
                dft_size = testCase.bank.supports(filt_idx,2) - testCase.bank.supports(filt_idx,1);
                Xh = testCase.bank.get_frequency_response(filt_idx, dft_size, true); 
                X = testCase.bank.get_frequency_response(filt_idx, dft_size, false);
                testCase.verifyEqual(X(1:length(Xh)), Xh, 'AbsTol', 1e-3);
            end 
        end
        
        function test_frequency_matches_impulse(testCase, ...
            num_filts, high_hz, low_hz, sampling_rate, analytic)
        
            testCase.bank = FBank(num_filts, high_hz, low_hz, ...
                sampling_rate, analytic);
            
            for filt_idx = 1:testCase.bank.num_filts
                left_hz = testCase.bank.supports_hz(filt_idx, 1);
                right_hz = testCase.bank.supports_hz(filt_idx, 2);
                left_samp = testCase.bank.supports(filt_idx,1);
                right_samp = testCase.bank.supports(filt_idx,2);
                required_freq_size = 2 * testCase.bank.sampling_rate / (right_hz - left_hz);
                required_temp_size = right_samp - left_samp;
                if required_temp_size < 5 || required_freq_size < 5
                    % TODO: this is a stopgap for when filters are
                    % *too* localized in time or frequency. This'll cause too
                    % much attenuation/gain in one domain or the other.
                    continue
                end
                % allow over- or under-sampling
                dft_size = floor(max((right_samp - left_samp), 2 * testCase.bank.sampling_rate / (right_hz - left_hz)));
                X = testCase.bank.get_frequency_response(filt_idx, dft_size);
                x = testCase.bank.get_impulse_response(filt_idx, dft_size);
                if analytic
                    testCase.assertEqual(ifft(X), x, 'AbsTol', 5e-3);
                else
                    testCase.assertEqual(real(ifft(X)), x, 'AbsTol', 5e-3);
                end
            end 
        end
        
        function test_truncated_matches_full(testCase, ...
            num_filts, high_hz, low_hz, sampling_rate, analytic)
            
            testCase.bank = FBank(num_filts, high_hz, low_hz, ...
                sampling_rate, analytic);
                        
            for filt_idx = 1:testCase.bank.num_filts
                left_hz = testCase.bank.supports_hz(filt_idx, 1);
                right_hz = testCase.bank.supports_hz(filt_idx, 2);
                left_samp = testCase.bank.supports(filt_idx,1);
                right_samp = testCase.bank.supports(filt_idx,2);

                % initialize the random number generator to make the 
                % results repeatable
                rng();
                
                dft_size = round(...
                        max([...
                        (right_samp - left_samp) * (1 + rand), ...
                        2 * testCase.bank.sampling_rate / (right_hz - left_hz), ...
                        1 ...
                        ]));
                    
                left_period = round(floor(left_hz / testCase.bank.sampling_rate));
                right_period = round(ceil(right_hz / testCase.bank.sampling_rate));
                full_response = testCase.bank.get_frequency_response(filt_idx, dft_size);
                [bin_idx, truncated] = testCase.bank.get_truncated_response(filt_idx, dft_size);
                challenge = zeros(dft_size, 1); % dtype=truncated.dtype
                wrap = min((bin_idx-1) + length(truncated), dft_size) - (bin_idx-1);
                challenge(bin_idx:bin_idx + wrap-1) = truncated(1:wrap);
                challenge(1:length(truncated) - wrap) = truncated(wrap+1:end);
                
                if testCase.bank.is_real
                    if bin_idx
                        challenge(...
                            length(challenge) - bin_idx - length(truncated) + 2:...
                            length(challenge) - bin_idx +1 ...
                        ) = flip(truncated(:));
                    else
                        temp = flip(truncated(:));
                        challenge(...
                            length(challenge) - bin_idx - length(truncated) + 2:...
                            length(challenge) - bin_idx +1 ...
                        ) = temp(1:end-1);
                    end
                end
                testCase.assertEqual(full_response, challenge, 'absTol', (right_period - left_period) * Config.EFFECTIVE_SUPPORT_THRESHOLD);
            end

        end % Function
        
    end
    
end