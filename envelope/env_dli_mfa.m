% Multi-Frame Analysis based on derivative of Linear Interpolation (dLI-MFA)
%  
% Input
%  af             : Array of structures with matrices containing sinusoidal
%                   parameters (as the output of sin_analysis.m).
%                   Each matrix is made of:
%                      The 1st line for the frequency of each sinusoid [Hz]
%                      The 2nd line for their amplitude (linear scale)
%                   The DC has to be included (in the first column)
%  fs             : [Hz] Signal's sampling frequency
%  [extrap_dcny] : If true (default), extrapolate sinusoidal components at DC
%                  and up to Nyquist.
%                  If false, use the sinusoidal components as they are.
%
% Output
%  E              : The amplitude cepstral envelope
%  
% Reference
%  [1] G. Degottex, "A Time Regularization Technique for Discrete Spectral 
%      Envelopes Through Frequency Derivative", Signal Processing Letters, IEEE, 
%      22(7):978-982, July 2015.
%
% Copyright (c) 2013 Foundation for Research and Technology-Hellas - Institute
%                    of Computer Science (FORTH-ICS)
%
% License
%  This file is part of libphoni. libphoni is free software: you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. libphoni is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function [E] = env_dli_mfa(af, fs, dftlen, extrap_dcny)

    if nargin<4; extrap_dcny=true; end

    % If asked, extrapolate components at DC and up to Nyquist
    if extrap_dcny; af = env_extrap_sins_dcny(af, fs); end

    % For each frame, compute the frequency derivative of the linear envelope
    F = fs*(0:dftlen/2)/dftlen;
    dA = ones(numel(af),dftlen/2);
    for fi=1:numel(af)
        A = interp1(af(fi).sins(1,:), log(af(fi).sins(2,:)), F, 'linear', 'extrap');
        dA(fi,:) = diff(A); % derivative approximation
    end

    % Smooth the frequency derivative across time
    win = hamming(size(dA,1));
    win = win./sum(win);
    dAw = dA.*repmat(win, 1, size(dA,2));
    mA = sum(dAw,1);

    % Retrieve the envelope
    E = cumsum(mA);
    E = [E(1), E];

    % Align the envelope on the harmonics of the central frame
    ci = floor((numel(af)-1)/2)+1;
    ek = interp1(F', E, af(ci).sins(1,:));
    ak = log(af(ci).sins(2,:));
    idx = find(af(ci).sins(1,:)>0 & af(ci).sins(1,:)<4000);
    E = E + mean(ak(idx)) - mean(ek(idx));

    E = exp(E);

    % Plot the final solution
    if 0
        subplot(211);
            hold off;
            plot(F(1:end-1), dA, 'r');
            hold on;
            plot(F(1:end-1), mA, 'b');
            xlim([0 fs/2]);
            xlabel('Frequency [Hz]');
        subplot(212);
            hold off;
            plot(F, mag2db(E), 'b');
            hold on;
            plot(af(ci).sins(1,:), mag2db(af(ci).sins(2,:)), 'xk');
            xlim([0 fs/2]);
            xlabel('Frequency [Hz]');
            ylabel('Amplitude [dB]');
            title(['dLI-MFA']);
        pause
%          keyboard
    end

return
