% Discrete All-Pole (DAP) envelope (AR model)
%
% Estimate an AR model based on given spectral peaks following the method in [1].
%
% Input
%  sins    : [Hz;amp] [2xN] Spectral peaks with frequency and linear amplitudes.
%            (as provided by sin_analysis.m)
%  fs      : [Hz] Sampling frequency
%  order   : Order of the AR model
% [dftlen] : DFT length of the envelope (if requested by the output)
% [opt]    : Additional options (see code below)
%
% Output
%  g  : gain factor
%  a  : AR coefficients
%  E  : Amplitude envelope (if dftlen is not empty)
%
% Reference
%  [1] El-Jaroudi, A., Makhoul, J.: Discrete All-Pole Modeling, IEEE Transactions
%      on Signal Processing 39(2), 411–423, 1991.
%
% Copyright (c) 2011 University of Crete - Computer Science Department
%
% License
%  This file is under the LGPL license,  you can
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
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function [g a E opt] = env_dap(sins, fs, order, dftlen, opt)

    if nargin<5
        opt.maxit     = 50;     % Maximum number of iterations allowed
        opt.alpha     = 0.5;    % Convergence speed
        opt.dISthresh = 10e-7;  % Stopping threshold for Itakura-Saito err diff
        opt.minbw     = [];     % [Hz] Minimum bandwidth of the poles.
                                %    Slow down the computation significantly !
                                %    It should be at least zero in order to ensure
                                %    the stability of the AR filter.
                                %    Advised to set it to 50Hz

        opt.debug = 0;          % 1: show only final solution
                                % 2: show each step
    end
    if nargin==0; g=opt; return; end

    if ~isempty(opt.minbw); minrr=bandwidth2rootradius(opt.minbw, fs); end

    wm = 2*pi*sins(1,:)'/fs; % Radian frequencies
    P = abs(sins(2,:));      % Magnitues
    N = size(sins,2);        % Number of Peaks

    B = exp(-j*wm*double(0:order));  % The Fourier basis
    Binv = exp(j*wm*double(0:order));% The inverse Fourier basis

    if opt.debug>1
        subplot(211);
            hold off;
            plot(fs*wm/(2*pi), log(P), 'xk');
            hold on;
    end

    % Compute the target AR from the power-spectrum P^2, eq. (5) in [1]
    % Surely a mistake in [1] because the autocorr is the inv Fourier transform
    % of the POWER spectrum and not that of the magnitude.
    R = (1/N)*real(P.^2*Binv); 
    Rmatinv = inv(toeplitz(R)); % Prepare the inverse of the AR matrix

    % Initial guess using the Levinson-Durbin recursion
    [a,e] = levinson(R,order);
    a = (1/sqrt(e))*a'; % Put the gain factor into the AR coefficients

    Phat = 1./abs(B*a); % Get the model amplitude response
    if opt.debug>1; subplot(211); plot(fs*wm/(2*pi), log(Phat), 'g'); end

    err = nan(1,opt.maxit);
    it=1;
    cont = true;
    while cont && it<opt.maxit
        if opt.debug>1
            subplot(211);
                plot(fs*wm/(2*pi), log(abs(Phat)));
                F = fs*(0:dftlen/2)/dftlen;
                E = gba2hspec(1, 1, a, dftlen);
                plot(F, log(abs(E)), 'r');
                xlim([0 fs/2]);
                xlabel('Frequency [Hz]');
                ylabel('Amplitude [log]');

            subplot(212);
                plot(log10(err));
                xlabel('Iterations');
                ylabel('Itakura-Saito error [log10]');
            pause
        end

        A = B*a; % Compute the frequency response of the model

        h = (1/N)*real(B.'*(1./A)); % Compute the reversed impulse response

        anew = a*(1-opt.alpha) + opt.alpha * (Rmatinv*h); % Update the AR coeffs

        % If asked, ensure stability of the AR filter
        % by limiting the bandwidth of the poles
        if ~isempty(opt.minbw); anew=polystab2(anew, minrr); end

        % Compute the Itakura-Saito error
        Phat = 1./abs(B*anew); % Get the model amplitude response
        err(it)  =  1/N *sum(P(:)./Phat-log(P(:)./Phat)-1); % eq. (14) in [1]

        if it>2 % Start at 3 in order to skip the influence of the initial guess
            cont = err(it-1)-err(it)>opt.dISthresh;
        end

        a = anew;
        
        it = it + 1;
    end

    % Put the gain back to a separate parameter in order to have a0=1
    g = 1./a(1);
    a = a./a(1);

    if nargout>2
        if isempty(dftlen); dftlen=4096; end
        E = gba2hspec(g, 1, a, dftlen);
    end

    if opt.debug==1
        hold off;
        plot(sins(1,:), log(sins(2,:)), 'xk');
        hold on;
        F = fs*(0:dftlen/2)/dftlen;
        E = gba2hspec(g, 1, a, dftlen);
        plot(F, log(abs(E)));
        xlim([0 fs/2]);
        keyboard
    end

return

