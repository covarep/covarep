% Compute the "True-Envelope" (TE) (cepstral model)
%
%  This method estimate an amplitude cepstral envelope according to [1,3,4].
%  This implementation does NOT provide the complexity optimizations
%  proposed in [2].
%
%  Note that, unfortunately, the "true-envelope" is not the true envelope ...
%
%  The interface of this implementation is similar to that of Ircam's
%  implementation:
%  [env,coef,envphas,freqmel,freqlin]=mtltrueenv(ampspec,order,winlen
%                                     ,mode,maxit,prec,presmooth,freqbreak)
%  Please follow this interface for any extension of this code.
% 
% Inputs
%  S       : Spectrum to be enveloped (full spectrum as computed by fft.m)
%            Currently only even size is supported
%  order   : Cepstral order
% [winlen] : Length of the window used to compute the spectrum
%            (Not yet implemented !, please use [])
% [mode]   : Binary flags encoded in a single value [def. 1+8+16+32]
%            1   : Smooth the cepstral coefficients with a hamming window.
%            512 : Correct the DC and Nyquist bins such as they do not create
%                  ripples in the envelope. It is done by estimating the
%                  envelope with a lower order than given (see presmooth_factor)
%                  and replace the DC and Nyquist bins by the values of this 
%                  oversmoothed envelope.
%            (bits 8+16+32 are related to the complexity optimization
%             which has still to be implemented !)
%            (def. value is 1+512 (advised for speech))
% [maxit]  : Maximum number of iterations [def. 200]
% [prec]   : [dB] Precision of the envelope. [def. 2dB]
%            Maximum distance between the estimated envelope and the spectral 
%            peaks.
% [presmooth_factor]: Oversmoothing factor to correct the DC and Nyquist bins.
%             (see mode=512)
%
% Outputs
%  E        : Amplitude spectral envelope (half spectrum)
%  cc       : Cepstral coefficients
%  n        : Resulting number of iterations
%
% References
%  [1] A. Roebel, and X. Rodet, "Efficient Spectral Envelope Estimation And Its
%      Application To Pitch Shifting And Envelope Preservation", Proc. Digital
%      Audio Effects (DAFx), DAFX, 30-35, 2005
%  [2] A. Roebel, and X. Rodet, "Real Time Signal Transposition With Envelope
%      Preservation In The Phase Vocoder", ICMC, 2005
%  [3] A. Roebel, F. Villavicencio, X. Rodet, "On Cepstral and All-Pole based
%      Spectral Envelope Modeling with unknown Model order", Pattern Recognition
%      Letters, vol. 28, no 11,2007.
%  [4] S. Imai, "Cepstral analysis synthesis on the mel frequency scale", Proc.
%      IEEE International Conference on Acoustics, Speech, and Signal Processing
%      (ICASSP), vol 8, pp 93-96, 1983.
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

function [E, cc, n] = env_te(S, order, winlen, mode, maxit, prec, presmooth_factor, gamma)

    if nargin<3; winlen = []; end
    if ~isempty(winlen); disp('Use of winlen not yet implemented'); end
    if nargin<4; mode = 1+512; end
    if nargin<5; maxit = 200; end
    if nargin<6; prec = 2; end
    if nargin<7; presmooth_factor = 2; end
    if nargin<8;
        fnlog = @log;
        fnexp = @exp;
    else
        fnlog = @(x) genlog(x,gamma);
        fnexp = @(x) igenlog(x,gamma);
    end
    dblim = log(10.^(prec/20));

    if bitand(mode,1) % Use smoothing window
        order = round(1.2*order); % [1] 1.66, [matmtl] 1.2
        win = hamming(2*order+1)';
        win = win((end-1)/2+1:end);
    end

    dftlen = length(S);
    order = min(order, dftlen/2-1);

    A = fnlog(abs(S));
    cc = zeros(1,1+order);

    if bitand(mode,512)
        pE = env_te(S, round((order)/presmooth_factor), [], 1);
        slim = round(0.25*dftlen/order);
        pE = hspec2spec(pE);
        A(1:slim) = fnlog(abs(pE(1:slim)));
        A(end/2+1-slim+2:end/2+1) = fnlog(abs(pE(end/2+1-slim+2:end/2+1)));
        A = hspec2spec(A(1:end/2+1));
    end
    A0 = A;

    n=0;
    max_diff = Inf;
    while n<maxit && max_diff>dblim

        cca = ifft(A);

        ccp = cca;
        ccp = [ccp(1), 2*ccp(2:1+order)];

        if bitand(mode,1)            % Use smoothing window
            Eo = sqrt(sum((2*cca(order+2:end/2)).^2));
            Ei = sqrt(sum(ccp(1:order+1).^2));
            lambda = (Ei+Eo)/Ei;
            cc = lambda.*win.*(ccp-cc) + cc; % Eq. (5) in [1]
        else
            cc = ccp;
        end

        lV = fft(cc, dftlen);
        A = max(A,real(lV));         % Max of log amplitudes

        max_diff = max(A0-real(lV)); % Can create over-shot

        n = n+1;

        if 0
            V3clear();
            V3spec(fnexp(A0), 1, 'k');
            V3spec(fnexp(A), 1, 'g');
            V3spec(fnexp(lV), 1, 'b'); subplot(312);xlim([0 0.04]);
            keyboard
        end
    end

    if mod(dftlen,2)==0; lV=lV(1:end/2+1);
    else;                lV=lV(1:(end-1)/2+1);   end

    E = fnexp(lV);

return

function y = genlog(x, gamma)

    if gamma==0
        y = log(x);
    else
        y = (x.^gamma - 1)./gamma;
    end

return

function x = igenlog(y, gamma)

    if gamma==0
        x = exp(y);
    else
        x = exp(log(gamma*y+1)./gamma);
    end

return
