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
%            Currently only even size is supported. If S is a vector, S
%            should have the size [1 m]. S can be matrix f of size [n m]
%  order   : Cepstral order
% [winlen] : Length of the window used to compute the spectrum
%            (Not yet implemented !, please use [])
% [mode]   : Binary flags encoded in a single value [def. 1+8+16+32]
%            1   : Smooth the cepstral coefficients with a hamming window
%                  To keep the same bandwith as the rectangular window
%                  the order has to be extend (factor 1.2).
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
%

function [E, cc, n] = env_te(S, order, winlen, mode, maxit, prec, presmooth_factor, gamma)

    if size(S,2)==1, S = S'; end

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

    dftlen = size(S,2);
    order = min(order, dftlen/2-1);

    if bitand(mode,1) % Use smoothing window
        order = round(1.2*order);
        order = min(order, dftlen/2-1);
        win = hamming(2*order+1)';
        win = win((end-1)/2+1:end);
    end

    A = fnlog(abs(S));
    cc = zeros(size(S,1),1+order);

    if bitand(mode,512)
        pE = env_te(S, round((order)/presmooth_factor), [], 1);
        slim = round(0.25*dftlen/order);
        pE = hspec2spec(pE);
        A(:,1:slim) = fnlog(abs(pE(:,1:slim)));
        A(:,end/2+1-slim+2:end/2+1) = fnlog(abs(pE(:,end/2+1-slim+2:end/2+1)));
        A = hspec2spec(A(:,1:end/2+1));
    end
    A0 = A;

    n=0;
    lV_final = nan(size(S,1),dftlen);
    idx_final = 1:size(S,1);
    while n<maxit

        cca = ifft(A,[],2,'symmetric');

        ccp = cca;
        ccp = [ccp(:,1), 2*ccp(:,2:1+order)];

        if bitand(mode,1)            % Use smoothing window
            Eo = sqrt(sum((2*cca(:,order+2:end/2)).^2,2));
            Ei = sqrt(sum(ccp(:,1:order+1).^2,2));
            lambda = (Ei+Eo)./Ei;
            cc = bsxfun(@times,lambda,win).*(ccp-cc) + cc; % Eq. (5) in [1]
        else
            cc = ccp;
        end

        lV = fft(cc, dftlen,2);

        A = max(A,real(lV));         % Max of log amplitudes

        max_diff = max(A0-real(lV),[],2); % Can create over-shot

        n = n+1;

        idx = max_diff<=dblim;
        if any(idx)
            % save finished lV
            lV_final(idx_final(idx),:) = lV(idx,:);
            
            % remove finished data
            idx_final = idx_final(~idx);
            A = A(~idx,:);
            cc = cc(~idx,:);
            A0 = A0(~idx,:);
        end
        
        % don't delete lV if we reach maxit
        if isempty(A)
            lV = [];
            break;
        end
    end

    if ~isempty(idx_final);
        lV_final(idx_final,:) = lV;
    end

    if mod(dftlen,2)==0,    lV_final=lV_final(:,1:end/2+1);
    else                    lV_final=lV_final(:,1:(end-1)/2+1);   end

    E = fnexp(lV_final);
end

function y = genlog(x, gamma)

    if gamma==0
        y = log(x);
    else
        y = (x.^gamma - 1)./gamma;
    end

end

function x = igenlog(y, gamma)

    if gamma==0
        x = exp(y);
    else
        x = exp(log(gamma*y+1)./gamma);
    end

end
