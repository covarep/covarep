% Discrete Cepstral Envelope (DCE)
%
% Input
%  af            : Array of structures with matrices containing sinusoidal
%                  parameters (as the output of sin_analysis.m).
%                  Each matrix is made of:
%                   The 1st line for the frequency of each sinusoid [Hz]
%                   The 2nd line for their amplitude (linear scale)
%                  The DC has to be included (in the first column)
%  fs            : [Hz] Signal's sampling frequency
%  order         : Cepstral order
%  [extrap_dcny] : If true (default), alleviate stability problems of the DCE 
%                     by replacing the DCE value and extrapolating sinusoidal
%                     components up to Nyquist.
%                  If false, use the sinusoidal components as they are.  
%  [scale]       : empty  : Frequency linear scale (default)
%                 'mel'  : see frq2mel (for MFCC computation)
%                 'bark' : see frq2bark
%                 'erb'  : see frq2erb
%  [Bw]          : [Hz] Standard-deviation of the Gaussian used as weighting
%                       function (for emphasizing the importance of the low
%                       frequencies in the solution).
%                       Bw As to be big enough to have the weights still
%                       significant up to Nyquist.
%  [lr]          : Regularization parameter (as in [2]) (def. 0)
%  [dftlen]      : DFT's length, if the 4th output argument is requested.
%
% Output
%  cc             : Cepstral coefficients
%  E              : The amplitude cepstral envelope
%  
% References
%  [1] T. Galas and X. Rodet, "Generalized discrete cepstral analysis for
%      deconvolution of source-filter system with discrete spectra,"
%      in IEEE Applications of Signal Processing to Audio and Acoustics (ASSP)
%      Workshop, pp. 71-72, 1991.
%  [2] M. Campedel-Oudot, O. Cappe and E. Moulines, "Estimation of the Spectral
%      Envelope of Voiced Sounds Using a Penalized Likelihood Approach", IEEE
%      Transactions on Speech and Audio Processing, vol. 9, no. 5, pp. 469-481,
%      2001.
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
% This function is part of the COVAREP project: http://covarep.github.io/covarep
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function [cc E] = env_dce_sfa(af, fs, order, extrap_dcny, scale, Bw, lr, dftlen)

    % Input parameters
    if nargin<4; extrap_dcny=true; end
    if nargin<5; scale = []; end
    if nargin<6 || isempty(Bw); Bw = 2000*(fs/16000); end % To make it as
                                           % similar as possible to the DCE-MFA
    if nargin<7 || isempty(lr); lr = 0.035; end % [2]
    if nargin<8; dftlen=4096; end

    if ~isempty(scale)
        Bw = 100*fs;
        lr = 0.035; % To make it as similar as possible to [2]
        eval(['fnscale=@frq2' scale ';']);
        af.f = 0.5*fs*fnscale(af(n).f)/fnscale(fs/2);
    end

    % If asked, extrapolate components at DC and up to Nyquist
    if extrap_dcny; af = env_extrap_sins_dcny(af, fs); end

    af.f = af.sins(1,:);
    af.a = log(af.sins(2,:));


    % Prepare the weighting functions
    fk = af.f(:);
    Nk = length(fk);
    wk = exp(-fk.^2 / (2 * Bw * Bw))';
    wk = wk./sum(wk);
    Wk = diag(wk);

    % Compute the DCE envelope
    fk = af.f(:);
    ak = af.a(:);
    Nk = length(fk);
    Bk = [ones(Nk,1) cos(2*pi*fk/fs*(1:order))];

    BWB = (Bk'*Wk*Bk);
    rt =  (Bk'*Wk*(ak));

    if lr>0
        cc = (BWB+lr*diag(ones(size(BWB,1),1)))\rt;
    else
        cc = BWB\rt;
    end

    % If asked, compute the envelope
    if nargout>1
        if isempty(scale)
            E = fft(cc, dftlen);
            E = exp(E(1:end/2+1));
        else
            if strcmp(scale,'bark')
                E = barkcc2spec(cc, fs, dftlen);
                E = E(1:end/2+1);
            elseif strcmp(scale,'mel')
                E = exp(fft(cc, dftlen));
                E = E(1:dftlen/2+1);
                E = cc2hspec(cc, fs, dftlen);
                E = fwcep2hspec(cc, fs, dftlen);
            end
        end
    end

    % Plot final solution
    if 0
        hold off;

        global S Erefg selfr;
        if ~isempty(S);
            Sbins = fs*(0:dftlen/2)/dftlen;
            plot(Sbins, mag2db(abs(S)), ':k');
            hold on;
        else; hold off; end
        if ~isempty(Erefg);
            Erefgbins = 0.5*fs*(0:length(Erefg)-1)/length(Erefg);
            plot(Erefgbins, mag2db(abs(Erefg)), 'k');
            hold on;
        end

        plot(af.f, mag2db(exp(af.a)),'xr');
        hold on;
        F = fs*(0:dftlen/2)/dftlen;
%          plot(F, mag2db(abs(L)), 'g');
        hold on;
        plot(F, mag2db(abs(E)), 'b', 'LineWidth', 2);
        keyboard
    end

return
