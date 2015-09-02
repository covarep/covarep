% Multi-Frame Analysis based on Discrete Cepstral Envelope (DCE-MFA)
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
%  Dk             : [log] Log energy corrections
%  af             : As in input, plus the aligned amplitudes values in the
%                   'a' field.
%  E              : The amplitude cepstral envelope
%  
% References
%  [1] Y. Shiga and S. King, "Estimation of voice source and vocal tract 
%      characteristics based on multi-frame analysis," EUROSPEECH, 2003.
%  [2] M. Campedel-Oudot, O. Cappe and E. Moulines, "Estimation of the Spectral
%      Envelope of Voiced Sounds Using a Penalized Likelihood Approach"
%
% Copyright (c) Yannis Stylianou, 2011, Bilbao
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
% This function is part of the COVAREP project: http://covarep.github.io/covarep
%
% Authors
%  Yannis Stylianou <yannis@csd.uoc.gr>
%  Gilles Degottex <degottex@csd.uoc.gr> (regularization term)
%

function [cc Dk af E] = env_dce_mfa(af, fs, order, extrap_dcny, scale, Bw, lr, dftlen)

    debug = 0; % 0:Do nothing; 1:Plot iterations info; 2:Plot results

    % Input parameters
    if nargin<4; extrap_dcny=true; end
    if nargin<5; scale = []; end
    if nargin<6 || isempty(Bw); Bw = 2000*(fs/16000); end% 2kHz with 16kHz[1]
    if nargin<7 || isempty(lr); lr = 0; end
    if nargin<8; dftlen=4096; end

    if ~isempty(scale)
        Bw = 100*fs; % To make it as similar as possible to the DCE-SFA
        lr = 3.5e-2; % To make it as similar as possible to the DCE-SFA [2]
        eval(['fnscale=@frq2' scale ';']);
        for n=1:numel(af)
            af(n).f = 0.5*fs*fnscale(af(n).f)/fnscale(fs/2);
        end
    end

    % If asked, extrapolate components at DC and up to Nyquist
    if extrap_dcny; af = env_extrap_sins_dcny(af, fs); end

    for n=1:numel(af)
        af(n).f = af(n).sins(1,:);
        af(n).a = log(af(n).sins(2,:));
    end

    if debug>1; subplot(211); hold off; end

    M = length(af);

    BWB = 0;  % [1](9)
    erI = 0;  % Error of previous step (initialization error)
    iter = 0; % Iteration number
    cc = zeros(order+1,1);  % initialization for ceps
    rt = 0;

    for k=1:M
        fk = af(k).f;
        ak = af(k).a;
        Nk = length(fk);
        fk = fk(:);
        ak = ak(:);
        Bk = [ones(Nk,1) 2*cos(2*pi*fk/fs*(1:order))];

        wk = exp(-fk.^2 / (2 * Bw * Bw))'; % Bottom-right paragraph p.3
        Wk = diag(wk)/Nk; % (8)

        % keep the main matrices
        BWB = BWB+(Bk'*Wk*Bk); % order by order
        af(k).Bk = Bk;        % keep Bk
        af(k).Wk = Wk;        % keep Wk
        uk = ones(Nk,1);

        dk = uk'*Wk*(ak)/(uk'*Wk*uk); % (10) with c=0
        erI = erI + ((ak-dk*uk)'*Wk*(ak-dk*uk));% (7) with c=0
        rt =  rt + (Bk'*Wk*(ak-dk*uk)); % Right-hand term of (9)
    end

    if debug>0; disp(['iter:' num2str(iter) ' error=' num2str(erI) 'log']); end

    % first estimation of ceps
    if lr==0
        cc = BWB\rt; % Solution of (9)
    else
        cc = (BWB+lr*diag(ones(size(BWB,1),1)))\rt; % Solution of (9) + Regul term
    end
    cc(1) = 0;

    iter = 1;
    while(1)
        er = 0;
        rt = 0; % Right-hand term of (9)
        Dk = zeros(M,1); % log energy corrections

        for k=1:M        
            ak = af(k).a;ak= ak(:);
            Bk = af(k).Bk;
            Wk = af(k).Wk;

            Nk = length(ak);
            uk = ones(Nk,1);

            dk = uk'*Wk*(ak-Bk*cc)/(uk'*Wk*uk);% (10)
            h = (ak-dk*uk-Bk*cc);
            er = er + (h'*Wk*h);% (7)
            rt =  rt + (Bk'*Wk*(ak-dk*uk)); % Right-hand term of (9)

            Dk(k) = dk;
        end
        if lr==0
            cc = BWB\rt; % Solution of (9)
        else
            cc = (BWB+lr*diag(ones(size(BWB,1),1)))\rt; % Solution of (9) + Regul term
        end
        cc(1) = 0;

        if debug>0; disp(['iter:' num2str(iter) ' error=' num2str(er) 'log']); end
    %      disp(['reldiff=' num2str(abs((erI-er)/er))]);

        if debug>1
            % plot
            subplot(211);
                plot(iter-1,log(erI), 'o');
                hold on
                plot(iter,log(er), 'x');
                title(num2str(iter));
            subplot(212);
                hold off;
                for k=1:M
                    plot(af(k).f,af(k).a-Dk(k),'+');
                    hold on;
                end
                dftlen = 2048;
                fv = fs*(0:dftlen/2)/dftlen;
                lF = 2*cos(2*pi*fv'/fs*(1:order))*cc(2:order+1);
                plot(fv,lF,'r', 'LineWidth', 2);
            keyboard
        end

        if( abs((erI-er)/er)<0.001 ) % Stop if error doesn't improve more than 0.1%
            break;
        else
            erI = er;
            iter = iter+1;
        end
    end


    % Align the gain corrections with respect to the central frame
    ci = floor((numel(af)-1)/2)+1;
    cc(1) = cc(1)+Dk(ci);
    Dk = Dk - Dk(ci);
    cc(2:end) = 2*cc(2:end);

    % Include the log energy corrections into the output
    for fi=1:numel(af)
        af(fi).a = af(fi).a - Dk(fi);
    end

    % If asked, compute the envelope
    if nargout>2
        if isempty(scale)
            E = exp(fft(cc, dftlen));
            E = E(1:end/2+1);
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

return
