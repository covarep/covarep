function [ss,zo]=ssubmmse(si,fsz,pp)
%SSUBMMSE performs speech enhancement using mmse estimate of spectral amplitude or log amplitude [SS,ZO]=(S,FSZ,P)
%
% Inputs:
%   si      input speech signal
%   fsz     sample frequency in Hz
%           Alternatively, the input state from a previous call (see below)
%   pp      algorithm parameters [optional]
%
% Outputs:
%   ss      output enhanced speech (length is rounded down to the nearest frame boundary)
%   zo      output state
%
% The algorithm operation is controlled by a small number of parameters:
%
%        pp.of          % overlap factor = (fft length)/(frame increment) [2]
%        pp.ti          % desired frame increment [0.016 seconds]
%        pp.ri          % set to 1 to round ti to the nearest power of 2 samples [0]
%        pp.ta          % Time const for smoothing SNR estimate [0.396 seconds]
%        pp.gx          % maximum posterior SNR as a power ratio [1000]
%        pp.gn          % min posterior SNR as a power ratio when estimating prior SNR [1]
%        pp.xn          % minimum prior SNR [0]
%        pp.xb          % bias compensation factor for prior SNR [1]
%        pp.lg          % optimal estimate for log spectrum rather than spectrum [1]
%        pp.bt          % threshold for binary gain or -1 for continuous gain [-1]
%        pp.mx          % input mixture gain [0]
%
% The applied gain is mx+(1-mx)*optgain where optgain is calculated according to [1] or [2].
% If pp.bt>=0 then optgain is first thresholded with pp.bt to produce a binary gain 0 or 1.
%
% The default parameters implement the original algorithm in [1,2].
%
% Several parameters relate to the estimation of xi, the so-called "prior SNR",
% 
%             xi=max(a*pp.xb*xu+(1-a)*max(gami-1,pp.gn-1),pp.xn);
%
% This is estimated as a smoothed version of gami, the "posterior SNR"
% which is the noisy speech power divided by the noise power. This is
% clipped to (pp.gn-1), smoothed using a factor "a" which corresponds to a
% time-constasnt of pp.ta and then clipped to a minimum of pp.xn. The
% previous value is taken to be pp.xb*xu where xu is the ratio of the
% estimated speech amplitude squared to the noise power.
% 
% In addition it is possible to specify parameters for the noise estimation algorithm
% which implements reference [3] from which equation numbers are given in parentheses.
% They are as follows:
%
%        pp.taca      % (11): smoothing time constant for alpha_c [0.0449 seconds]
%        pp.tamax     % (3): max smoothing time constant [0.392 seconds]
%        pp.taminh    % (3): min smoothing time constant (upper limit) [0.0133 seconds]
%        pp.tpfall    % (12): time constant for P to fall [0.064 seconds]
%        pp.tbmax     % (20): max smoothing time constant [0.0717 seconds]
%        pp.qeqmin    % (23): minimum value of Qeq [2]
%        pp.qeqmax    % max value of Qeq per frame [14]
%        pp.av        % (23)+13 lines: fudge factor for bc calculation  [2.12]
%        pp.td        % time to take minimum over [1.536 seconds]
%        pp.nu        % number of subwindows to use [3]
%        pp.qith      % Q-inverse thresholds to select maximum noise slope [0.03 0.05 0.06 Inf ]
%        pp.nsmdb     % corresponding noise slope thresholds in dB/second   [47 31.4 15.7 4.1]
%
%
% If convenient, you can call specsub in chunks of arbitrary size. Thus the following are equivalent:
%
%                   (a) y=ssubmmse(s,fs);
%
%                   (b) [y1,z]=ssubmmse(s(1:1000),fs);
%                       [y2,z]=ssubmmse(s(1001:2000),z);
%                       y3=ssubmmse(s(2001:end),z);
%                       y=[y1; y2; y3];
%
% Note that in all cases the number of output samples will be less than the number of input samples if
% the input length is not an exact number of frames. In addition, if two output arguments
% are specified, the last partial frame samples will be retained for overlap adding
% with the output from the next call to ssubmmse().
%
% See also specsub() for an alternative gain function
%
% Refs:
%    [1] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 32(6):1109-1121, Dec 1984
%    [2] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum mean-square error log-spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 33(2):443-445, Apr 1985
%    [3] Rainer Martin.
%        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.
%    [4] O. Cappe.
%        Elimination of the musical noise phenomenon with the ephraim and malah noise suppressor.
%        IEEE Trans Speech Audio Processing, 2 (2): 345–349, Apr. 1994. doi: 10.1109/89.279283.
%    [5] J. Erkelens, J. Jensen, and R. Heusdens.
%        A data-driven approach to optimizing spectral speech enhancement methods for various error criteria.
%        Speech Communication, 49: 530–541, 2007. doi: 10.1016/j.specom.2006.06.012.
%    [6] R. Martin.
%        Statistical methods for the enhancement of noisy speech.
%        In J. Benesty, S. Makino, and J. Chen, editors,
%        Speech Enhancement, chapter 3, pages 43–64. Springer-Verlag, 2005.

%      Copyright (C) Mike Brookes 2004-2011
%      Version: $Id: ssubmmse.m,v 1.5 2011/07/18 07:40:53 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isstruct(fsz)
    fs=fsz.fs;
    qq=fsz.qq;
    qp=fsz.qp;
    ze=fsz.ze;
    s=zeros(length(fsz.si)+length(si(:)),1); % allocate space for speech
    s(1:length(fsz.si))=fsz.si;
    s(length(fsz.si)+1:end)=si(:);
else
    fs=fsz;     % sample frequency
    s=si(:);
    % default algorithm constants

    qq.of=2;        % overlap factor = (fft length)/(frame increment)
    qq.ti=16e-3;    % desired frame increment (16 ms)
    qq.ri=0;        % round ni to the nearest power of 2
    qq.ta=0.396;    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [1]
    qq.gx=1000;     % maximum posterior SNR = 30dB
    qq.gn=1;        % min posterior SNR as a power ratio when estimating prior SNR [1]
    qq.xn=0;        % minimum prior SNR = -Inf dB
    qq.xb=1;        % bias compensation factor for prior SNR [1]
    qq.lg=1;        % use log-domain estimator by default
    qq.bt=-1;       % suppress binary masking
    qq.mx=0;        % no input mixing
    if nargin>=3 && ~isempty(pp)
        qp=pp;      % save for estnoisem call
        qqn=fieldnames(qq);
        for i=1:length(qqn)
            if isfield(pp,qqn{i})
                qq.(qqn{i})=pp.(qqn{i});
            end
        end
    else
        qp=struct;  % make an empty structure
    end
end
% derived algorithm constants
if qq.ri
    ni=pow2(nextpow2(ti*fs*sqrt(0.5)));
else
    ni=round(qq.ti*fs);    % frame increment in samples
end
tinc=ni/fs;          % true frame increment time
a=exp(-tinc/qq.ta); % SNR smoothing coefficient
gx=qq.gx;           % max posterior SNR = 20 dB
kk=sqrt(2*pi);   % sqrt(8)*Gamma(1.5) - required constant
xn=qq.xn;           % floor for prior SNR, xi
gn1=max(qq.gn-1,0);        % floor for posterior SNR when estimating prior SNR
xb=qq.xb;

% calculate power spectrum in frames

no=round(qq.of);                                   % integer overlap factor
nf=ni*no;           % fft length
w=sqrt(hamming(nf+1))'; w(end)=[]; % for now always use sqrt hamming window
w=w/sqrt(sum(w(1:ni:nf).^2));       % normalize to give overall gain of 1
y=enframe(s,w,ni);
yf=rfft(y,nf,2);
yp=yf.*conj(yf);        % power spectrum of input speech
[nr,nf2]=size(yp);              % number of frames
if isstruct(fsz)
    [dp,ze]=estnoisem(yp,ze);   % estimate the noise using minimum statistics
    ssv=fsz.ssv;
    xu=fsz.xu;                  % saved unsmoothed SNR
else
    [dp,ze]=estnoisem(yp,tinc,qp);   % estimate the noise using minimum statistics
    ssv=zeros(ni*(no-1),1);             % dummy saved overlap
    xu=1;               % dummy unsmoothed SNR from previous frame
end
if ~nr                                  % no data frames
    ss=[];
else
    gam=min(yp./dp,gx);               % gamma = posterior SNR
    g=zeros(nr,nf2);   % create space for gain matrix
    if qq.lg            % use log domain estimator
        for i=1:nr
            gami=gam(i,:);
            xi=max(a*xb*xu+(1-a)*max(gami-1,gn1),xn);
            xir=xi./(1+xi);
            gi=xir.*exp(0.5*expint(xir.*gami));
            g(i,:)=gi;                  % save gain for later
            xu=gami.*gi.^2;         % unsmoothed prior SNR
        end
    else
        for i=1:nr
            gami=gam(i,:);
            xi=max(a*xb*xu+(1-a)*max(gami-1,gn1),xn);
            v=0.5*xi.*gami./(1+xi);     % note that this is 0.5*vk in [1]
            gi=(0.277+2*v)./gami; % accurate to 0.02 dB for v>0.5
            mv=v<0.5;
            if any(mv)
                vmv=v(mv);
                gi(mv)=kk*sqrt(vmv).*((0.5+vmv).*besseli(0,vmv)+vmv.*besseli(1,vmv))./(gami(mv).*exp(vmv));
            end
            g(i,:)=gi;                  % save gain for later
            xu=gami.*gi.^2;         % unsmoothed prior SNR
        end
    end
    if qq.bt>=0
        g=g>qq.bt;
    end
    g=qq.mx+(1-qq.mx)*g;   % mix in some of the input
    se=(irfft((yf.*g).',nf).').*repmat(w,nr,1);   % inverse dft and apply output window
    ss=zeros(ni*(nr+no-1),no);                      % space for overlapped output speech
    ss(1:ni*(no-1),end)=ssv;
    for i=1:no
        nm=nf*(1+floor((nr-i)/no));  % number of samples in this set
        ss(1+(i-1)*ni:nm+(i-1)*ni,i)=reshape(se(i:no:nr,:)',nm,1);
    end
    ss=sum(ss,2);
end
if nargout>1
    if nr
        zo.ssv=ss(end-ni*(no-1)+1:end);    % save the output tail for next time
        ss(end-ni*(no-1)+1:end)=[];
    else
        zo.ssv=ssv;  %
    end
    zo.si=s(length(ss)+1:end);      % save the tail end of the input speech signal
    zo.fs=fs;                       % save sample frequency
    zo.qq=qq;                       % save loval parameters
    zo.qp=qp;                       % save estnoisem parameters
    zo.ze=ze;                       % save state of noise estimation
    zo.xu=xu;
end
if ~nargout && nr>0
    ttax=(1:nr)*tinc;
    ffax=(0:size(g,2)-1)*fs/nf/1000;    ax=zeros(4,1);
    ax(1)=subplot(223);
    imagesc(ttax,ffax,20*log10(g)');
    colorbar;
    axis('xy');
    title(sprintf('Filter Gain (dB): ta=%.2g',qq.ta));
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');

    ax(2)=subplot(222);
    imagesc(ttax,ffax,10*log10(yp)');
    colorbar;
    axis('xy');
    title('Noisy Speech (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');

    ax(3)=subplot(224);
    imagesc(ttax,ffax,10*log10(yp.*g.^2)');
    colorbar;
    axis('xy');
    title('Enhanced Speech (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');

    ax(4)=subplot(221);
    imagesc(ttax,ffax,10*log10(dp)');
    colorbar;
    axis('xy');
    title('Noise Estimate (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');
    linkaxes(ax);
end