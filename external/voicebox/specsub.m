function [ss,zo]=specsub(si,fsz,pp)
%SPECSUB performs speech enhancement using spectral subtraction [SS,ZO]=(S,FSZ,P)
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
%        pp.g           % subtraction domain: 1=magnitude, 2=power [1]
%        pp.e           % gain exponent [1]
%        pp.am          % max oversubtraction factor [3]
%        pp.b           % max noise attenutaion in power domain [0.01]
%        pp.al          % SNR for oversubtraction=am (set this to Inf for fixed a) [-5 dB]
%        pp.ah          % SNR for oversubtraction=1 [20 dB]
%        pp.bt          % threshold for binary gain or -1 for continuous gain [-1]
%        pp.mx          % input mixture gain [0]
%        pp.gh          % maximum gain for noise floor [1]
%
% Following [1], the magnitude-domain gain in each time-frequency bin is given by
%                          gain=mx+(1-mx)*max((1-(a*N/X)^(g/2))^(e/g),min(gh,(b*N/X)^(e/2)))
% where N and X are the powers of the noise and noisy speech respectively.
% The oversubtraction factor varies linearly between a=am for a frame SNR of al down to
% a=1 for a frame SNR of ah. To obtain a fixed value of a for all values of SNR, set al=Inf.
% Common exponent combinations are:
%                      g=1  e=1    Magnitude Domain spectral subtraction
%                      g=2  e=1    Power Domain spectral subtraction
%                      g=2  e=2    Wiener filtering
% Many authors use the parameters alpha=a^(g/2), beta=b^(g/2) and gamma2=e/g instead of a, b and e
% but this increases interdependence amongst the parameters.
% If bt>=0 then the max(...) expression above is thresholded to become 0 or 1.
%
% In addition it is possible to specify parameters for the noise estimation algorithm
% which implements reference [2] from which equation numbers are given in parentheses.
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
%                   (a) y=specsub(s,fs);
%
%                   (b) [y1,z]=specsub(s(1:1000),fs);
%                       [y2,z]=specsub(s(1001:2000),z);
%                       y3=specsub(s(2001:end),z);
%                       y=[y1; y2; y3];
%
% Note that the number of output samples will be less than the number of input samples if
% the input length is not an exact number of frames. In addition, if two output arguments
% are specified, the last partial frame samples will be retained for overlap adding
% with the output from the next call to specsub().
%
% See also ssubmmse() for an alternative gain function
%
% Refs:
%    [1] M. Berouti, R. Schwartz and J. Makhoul
%        Enhancement of speech corrupted by acoustic noise
%        Proc IEEE ICASSP, 1979, 4, 208-211
%    [2] Rainer Martin.
%        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.

%      Copyright (C) Mike Brookes 2004
%      Version: $Id: specsub.m,v 1.4 2010/11/12 14:30:44 dmb Exp $
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

    qq.of=2;   % overlap factor = (fft length)/(frame increment)
    qq.ti=16e-3;   % desired frame increment (16 ms)
    qq.ri=0;       % round ni to the nearest power of 2
    qq.g=1;        % subtraction domain: 1=magnitude, 2=power
    qq.e=1;        % gain exponent
    qq.am=3;      % max oversubtraction factor
    qq.b=0.01;      % noise floor
    qq.al=-5;       % SNR for maximum a (set to Inf for fixed a)
    qq.ah=20;       % SNR for minimum a
    qq.bt=-1;       % suppress binary masking
    qq.mx=0;        % no input mixing
    qq.gh=1;        % maximum gain
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
    ni=pow2(nextpow2(qq.ti*fs*sqrt(0.5)));
else
    ni=round(qq.ti*fs);    % frame increment in samples
end
tinc=ni/fs;          % true frame increment time

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
else
    [dp,ze]=estnoisem(yp,tinc,qp);   % estimate the noise using minimum statistics
    ssv=zeros(ni*(no-1),1);             % dummy saved overlap
end
if ~nr                                  % no data frames
    ss=[];
else

    mz=yp==0;   %  mask for zero power time-frequency bins (unlikely)
    if qq.al<Inf
        ypf=sum(yp,2);
        dpf=sum(dp,2);
        mzf=dpf==0;     % zero noise frames = very high SNR
        af=1+(qq.am-1)*(min(max(10*log10(ypf./(dpf+mzf)),qq.al),qq.ah)-qq.ah)/(qq.al-qq.ah);
        af(mzf)=1;      % fix the zero noise frames
    else
        af=repmat(qq.am,nr,1);
    end
    switch qq.g
        case 1   % magnitude domain subtraction
            v=sqrt(dp./(yp+mz));
            af=sqrt(af);
            bf=sqrt(qq.b);
        case 2   % power domain subtraction
            v=dp./(yp+mz);
            bf=qq.b;
        otherwise % arbitrary subtraction domain
            v=(dp./(yp+mz)).^(0.5*qq.g);
            af=af.^(0.5*qq.g);
            bf=qq.b^(0.5*qq.g);
    end
    af =repmat(af,1,nf2);       % replicate frame oversubtraction factors for each frequency
    mf=v>=(af+bf).^(-1);        % mask for noise floor limiting
    g=zeros(size(v));           % reserve space for gain matrix
    eg=qq.e/qq.g;               % gain exponent relative to subtraction domain
    gh=qq.gh;
    switch eg
        case 1                          % Normal case
            g(mf)=min(bf*v(mf),gh);      % never give a gain > 1
            g(~mf)=1-af(~mf).*v(~mf);
        case 0.5                        
            g(mf)=min(sqrt(bf*v(mf)),gh);
            g(~mf)=sqrt(1-af(~mf).*v(~mf));
        otherwise
            g(mf)=min((bf*v(mf)).^eg,gh);
            g(~mf)=(1-af(~mf).*v(~mf)).^eg;
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
end
if ~nargout && nr>0
    ttax=(1:nr)*tinc;
    ffax=(0:size(g,2)-1)*fs/nf/1000;    ax=zeros(4,1);
    ax(1)=subplot(223);
    imagesc(ttax,ffax,20*log10(g)');
    colorbar;
    axis('xy');
    if qq.al==Inf
        title(sprintf('Filter Gain (dB): a=%.2g, b=%.3g',qq.am,qq.b));
    else
        title(sprintf('Filter Gain (dB): a=%.2g (%.0f to %.0fdB), b=%.3g',qq.am,qq.al,qq.ah,qq.b));
    end
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
    title(sprintf('Enhanced Speech (dB): g=%.2g, e=%.3g',qq.g,qq.e));
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