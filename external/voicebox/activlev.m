function [lev,af,fso,vad]=activlev(sp,fs,mode)
%ACTIVLEV Measure active speech level as in ITU-T P.56 [LEV,AF,FSO]=(sp,FS,MODE)
%
%Usage: (1) lev=activlev(s,fs);     % speech level in units of power
%       (2) db=activlev(s,fs,'d');  % speech level in dB
%       (3) s=activlev(s,fs,'n');   % normalize active level to 0 dB
%
%Inputs: sp is the speech signal (with better than 20dB SNR)
%        FS is the sample frequency in Hz (see also FSO below)
%        MODE is a combination of the following:
%            r - raw omit input filters (default is 200 Hz to 5.5 kHz)
%            0 - no high pass filter (i.e. include DC)
%            4 - high pass filter at 40 Hz (but allows mains hum to pass)
%            1 - use cheybyshev 1 filter
%            2 - use chebyshev 2 filter (default)
%            e - use elliptic filter
%            h - omit low pass filter at 5.5 kHz
%            d - give outputs in dB rather than power
%            n - output a normalized speech signal as the first argument
%            N - output a normalized filtered speech signal as the first argument
%            l - give both active and long-term power levels
%            a - include A-weighting filter
%            i - include ITU-R-BS.468/ITU-T-J.16 weighting filter

%Outputs:
%    If the "n" option is specified, a speech signal normalized to 0dB will be given as
%    the first output followed by the other outputs.
%        LEV gives the speech level in units of power (or dB if mode='d')
%            if mode='l' is specified, LEV is a row vector with the "long term
%            level" as its second element (this is just the mean power)
%        AF  is the activity factor (or duty cycle) in the range 0 to 1
%        FSO is a column vector of intermediate information that allows
%            you to process a speech signal in chunks. Thus:
%
%            fso=fs; for i=1:inc:nsamp, [lev,fso]=activlev(sp(i:i+inc-1),fso,mode); end
%
%            is equivalent to:                lev=activlev(sp(1:nsamp),fs,mode)
%
%            but is much slower. The two methods will not give identical results
%            because they will use slightly different thresholds.
%        VAD is a boolean vector the same length as sp that acts as an approximate voice activity detector

%For completeness we list here the contents of the FSO structure:
%
%   ffs : sample frequency
%   fmd : mode string
%    nh : hangover time in samples
%    ae : smoothing filter coefs
%    bl : 200Hz HP filter numerator
%    al : 200Hz HP filter denominator
%    bh : 5.5kHz LP filter numerator
%    ah : 5.5kHz LP filter denominator
%    ze : smoothing filter state
%    zl : 200Hz HP filter state
%    zh : 5.5kHz LP filter state
%    zx : hangover max filter state
%  emax : maximum envelope exponent + 1
%   ssq : signal sum of squares
%    ns : number of signal samples
%    ss : sum of speech samples (not actually used here)
%    kc : cumulative occupancy counts
%    aw : weighting filter denominator
%    bw : weighting filter numerator
%    zw : weighting filter state

%      Copyright (C) Mike Brookes 2008-2012
%      Version: $Id: activlev.m 2516 2012-11-19 08:19:45Z dmb $
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

persistent nbin thresh c25zp c15zp e5zp
if isempty(nbin)
    nbin=20;    % 60 dB range at 3dB per bin
    thresh=15.9;    % threshold in dB
    % High pass s-domain zeros and poles of filters with passband ripple<0.25dB, stopband<-50dB, w0=1
    %    w0=fzero(@ch2,0.5); [c2z,c2p,k]=cheby2(5,50,w0,'high','s');
    %    function v=ch2(w); [c2z,c2p,k]=cheby2(5,50,w,'high','s'); v= 20*log10(prod(abs(1i-c2z))/prod(abs(1i-c2p)))+0.25;
    c25zp=[0.37843443673309i 0.23388534441447i; -0.20640255179496+0.73942185906851i -0.54036889596392+0.45698784092898i];
    c25zp=[[0; -0.66793268833792] c25zp conj(c25zp)];
    %       [c1z,c1p,c1k] = cheby1(5,0.25,1,'high','s');
    c15zp=[-0.659002835294875+1.195798636925079i -0.123261821596263+0.947463030958881i];
    c15zp=[zeros(1,5); -2.288586431066945 c15zp conj(c15zp)];
    %      [ez,ep,ek] = ellip(5,0.25,50,1,'high','s')
    e5zp=[0.406667680649209i 0.613849362744881i; -0.538736390607201+1.130245082677107i -0.092723126159100+0.958193646330194i];
    e5zp=[[0; -1.964538608244084]  e5zp conj(e5zp)];
    %    w=linspace(0.2,2,100);
    %    figure(1); plot(w,20*log10(abs(freqs(real(poly(c15zp(1,:))),real(poly(c15zp(2,:))),w)))); title('Chebyshev 1');
    %    figure(2); plot(w,20*log10(abs(freqs(real(poly(c25zp(1,:))),real(poly(c25zp(2,:))),w)))); title('Chebyshev 2');
    %    figure(3); plot(w,20*log10(abs(freqs(real(poly(e5zp(1,:))),real(poly(e5zp(2,:))),w)))); title('Elliptic');
end

if ~isstruct(fs)                        % no state vector given
    if nargin<3
        mode=' ';
    end
    fso.ffs=fs;                       	% sample frequency
    if any(mode=='r')                   % included for backward compatibility
        mode=['0h' mode];               % abolish both filters
    elseif fs<14000
        mode=['h' mode];               % abolish lowpass filter at low sample rates
    end
    fso.fmd=mode;                      % save mode flags
    ti=1/fs;
    g=exp(-ti/0.03);                    % pole position for envelope filter
    fso.ae=[1 -2*g g^2]/(1-g)^2;        % envelope filter coefficients (DC gain = 1)
    fso.ze=zeros(2,1);
    fso.nh=ceil(0.2/ti)+1;              % hangover time in samples
    fso.zx=-Inf;                        % initial value for maxfilt()
    fso.emax=-Inf;                      % maximum exponent
    fso.ns=0;
    fso.ssq=0;
    fso.ss=0;
    fso.kc=zeros(nbin,1);               % cumulative occupancy counts
    % s-plane zeros and poles of high pass 5'th order filter -0.25dB at w=1 and -50dB stopband
    if any(mode=='1')
        szp=c15zp;            % Chebyshev 1
    elseif any(mode=='e')
        szp=e5zp;             % Elliptic
    else
        szp=c25zp;            % Chebyshev 2
    end
    if all(mode~='0')
        if any(mode=='4')
            fl=40;               % 40 Hz cutoff
        else
            fl=200;              % 200 Hz cutoff
        end
        zl=2./(1-szp*tan(fl*pi/fs))-1;      % 200 Hz LF limit
        al=real(poly(zl(2,:)));              % high pass filter
        bl=real(poly(zl(1,:)));
        sw=1-2*rem(0:5,2).';
        fso.bl=bl*(al*sw)/(bl*sw);                  % scale to give HF gain of 1
        fso.al=al;
        fso.zl=zeros(5,1);                   % LF filter state
    end
    if all(mode~='h')
        zh=2./(szp/tan(5500*pi/fs)-1)+1;
        ah=real(poly(zh(2,:)));
        bh=real(poly(zh(1,:)));
        fso.bh=bh*sum(ah)/sum(bh);
        fso.ah=ah;
        fso.zh=zeros(5,1);
    end
    if any(mode=='a')
        [fso.bw fso.aw]=stdspectrum(2,'z',fs);
        fso.zw=zeros(length(fso.aw)-1,1);
    elseif any(mode=='i')
        [fso.bw fso.aw]=stdspectrum(8,'z',fs);
        fso.zw=zeros(length(fso.aw)-1,1);
    end
else
    fso=fs;             % use existing structure
end
md=fso.fmd;
ns=length(sp);
if ns                       % process this speech chunk
    % apply the input filters to the speech
    if all(md~='0')
        [sq,fso.zl]=filter(fso.bl,fso.al,sp(:),fso.zl);     % highpass filter
    else
        sq=sp(:);
    end
    if all(md~='h')
        [sq,fso.zh]=filter(fso.bh,fso.ah,sq(:),fso.zh);     % lowpass filter
    end
    if any(md=='a') || any(md=='i')
        [sq,fso.zw]=filter(fso.bw,fso.aw,sq(:),fso.zw);     % weighting filter
    end
    fso.ns=fso.ns+ns;      % count the number of speech samples
    fso.ss=fso.ss+sum(sq);  % sum of speech samples
    fso.ssq=fso.ssq+sum(sq.*sq);    % sum of squared speech samples
    [s,fso.ze]=filter(1,fso.ae,abs(sq(:)),fso.ze);     % envelope filter
    [qf,qe]=log2(s.^2);         % take efficient log2 function, 2^qe is upper limit of bin
    qe(qf==0)=-Inf;           % fix zero values
    [qe,qk,fso.zx]=maxfilt(qe,1,fso.nh,1,fso.zx);       % apply the 0.2 second hangover
    oemax=fso.emax;
    fso.emax=max(oemax,max(qe)+1);
    if fso.emax==-Inf
        fso.kc(1)=fso.kc(1)+ns;
    else
        qe=min(fso.emax-qe,nbin);   % force in the range 1:nbin. Bin k has 2^(emax-k-1)<=s^2<=2^(emax-k)
        wqe=ones(length(qe),1);
        % below: could use kc=cumsum(accumarray(qe,wqe,nbin)) but unsure about backwards compatibility
        kc=cumsum(full(sparse(qe,wqe,wqe,nbin,1)));    % cumulative occupancy counts
        esh=fso.emax-oemax;         % amount to shift down previous bin counts
        if esh<nbin-1               % if any of the previous bins are worth keeping
            kc(esh+1:nbin-1)=kc(esh+1:nbin-1)+fso.kc(1:nbin-esh-1);
            kc(nbin)=kc(nbin)+sum(fso.kc(nbin-esh:nbin));
        else
            kc(nbin)=kc(nbin)+sum(fso.kc); % otherwise just add all old counts into the last (lowest) bin
        end
        fso.kc=kc;
    end
end
if fso.ns                       % now calculate the output values
    if fso.ssq>0
        aj=10*log10(fso.ssq*(fso.kc).^(-1));
        % equivalent to cj=20*log10(sqrt(2).^(fso.emax-(1:nbin)-1));
        cj=10*log10(2)*(fso.emax-(1:nbin)-1);   % lower limit of bin j in dB
        mj=aj'-cj-thresh;
        %  jj=find(mj*sign(mj(1))<=0); % Find threshold
        jj=find(mj(1:end-1)<0 &  mj(2:end)>=0,1); % find +ve transition through threshold
        if isempty(jj)  % if we never cross the threshold
            if mj(end)<=0 % if we end up below if
                jj=length(mj)-1; % take the threshold to be the bottom of the last (lowest) bin
                jf=1;
            else          % if we are always above it
                jj=1;     % take the threshold to be the bottom of the first (highest) bin
                jf=0;
            end
        else
            jf=1/(1-mj(jj+1)/mj(jj));   % fractional part of j using linear interpolation
        end
        lev=aj(jj)+jf*(aj(jj+1)-aj(jj)); % active level in decibels
        lp=10.^(lev/10);
        if any(md=='d')
            lev=[lev 10*log10(fso.ssq/fso.ns)];
        else
            lev=[lp fso.ssq/fso.ns];
        end
        af=fso.ssq/(fso.ns*lp);
    else
        af=0;
        if all(md~='d')
            lev=[0 0];
        else
            lev=[-200 -200];
        end
    end
    if all(md~='l')
        lev=lev(1);         % only output the first element of lev normally
    end
end
if nargout>3
    vad=maxfilt(s,1,fso.nh,1);
    vad=vad>(sqrt(lp)/10^(thresh/20));
end
if ~nargout
    levdb=10*log10(lp);
    clf;
    subplot(211);
    plot((1:ns)/fso.ffs,[sp s (qe<=jj)*sqrt(lp)]);
    xlabel('Time (s)');
    title(sprintf('Active Level = %.2g dB, Activity = %.0f%% (ITU-T P.56)',levdb,100*af));
    ylabel('Amplitude');
    legend('Signal','Smoothed envelope','Active Level','Location','SouthEast');
    subplot(212);
    plot(cj,repmat(levdb,nbin,1),'k:',cj,[aj(:) cj(:)+thresh cj(:)]);
    xlabel('Threshold (dB)');
    ylabel('Active Level (dB)');
    legend('Active Level','Speech Level',sprintf('Threshold+%.1f dB',thresh),'Threshold','Location','SouthEast');
elseif any(md=='n') || any(md=='N')
    fsx=fso;
    fso=af;
    af=lev;
    if any(md=='n')
        sq=sp;
    end
    if fsx.ns>0 && fsx.ssq>0
        lev=sq/sqrt(lp);
    else
        lev=sq;
    end
end
