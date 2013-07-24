% RESKEW is an efficient technique to determine the polarity of the speech
% signal.
%
% Octave compatible
%
% Description
% The Residual Excitation Skewness (RESKEW) method is described in [1].
% This algorithm determines the polarity of the speech signal by inspecting
% the skewness of two residual excitation signals. Its advantages (shown in
% [1]) are: i) its high performance; ii) its robustness to an additive
% noise; iii) the fact that it does not depend any voicing or pitch
% estimate.
%
%
% Inputs
%  s               : [samples] [Nx1] input signal (speech signal)
%  fs              : [Hz]      [1x1] sampling frequency
%
% Outputs
%  polarity              : [1x1] the speech polarity
%
% Example
%  polarity = polarity_reskew(s,fs);
%
% References
%  [1] T.Drugman, Residual Excitation Skewness for Automatic Speech Polarity
%  Detection, IEEE Signal Processing Letters, vol. 20, issue 4, pp.
%  387-390, 2013.
%  Publication available at the following link:
%  http://tcts.fpms.ac.be/~drugman/files/SPL-Polarity.pdf
%
% Copyright (c) 2013 University of Mons, FNRS
%
% License
%  This code is a part of the Voice Analysis Toolkit with the following
%  licence:
%  The software product (and any modifications thereof) is distributed under 
%  a dual licence: an open source license for individual, non-commercial 
%  purposes, and a commercial license. The opensource licence under which 
%  the product is distributed is GNU GPL v2. For individual users, this 
%  licence suits their use as these are not distributing proprietary 
%  modifications, additions to, or derivatives of the product and don't 
%  require legal protection of a commercial licence. For commercial users, 
%  where open source does not meet their requirements, we offer commercial 
%  licensing of the product. A commercial license permits customers to 
%  modify, add or produce derivative products without the obligation of 
%  making the subsequent code open source. For more information regarding 
%  our commercial licence, please contact john.whelan@tcd.ie
%
% This function is part of the Common Speech Processing Repository (TODO)
% TODO URL
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be
%
% $Id <info set by the versioning system> $

function polarity = polarity_reskew(s,fs)

Wp = 480/(fs/2);
Ws = 500/(fs/2);
Rp = 3; Rs = 60;
[n,Wn] = ellipord(Wp,Ws,Rp,Rs);
[b_h,a_h] = ellip(n,Rp,Rs,Wn,'high');

s_h=filtfilt(b_h,a_h,s);


[res] = GetLPCresidual_WithTwoSignals(s,s_h,round(25/1000*fs),round(5/1000*fs),round(fs/1000+2));
[res2] = GetLPCresidual_WithTwoSignals(s,s,round(25/1000*fs),round(5/1000*fs),round(fs/1000+2));

Val1=skewness(res);
Val2=skewness(res2);

polarity=sign(Val2-Val1);



%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,LPCcoeff] = GetLPCresidual_WithTwoSignals(wave,waveToComputeLPC,L,shift,order)

% %%%
%  
% Use: [res] = GetLPCresidual(wave,L,shift,order,gci,type,t0)
% 
% 
% L=window length (samples) (typ.25ms)
% shift=window shift (samples) (typ.5ms)
% order= LPC order
% gci=gci position (samples)
% type=vector of voicing decisions (=0 if Unvoiced, =1 if Voiced)
% t0=vector of period values (in samples)
% 
% Written by Thomas Drugman, TCTS Lab.
% 
% %%%


start=1;
stop=start+L;

res=zeros(1,length(wave));
LPCcoeff=zeros(order+1,round(length(wave)/shift));
n=1;
while stop<length(wave)
    
    segment=waveToComputeLPC(start:stop);
    segment=segment.*hanning(length(segment));
    
    segment2=wave(start:stop);
    segment2=segment2.*hanning(length(segment2));    
    
    [A]=lpc(segment,order);
    LPCcoeff(:,n)=A(:);
    
    inv=filter(A,1,segment2);
    
    inv=inv*sqrt(sum(segment2.^2)/sum(inv.^2));
    
    res(start:stop)=res(start:stop)+inv';
    
    start=start+shift;
    stop=stop+shift;
    n=n+1;
end

res=res/max(abs(res));
