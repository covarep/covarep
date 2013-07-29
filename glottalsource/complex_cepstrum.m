% This script achieves an estimation of the glottal flow derivative based
% on the complex cesptrum.
%
% Octave compatible
%
% Description
%  The complex cepstrum-based technique, fully described in [1] and [2],
%  allows an estimation of the glottal flow derivative based on the
%  mixed-phase properties of the speech signal. For each voiced glottal
%  cycle, a GCI-synchronous analysis is performed and the maximum-phase
%  (i.e. anti-causal) component of speech is isolated. This component
%  corresponds to the glottal open-phase.
%
% Inputs
%  wave             : [samples] [Nx1] input signal (speech signal)
%  fs               : [Hz]      [1x1] sampling frequency
%  gci              : [s] vector containing the location of the estimated GCIs.%                     
%                     (you can use the "gci_sedreams" function included in this 
%                     toolbox for this).
%  F0s              : [Hz] vector containing the F0 estimates. It is required to
%                     use a frame shift of 10 ms.
%                     (you can use the "pitch_srh" function included in this 
%                     toolbox for this).
%  VUVDecisions     : vector containing the binary voicing decisions. It is
%                     required to use a frame shift of 10 ms.
%                     (you can use the "pitch_srh" function included in this 
%                     toolbox for this).
%
% Outputs
%  GlottalSource    : [samples] [Nx1] estimate of the glottal flow
%                     derivative using the complex cepstrum-based method.
%
%
% Example
%  And see the HOWTO_glottalsource.m example file.
%
% References
% [1] T.Drugman, B.Bozkurt, T.Dutoit, Complex Cepstrum-based Decomposition
%     of Speech for Glottal Source Estimation, Interspeech09, Brighton, U.K, 2009
%
% [2] T.Drugman, B.Bozkurt, T.Dutoit, Causal-anticausal Decomposition of
%     Speech using Complex Cepstrum for Glottal Source Estimation, Speech
%     Communication, Volume 53, Issue 6, pp. 855-866, July 2011.

%  Publications available at the following link:
%  http://tcts.fpms.ac.be/~drugman/files/SPECOM11-ComplexCepstrum.pdf
%  http://tcts.fpms.ac.be/~drugman/files/IS09-ComplexCepstrum.pdf
%
% Copyright (c) 2011 University of Mons, FNRS
%
% License
%  This code is a part of the GLOAT toolbox with the following
%  licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be
%

function [GlottalSource] = complex_cepstrum(wave,fs,gci,f0,VUVDecisions)

% Have a VUV decision and F0 estimation for each sample
VUVDecisions2=zeros(1,length(wave));
f02=zeros(1,length(wave));
HopSize=round(10/1000*fs);
for k=1:length(VUVDecisions)
    VUVDecisions2((k-1)*HopSize+1:k*HopSize)=VUVDecisions(k);
    f02((k-1)*HopSize+1:k*HopSize)=f0(k);
end

% Keep only GCIs for voiced segments
gci = round(gci*fs+1);
gindex = VUVDecisions2(gci)~=0;
gcis = gci(gindex);

% For each glottal cycle, apply a specific windowing satisfying the
% conditions explained in the paper mentioned above, and compute the
% complex cepstrum-based decomposition
GlottalSource=zeros(1,length(wave));

% Low-pass filtering because after separation, an irrelevant
% high-frequency noise (>4000Hz) may appear
Wp = 4000/(fs/2); Ws = 4100/(fs/2);
Rp = 3; Rs = 60;
[n,Wn] = ellipord(Wp,Ws,Rp,Rs);
[b_l,a_l] = ellip(n,Rp,Rs,Wn);

for k=1:length(gcis)
    
    T0=round(fs/f02(gcis(k)));
    
    Seg=wave(gcis(k)-round(0.9*T0):gcis(k)+round(0.9*T0));
    Seg=Seg.*blackman(length(Seg));
    
    [minPhase,maxPhase] = CCD_FrameLevel_GlottalFlowEstimation(Seg);
    [maxPhase] = FilterAndNormalize(maxPhase,b_l,a_l);
    
    maxPhase=maxPhase*sqrt(sum(Seg.^2)/sum(maxPhase.^2));
    
    GlottalSource(gcis(k)-length(maxPhase)+1:gcis(k))=GlottalSource(gcis(k)-length(maxPhase)+1:gcis(k))+maxPhase';    
end



function [sig] = FilterAndNormalize(sig,b_l,a_l)

sig=filtfilt(b_l,a_l,sig);

msig1 = max(sig);
msig2 = abs(min(sig));

if(msig1>msig2)
    sig=-sig;
    msig = msig1;
else
    msig = msig2;
end

sig=sig/(msig);


function [minPhase,maxPhase] = CCD_FrameLevel_GlottalFlowEstimation(inSignal)

halfLength=round(length(inSignal)/2);

[xhat]=ComputeComplexCepstrum(inSignal);   

antiCausalPortionCepstrum=xhat;
antiCausalPortionCepstrum(2:round(length(xhat)/2))=0;
antiCausalPortionCepstrum(1)=antiCausalPortionCepstrum(1)/2;
           
causalPortionCepstrum=xhat;
causalPortionCepstrum(round(length(xhat)/2)+1:end)=0;
causalPortionCepstrum(1)=causalPortionCepstrum(1)/2;

minPhase=inverse_CC(causalPortionCepstrum);
maxPhase=inverse_CC(antiCausalPortionCepstrum);
    
maxPhase=diff(maxPhase);
minPhase=filter([1 0],[1 -1],minPhase);

maxPhase=maxPhase(end-halfLength+1:end);
if(max(maxPhase)>max(-maxPhase))
    maxPhase=-maxPhase;    
end

minPhase=minPhase(1:length(inSignal)-length(maxPhase));
if(max(minPhase)>max(-minPhase))
    minPhase=-minPhase;    
end

maxPhase=maxPhase/(max(abs(maxPhase)));
minPhase=minPhase/(max(abs(minPhase)));



function [xhat] = ComputeComplexCepstrum(x)

n=2048;

h = fft(x,n);
[ah] = rcunwrap(angle(h));

logh = log(abs(h))+1i*ah;
xhat = real(ifft(logh));

function [y,nd] = rcunwrap(x)
%RCUNWRAP Phase unwrap utility used by CCEPS.
%   RCUNWRAP(X) unwraps the phase and removes phase corresponding
%   to integer lag.  See also: UNWRAP, CCEPS.

%   Author(s): L. Shure, 1988
%   	   L. Shure and help from PL, 3-30-92, revised

n = length(x);
y = unwrap(x);

nh = fix((n+1)/2);
y(:) = y(:)' - y(nh+1)*(0:(n-1))/nh;


function x = inverse_CC(xhat)

logh = fft(xhat);
h = exp(real(logh)+1i*imag(logh));
x = real(ifft(h));