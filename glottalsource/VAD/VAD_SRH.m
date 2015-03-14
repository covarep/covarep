function [SRHVal1,SRHVal2] = VAD_SRH(res,Fs,F0min,F0max)

% Periodicity measures based on the SRH robust pitch tracker.
%
% Octave compatible
%
% Description
%  The Summation of the Residual Harmonics (SRH) method is described in [1].
%  This algorithm exploits a criterion taking into the strength of the
%  harmonics and subharmonics of the residual excitation signal in order to
%  determine both voicing decision and F0 estimates. It is shown in [1] to
%  be particularly interesting in adverse conditions (low SNRs with various
%  types of additive noises).
%
%
% Inputs
%  res            : [samples] [Nx1] residual signal (as obtained by the
%                   "lpcresidual.m" function
%  Fs              : [Hz]      [1x1] sampling frequency
%  F0min           : [Hz] [1x1] minimum possible F0 value
%  F0max           : [Hz]  [1x1] maximum possible F0 value
%
% Outputs
%  SRHVal1,SRHVal2   : vectors containing the SRH-based periodicity values.
%
% 
% References
%  [1] T.Drugman, A.Alwan, "Joint Robust Voicing Detection and Pitch Estimation
%      Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011.
%      Publication available at the following link:
%      http://tcts.fpms.ac.be/~drugman/files/IS11-Pitch.pdf
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
%  Thomas Drugman thomas.drugman@umons.ac.b


if Fs>16000
    res=resample(res,16000,Fs);
    Fs=16000;
end

Nfft=4096;

F0min=round(F0min/Fs*Nfft);
F0max=round(F0max/Fs*Nfft);

start=1;
stop=round(100/1000*Fs);
shift=round(10/1000*Fs);

sig=res';

Nframes=floor((length(sig)-stop)/shift) + 1;
SRHVal1 = zeros(1,Nframes);
SRHVal2 = zeros(1,Nframes);

BlackWin=blackman(stop-start+1);

index=1;
while stop<=length(sig)

    seg=sig(start:stop);
    seg=seg.*BlackWin;    
    seg=seg-mean(seg);

    Spec=fft(seg,Nfft);
    Spec=abs(Spec(1:Nfft/2));
    
    Spec2=Spec/sqrt(sum(Spec.^2));
        
    SRHs=zeros(1,F0max);    
    SRHs2=zeros(1,F0max);    
    
    % SRH spectral criterion
    for freq=F0min:F0max
        SRHs(freq)=(Spec(freq)+Spec(2*freq)+Spec(3*freq)+Spec(4*freq)+Spec(5*freq))-(Spec(round(1.5*freq))+Spec(round(2.5*freq))+Spec(round(3.5*freq))+Spec(round(4.5*freq)));
        SRHs2(freq)=(Spec2(freq)+Spec2(2*freq)+Spec2(3*freq)+Spec2(4*freq)+Spec2(5*freq))-(Spec2(round(1.5*freq))+Spec2(round(2.5*freq))+Spec2(round(3.5*freq))+Spec2(round(4.5*freq)));
    end
    
    [maxi,posi]=max(SRHs);
    F0frame=posi;    
    SRHVal2(index)=SRHs(F0frame);
    
    [maxi,posi]=max(SRHs2);
    F0frame=posi;    
    SRHVal1(index)=SRHs2(F0frame);
    
    start=start+shift;
    stop=stop+shift;
    index=index+1;
end

% Zero-padding necessary as the frames are 100ms-long
SRHVal1=[0 0 0 0 0 SRHVal1 0 0 0 0 0];
SRHVal2=[0 0 0 0 0 SRHVal2 0 0 0 0 0];
