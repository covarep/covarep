% Feature extraction for the VAD system and Robust DNN-based ASR.
%
% Description
%  This is the module of feature extraction for the VAD system.
%
% Inputs
%  wave            : [samples] [Nx1] input signal (speech signal)
%  Fs              : [Hz]      [1x1] sampling frequency
%   
% Outputs
%  MatFeat      :  Matrix containg the feature vectors. 
%  t            : [seconds] Analysis instants.
%
% Example
%  Please see the "HOWTO_VAD.m" example file, or the "VAD_Drugman.m" function
%
% References
%  [1] T.Drugman, Y. Stylianou, Y. Kida, M. Akamine: "Voice Activity Detection:
%  Merging Source and Filter-based Information", IEEE Signal Processing Letters,
%  2014.
%  [2] S.O. Sadjadi, J. Hansen: "Unsupervised Speech Activity Detection Using
%  Voicing Measures and Perceptual Spectral Flux", IEEE Sig. Pro. Letters,
%  vol. 20, pp. 197-200, 2013.
%  [3] T. Drugman, Y. Stylianou, L. Chen, X. Chen, M. Gales, Robust
%  Excitation-based Features for Automatic Speech Recognition, IEEE 
%  International Conference on Acoustics, Speech and Signal Processing (ICASSP),
%   2015: http://tcts.fpms.ac.be/~drugman/files/ICASSP15_ASR.pdf
%
% Copyright (c) 2014 Toshiba Cambridge Research Laboratory
%
% License
%  This code will be part of the GLOAT toolbox
%    (http://tcts.fpms.ac.be/~drugman/Toolbox/)
%  with the following licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is also be part of the Covarep project:
%  http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be


function [MatFeat,t] = VAD_Drugman_FeatureExtraction(wave,Fs)


% Consider only signals sampled at 16kHz
if Fs~=16000
    wave=resample(wave,16000,Fs);
    Fs=16000;
end

% Windows of 30ms spaced every 10ms.
WinL=30/1000*Fs;
Hopsize=10/1000*Fs;
Win=hanning(WinL);
CorrWin=xcorr(Win,Win);
                                
Start=1;
Stop=WinL;
Ind=1;

display('Extracting Sadjadi''s features')
pause(0.0001)

% Extract the harmonicity, clarity, LP error and Harmonic Product Spectrum
% features as described in Sadjadi's paper [2]

Harm=[];
Clarity=[];
LPerr=[];
HPS=[];
t=[];

while Stop<length(wave)
    
    Seg=wave(Start:Stop);
    Seg=Seg.*Win;
    
    Rxx=xcorr(Seg,Seg)./CorrWin;
    Rxx(1:length(Seg)-1)=[];
    
    Lmin=round(2/1000*Fs);
    Lmax=round(16/1000*Fs);
    Vec=Rxx(Lmin:Lmax);
    
    [MaxiTmp,PosiTmp]=max(Vec);
    Harm(Ind)=log10(abs(MaxiTmp/(Rxx(1)-MaxiTmp)));
    
    Dxx=0.8*sqrt(2*(abs(Rxx(1)-Rxx)));
    Vec=Dxx(Lmin:Lmax);
    
    MiniTmp=min(Vec);
    MaxiTmp=max(Vec);
    Clarity(Ind)=1-MiniTmp/MaxiTmp;
    
    [a,e]=lpc(Seg,12);
    LPerr(Ind)=log(Rxx(1)/e);
    
    % It is not clear in Sadjadi's paper whether this has to be normalized
    % or not. I tried to contact the author, but did not get any reply.
    %                 Seg=Seg/sqrt(sum(Seg.^2));
    
    XSpec=fft(Seg,2048);
    XSpec=20*log10(abs(XSpec(1:1024)));
    
    Lmin=round(62.5/Fs*2048);
    Lmax=round(500/Fs*2048);
    
    VecTmp=[];
    for kk=Lmin:Lmax
        S=0;
        for hrm=1:8
            S=S+XSpec(hrm*kk);
        end
        VecTmp(kk-Lmin+1)=S;
    end
    
    HPS(Ind)=max(VecTmp);
    
    t(Ind)=round((Start+Stop)/2)/Fs;
    
    %%%%%%%%%%%%
    Start=Start+Hopsize;
    Stop=Stop+Hopsize;
    Ind=Ind+1;
    
end

display('Extracting MFCCs')
pause(0.0001)
[MatMFCC] = VAD_MFCC(wave,Fs);

display('Extracting the new features')
% Extract the features proposed in [1]


% CPP is computed on a low-passed version of the signal, as the harmonic
% part of the spectrum rarely goes beyond that limit [3]
wave2=resample(wave,8000,16000);
CPP = cpp(wave2,8000);
CPP=CPP(:,1)';
% This is just to have synchronous streams
CPP=[CPP(1) CPP(1) CPP(1) CPP CPP(end) CPP(end) CPP(end)];
Delta=length(CPP)-length(Harm);
if Delta>0
    CPP(end-Delta+1:end)=[];
end

% Computing the residual signal
res = lpcresidual(wave,25/1000*Fs,5/1000*Fs,12);

% Extract both SRH-based periodicity measures
[SRHVal1,SRHVal2] = VAD_SRH(res,Fs,80,350);

% This is just to have synchronous streams
SRHVal1(end-1:end)=[];
SRHVal2(end-1:end)=[];
Delta=length(SRHVal1)-length(Harm);
if Delta>0
    SRHVal1(1:Delta)=[];
    SRHVal2(1:Delta)=[];
end

MatFeat=[MatMFCC' Harm' Clarity' LPerr' HPS' CPP' SRHVal1' SRHVal2']';