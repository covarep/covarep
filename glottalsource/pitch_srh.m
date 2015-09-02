% SRH is a robust pitch tracker.
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
%  wave            : [samples] [Nx1] input signal (speech signal)
%  fs              : [Hz]      [1x1] sampling frequency
%  f0min           : [Hz]      [1x1] minimum possible F0 value
%  f0max           : [Hz]      [1x1] maximum possible F0 value
%  hopsize         : [ms]      [1x1] time interval between two consecutive
%                    frames (i.e. defines the rate of feature extraction).
%
% Outputs
%  F0s             : vector containing the F0 estimates (values are
%                    provided even in unvoiced parts).
%  VUVDecisions    : vector containing the binary voicing decisions.
%  SRHVal          : vector containing the SRH values (according the
%                    harmonic criterion - voicing decision are derived from
%                    these values by simple thresholding).
%  time            : [s] Analysis instants of the features described above.
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
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
%  Thomas Drugman thomas.drugman@umons.ac.be
%
% Modified
%  John Kane kanejo@tcd.ie September 27th 2014 - Bug fix and efficiency


function [F0s,VUVDecisions,SRHVal,time] = pitch_srh(wave,fs,f0min,f0max,hopsize)

if length(wave)/fs<0.1
    display('SRH error: the duration of your file should be at least 100ms long');
    pause(0.001)
    F0s=0;VUVDecisions=0;SRHVal=0;time=0;
end

if f0max<=f0min
    display('You look funny! Your f0min should be lower than f0max!!')
    pause(0.001)
end
    

if fs>16000
    display('Sample rate is bigger than 16kHz. Audio is resampled.')
    wave=resample(wave,16000,fs);
    fs=16000;
end

if nargin < 5
    hopsize=10;
end

%% Setings
nHarmonics = 5;
SRHiterThresh = 0.1;
SRHstdThresh = 0.05;
VoicingThresh = 0.07;
VoicingThresh2 = 0.085;
LPCorder=round(3/4*fs/1000);
Niter=2;

%% Compute LP residual
[res] = lpcresidual(wave,round(25/1000*fs),round(5/1000*fs), ...
                    LPCorder);

%% Create frame matrix
waveLen = length(wave);
clear wave;
frameDuration = round(100/1000*fs)-2; % Minus 2 to make equivalent
                                      % to original
frameDuration = round(frameDuration/2)*2; % Enforce evenness of the frame's length
shift = round(hopsize/1000*fs);
halfDur = round(frameDuration/2);
time = halfDur+1:shift:waveLen-halfDur;
N = length(time);
frameMat=zeros(frameDuration,N);
for n=1:N
    frameMat(:,n) = res(time(n)-halfDur:time(n)+halfDur-1);
end
clear res;

%% Create window matrix and apply to frames
win = blackman(frameDuration);
frameMatWin = bsxfun(@times, frameMat, win);
clear frameMat;

%% Do mean subtraction
frameMean = mean(frameMatWin,1);
frameMatWinMean = bsxfun(@minus, frameMatWin, frameMean);
clear frameMean frameMatWin frameMat;

%% Compute spectrogram matrix
specMat = zeros(floor(fs/2), size(frameMatWinMean,2));
idx = 1:floor(fs/2);
for i = 1:size(frameMatWinMean,2)
    tmp = abs( fft(frameMatWinMean(:,i),fs) );
    specMat(:,i) = tmp(idx);
end
clear frameMatWinMean tmp idx;
% specMat = abs( fft(frameMatWinMean,fs) );
% specMat = specMat(1:fs/2,:);
specDenom = sqrt( sum( specMat.^2, 1 ) );
specMat = bsxfun(@rdivide, specMat, specDenom);
clear specDenom;

%% Estimate the pitch track in 2 iterations
for Iter=1:Niter   

    [F0s,SRHVal] = SRH( specMat, nHarmonics, f0min, f0max );
    
    
    if max(SRHVal) > SRHiterThresh 
        F0medEst = median( F0s( SRHVal > SRHiterThresh ) );
        
        % Only refine F0 limits if within the original limits
        if round(0.5*F0medEst) > f0min
            f0min=round(0.5*F0medEst);
        end
        if round(2*F0medEst) < f0max
            f0max=round(2*F0medEst);    
        end
    end
    
end
clear specMat;
time=time/fs;

%% Voiced-Unvoiced decisions are derived from the value of SRH (Summation of
%% Residual Harmonics)
VUVDecisions = zeros(1,N);

if std(SRHVal) > SRHstdThresh
   VoicingThresh = VoicingThresh2;
end

VUVDecisions( SRHVal > VoicingThresh ) = 1;

end


function [F0,SRHVal] = SRH( specMat, nHarmonics, f0min, f0max )

% Function to compute Summation of Residual harmonics function
% on a spectrogram matrix, with each column corresponding to one
% spectrum.

% Initial settings
N = size(specMat, 2);
SRHmat = zeros(f0max,N);

fSeq = f0min:f0max;
fLen = length(fSeq);

% Prepare harmonic indeces matrices.
plusIdx = repmat( (1:nHarmonics)',1,fLen) .* repmat(fSeq,nHarmonics,1);
subtrIdx = round( repmat( (1:nHarmonics-1)'+.5,1,fLen) .* ...
                  repmat(fSeq,nHarmonics-1,1) );

% avoid costly repmat operation by adjusting indices
plusIdx = mod(plusIdx-1,size(specMat,1))+1;
subtrIdx = mod(subtrIdx-1,size(specMat,1))+1;

% Do harmonic summation
for n=1:N
    specMatCur = specMat(:,n);
    SRHmat(fSeq,n) = ( sum( specMatCur(plusIdx), 1 ) - ...
                       sum( specMatCur(subtrIdx), 1 ) )';
end

% Retrieve f0 and SRH value
[SRHVal,F0] = max( SRHmat );

end


