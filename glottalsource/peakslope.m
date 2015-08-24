% A parameter which is essentially a correlate of spectral tilt, derived
% following wavelet analysis. This parameter is effective at discriminating
% lax-tense phonation types.
%
% Octave compatible
%
% Description
%  The peakslope parameter is derived by using wavelet analysis to decompose
%  the speech signal (or glottal source signal) into octave bands. A sliding
%  window is then used for measuring time-domain maxima in the different
%  frequency bands. A regression line is fit to log10 of these peaks, and
%  the slope coefficient is used as the peakSlope parameter. For the wavelet
%  analysis, different scaled versions of a cosine modulated gaussian pulse
%  are convolved with the input signal (as is done in d'Alessandro 2011 and
%  Tuan & d'Alessandro 1999).
%
% Inputs
%  s               : [samples] [Nx1] input signal
%                    (speech signal or glottal source)
%  fs              : [Hz]      [1x1] sampling frequency
%
% Outputs
%  PS              : [time,slope] [Mx2] outputted peakSlope parameter,
%                    and time instance [s]
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] Kane, J., Gobl, C., (2011) ``Identifying regions of non-modal 
%      phonation using features of the wavelet transform'', Proceedings of
%      Interspeech, Florence, Italy, pp. 177-180, 2011.
%
% Copyright (c) 2011 Trinity College Dublin
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
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  John Kane kanejo@tcd.ie

function PS = peakslope(s,fs)

%% Initial settings
frameLen_ms = 40; % Frame length chosen to ensure one pulse length down to f0=25 Hz
frameShift_ms = 10; % Frame shift set to 10 ms
frameLen = (frameLen_ms/1000)*fs; % Convert frame length to samples
frameShift = (frameShift_ms/1000)*fs; % Convert frame shift to samples

PS=zeros(round((length(s)-frameLen)/frameShift),2);

%% Do wavelet decomposition
i=0:6; % i=0:6 => 8 kHz, 4 kHz, 2 kHz, 1 kHz, 500 Hz, 250 Hz, 125 Hz

y=zeros(length(i),length(s)); % Allocate space for the different frequency bands

for n=1:length(i) 
    h_i = daless_MW(i,n,fs); % Generate mother wavelet
    y(n,:) = do_daless_filt(s,h_i); % Carry out zero-phase filtering
end

%% Measure peakSlope per frame
start=1;
finish = start+frameLen-1;
m=1;

while finish <= length(s)
    maxima = max(abs(y(:,start:finish)),[],2)';
    maxima = log10(maxima(end:-1:1)); % reverse order to follow frequency order and convert to dB
    t=1:length(maxima);
    p=polyfit(t,maxima,1); % do straight line regression fitting
    PS(m,2) = p(1); % take slope coefficient from regression line
    PS(m,1) = (round((start+finish)/2)-1)/fs;
    
    m=m+1;
    start = start+frameShift;
    finish = start+frameLen-1;
end
PS(isnan(PS(:,2)),2)=0;

end

% Changes
% 2012-11-01> John Kane kanejo@tcd.ie
% A modification to description in the 2011 Interspeech paper (reference
% [1] above) was made. Now before fitting the regression line maxima are
% converted to log10.
%
