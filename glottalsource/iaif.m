% IAIF Glottal Inverse Filtering
%  [g,dg,a,ag] = iaif(x,fs,p_vt,p_gl,d,hpfilt)
%
% Description
%  This function estimates vocal tract linear prediction coefficients and
%  the glottal volume velocity waveform from a speech signal frame using
%  Iterative Adaptive Inverse Filtering (IAIF) method.
%
% Inputs
%  x       : Speech signal frame [samples]
%  fs      : Sampling frequency [Hz]
%  p_vt    : Order of LPC analysis for vocal tract
%  p_gl    : Order of LPC analysis for glottal source
%  d       : Leaky integration coefficient (e.g. 0.99)
%  hpfilt  : High-pass filter flag (0: do not apply, 1...N: apply N times)
%
% Outputs
%  g       : Glottal volume velocity waveform
%  dg      : Glottal volume velocity derivative waveform
%  a       : LPC coefficients of vocal tract
%  ag      : LPC coefficients of source spectrum
%
% Notes
%  This function does not perform pitch synchronous analysis. This ensures
%  the robustness of the method regardless of the GCI estimation
%  performance.
%
% Example
%  Simplest, type g = iaif(x,fs) to estimate the glottal flow of a speech
%  frame.
%  And see the HOWTO_glottalsource.m example file.
%
% References
%  [1] P. Alku, "Glottal wave analysis with pitch synchronous iterative
%      adaptive inverse filtering", Speech Communication, vol. 11, no. 2-3,
%      pp. 109–118, 1992.
%
% Copyright (c) 2013 Aalto University
%
% License
%  Please contact the author or the Copyright owner for any IP related question.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% TODO Octave compatibility?
% 
% Author
%  Tuomo Raitio <tuomo.raitio@aalto.fi>

function [g,dg,a,ag] = iaif(x,fs,p_vt,p_gl,d,hpfilt)

% Set default parameters
if nargin < 6
    hpfilt = 1;
    if nargin < 5
        d = 0.99;
        if nargin < 4
            p_gl = 2*round(fs/4000);
            if nargin < 3
                p_vt = 2*round(fs/2000)+4;
                if nargin < 2
                    disp('Error: Not enough input arguments.');
                end
            end
        end
    end
end
preflt = p_vt+1;

% Ensure column vector form
x = x(:);

% High-pass filter speech in order to remove possible low frequency
% fluctuations (Linear-phase FIR, Fc = 70 Hz)
if hpfilt > 0
    Fstop = 40;                 % Stopband Frequency
    Fpass = 70;                 % Passband Frequency
    Nfir = round(300/16000*fs); % FIR numerator order
    if mod(Nfir,2) == 1
        Nfir = Nfir + 1;
    end
    B = hpfilter_fir(Fstop,Fpass,fs,Nfir);
    for i = 1:hpfilt
        x = filter(B,1,[x ; zeros(round(length(B)/2)-1,1)]);
        x = x(round(length(B)/2):end);
    end
end

% Estimate the combined effect of the glottal flow and the lip radiation
% (Hg1) and cancel it out through inverse filtering. Note that before
% filtering, a mean-normalized pre-frame ramp is appended in order to
% diminish ripple in the beginning of the frame. The ramp is removed after
% filtering.
Hg1 = lpc(x.*hanning(length(x)),1);
y = filter(Hg1,1,[linspace(-x(1),x(1),preflt)' ; x]);
y = y(preflt+1:end);

% Estimate the effect of the vocal tract (Hvt1) and cancel it out through
% inverse filtering. The effect of the lip radiation is canceled through
% intergration. Signal g1 is the first estimate of the glottal flow.
Hvt1 = lpc(y.*hanning(length(x)),p_vt);
g1 = filter(Hvt1,1,[linspace(-x(1),x(1),preflt)' ; x]);
g1 = filter(1,[1 -d],g1);
g1 = g1(preflt+1:end);

% Re-estimate the effect of the glottal flow (Hg2). Cancel the contribution
% of the glottis and the lip radiation through inverse filtering and
% integration, respectively.
Hg2 = lpc(g1.*hanning(length(x)),p_gl);
y = filter(Hg2,1,[linspace(-x(1),x(1),preflt)' ; x]);
y = filter(1,[1 -d],y);
y = y(preflt+1:end);

% Estimate the model for the vocal tract (Hvt2) and cancel it out through
% inverse filtering. The final estimate of the glottal flow is obtained
% through canceling the effect of the lip radiation.
Hvt2 = lpc(y.*hanning(length(x)),p_vt);
dg = filter(Hvt2,1,[linspace(-x(1),x(1),preflt)' ; x]);
g = filter(1,[1 -d],dg);
g = g(preflt+1:end);
dg = dg(preflt+1:end);

% Set vocal tract model to 'a' and glottal source spectral model to 'ag'
a = Hvt2;
ag = Hg2;

return


function B = hpfilter_fir(Fstop,Fpass,fs,N)
% FIR least-squares Highpass filter design using the FIRLS function
%
% Tuomo Raitio
% 10.7.2012

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fstop Fpass fs/2]/(fs/2), [0 0 1 1], [1 1]);
Hd = dfilt.dffir(b);

% Save to B
B = Hd.Numerator;

return
