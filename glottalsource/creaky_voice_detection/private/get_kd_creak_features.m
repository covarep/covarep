% Function to extract Kane-Drugman (KD) features used for creaky voice
% detection.
%
% Description
%
%  This version is a slightly later one than that described in the
%  above published 2013 CSL paper [2]. The algorithm here rather than using binary
%  decision trees using artificial neural networks and combines the
%  features used in the CSL paper with those proposed in Ishi et al.
%  (2008). This updated version has been submitted to CSL for a special
%  issue on glottal source processing on April 14th 2013. It will have
%  reference [1].
%
% Inputs
%  x        : [samples] [Nx1] Speech signal
%  fs       : [Hz]      [1x1] sampling frequency
%
% Outputs
%  H2H1     : [dB]      [Mx1] H2-H1 parameter
%  res_p    : [samples] [Mx1] Residual peak prominence
%  ZCR      : [samples] [Mx1] Zero-crossing rate
%  F0       : [Hz]      [Mx1] Fundamental frequency (from SRH algorithm)
%  F0mean   : [Hz]      [Mx1] Mean Fundamental frequency (from SRH algorithm)
%  enerN    : [dB]      [Mx1] Normalised energy contour
%  powStd   : [samples] [Mx1] Standard deviation of power contour frames
%  creakF0  : [Hz]      [Mx1] F0 as outputted by the H2-H1 algorithm
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%  [2] Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', Computer Speech and Language 27(4), pp.
%       1028-1047.
%  [3] Drugman, T., Kane, J., Gobl, C., (2012) `Resonator-based creaky 
%       voice detection', Interspeech 2012, Portland, Oregon, USA.
%
% Octave Compatibility
%  As this algorithm relies on the Matlab neural networks toolbox it is not
%  compatible with Octave. Furthermore, this algorithm requires Matlab
%  R2011 or later
%
% Copyright (c) 2013 University of Mons, FNRS & 2013 Trinity College Dublin
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
% Authors
%  Thomas Drugman <thomas.drugman@umons.ac.be> & John Kane <kanejo@tcd.ie>
%

function [H2H1,res_p,ZCR,F0,F0mean,enerN,pow_std,creakF0] = get_kd_creak_features(x,fs)

%% Initial settings
F0min=20;
F0max=500;
winShift=round(10/1000*fs); % Sample every 10 ms
res = lpcresidual(x,round(25/1000*fs),round(5/1000*fs),round(fs/1000)+2);
[f0,VUV] = pitch_srh(x,fs,F0min,F0max,10);
F0mean=median(f0(VUV==1&f0>F0min&f0<F0max));

%% Extract features for creak detection
time=winShift:winShift:length(x);

[~,Es,ZC_ms,Xpos] = sil_unv_features(x,fs,32);
ZC_inter = interp1(Xpos,ZC_ms,1:length(x));
ener_norm = Es-max(Es);
Es_inter=interp1(linspace(1,length(x),length(ener_norm)),ener_norm,1:length(x));
[~,~,pow_std_inter] = get_short_pow(x,fs);


[H2H1,F0_creak] = get_creak_h2h1(res,fs,F0mean);
peak_inter = res_peak(x,fs,F0mean,res,Es);

%% Do resampling
peak_inter=peak_inter(:,1);
Es_re=Es_inter(time);
pow_re=pow_std_inter(time);
H2H1_re=H2H1(time);
peak_re=peak_inter(time);
peak_re(isnan(peak_re))=0;
ZC_re=ZC_inter(time);
ZC_re(isnan(ZC_re))=0;
F0_creak=F0_creak(time);


%% Save to structs
ZCR=ZC_re(:);
enerN=Es_re(:);
creakF0=F0_creak(:);
pow_std=pow_re(:);
H2H1=H2H1_re(:);
res_p=peak_re(:);
F0mean=zeros(size(res_p))+F0mean;
F0=interp1(linspace(1,length(F0mean),length(f0)),f0,1:length(F0mean));




