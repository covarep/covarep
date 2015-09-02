% Function to derive short term power contour as is done in Ishi et al 2008
%
% Description
% Function to derive short term power contour as is done in Ishi et al 2008
%
%
% Inputs
%  x        : [samples] [Nx1]  Speech signal
%  fs       : [Hz]      [1x1]  Sampling frequency
%
% Outputs
%  pow      : [dB]      [Mx1]  Power contour
%  pow_std  : [samples] [Px1]  Standard deviation of power contour
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] Ishi, C., Sakakibara, K-I, Ishiguro, H., (2008) `A method for 
%       automatic detection of vocal fry', IEEE TASLP, 16(1), 47-56.
%  [2] Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%  [3] Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', Computer Speech and Language 27(4), pp.
%       1028-1047.
%  [4] Drugman, T., Kane, J., Gobl, C., (2012) `Resonator-based creaky 
%       voice detection', Interspeech 2012, Portland, Oregon, USA.
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

function [pow,pow_std,pow_std_inter] = get_short_pow(x,fs)

% Get very short term power contour
veryShort_len = 4*(fs/1000); % 4ms frame length for "very short-term" analysis
veryShort_shift = 2*(fs/1000); % 2ms shift for "very short-term" analysis
veryShort_powCont = zeros(1,ceil((length(x)-veryShort_len)/veryShort_shift));
start=1;
finish=start+veryShort_len-1;

n=1;
x2 = x.^2;
while finish <= length(x)
    veryShort_powCont(n) = mean(x2(start:finish));
    start = start + veryShort_shift;
    finish=start+veryShort_len-1;
    n=n+1;
end
clear x2;

pow = 20*log10(veryShort_powCont);
inf_idx=isinf(pow);
pow(inf_idx)=min(pow(~inf_idx));

pow_std=zeros(1,length(pow));
std_len=16;

for n=std_len+1:length(pow)-std_len
    
    pow_std(n) = std(pow(n-std_len:n+std_len));
end
pow_std=medfilt1(pow_std,13);


pow_std_inter=interp1(linspace(1,length(x),length(pow)),pow_std,1:length(x));