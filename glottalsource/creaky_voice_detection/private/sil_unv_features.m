% Function to derive audio features for silence and unvoiced regions
%
% Description
% Function to derive audio features for silence and unvoiced regions
%
%
% Inputs
%  x        : [samples] [Nx1]  Speech signal
%  fs       : [Hz]      [1x1]  Sampling frequency
%  frame_len_ms: [ms] [1x1]  Window length
%
% Outputs
%  Ts       : [samples] [Px1]  Locations of parameter values
%  Es       : [dB]      [Mx1]  Energy contour
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

function [Ts,Es,ZCs_ms,Xpos,pos] = sil_unv_features(x,fs,frame_len_ms)

if nargin < 3
    frame_len_ms=10;
end
frame_len=frame_len_ms/1000*fs;
frame_shift_ms=5;
Shift=frame_shift_ms/1000*fs;
 
Start=1; % intialise
Stop=Start+frame_len-1;
pos=mean([Start Stop]);

[ZCs,Xpos] = get_zero_x_rate(x,frame_len,Shift);
ZCs_ms=ZCs/1000*fs;

Ind=1;

Es=zeros(1,floor((length(x)-Stop+1)/Shift)+1);
Ts=zeros(size(Es));

win = hanning(Stop-Start+1);
 
while Stop<length(x)
   
    Mid=round((Start+Stop)/2);
    Ts(Ind)=Mid;
   
    Sig=x(Start:Stop); % Select frame segment
    Sig=Sig(:).*win; % Window segment

   
    Es(Ind)=mean(Sig.^2); % Get energy value
   
    Start=Start+Shift;
    Stop=Stop+Shift;
    pos(Ind)=mean([Start Stop]);
    Ind=Ind+1;
end 

Es=log(Es);