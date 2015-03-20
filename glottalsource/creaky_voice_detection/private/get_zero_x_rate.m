% Function to derive zero-crossing rate parameter
%
% Description
% Function to derive zero-crossing rate parameter
%
%
% Inputs
%  x        : [samples] [Nx1]  Speech signal
%  winLen   : [samples] [1x1]  Window length
%  shift    : [samples] [1x1]  Window shift
%
% Outputs
%  zeroXingRate: [samples] [Mx1] Zero-crossing rate parameter
%  pos:          [samples] [Mx1] Time locations of parameter
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

function [zeroXingRate,pos] = get_zero_x_rate(x,winLen,shift)

% Function to measure zero crossing rate to be used as a measure in creaky
% voice detection

%% Calculate Intraframe Periodicity Measure (IFP)
start=1;
finish=start+winLen-1;
n=1;

while finish <= length(x)
    frame = x(round(start):round(finish)); % select frame
    
    xing_idx=zerocros(frame,'b');
    zeroXingRate(n) = length(xing_idx)/length(frame);
    
    start=start+shift;
    finish=start+winLen-1;
    pos(n) = mean([start finish]);
    n=n+1;
end

%%%%%%%
% INTERNAL FUNCTIONS
%%%%%%%

 function [t,s]=zerocros(x,m)
%ZEROCROS finds the zeros crossings in a signal [T,S]=(X,M)% find zero
%crossings in a signal

if nargin<2
   m='b';
end
s=x>=0;
k=s(2:end)-s(1:end-1);
if any(m=='p')
   f=find(k>0);
elseif any(m=='n')
   f=find(k<0);
else
   f=find(k~=0);
end
s=x(f+1)-x(f);
t=f-x(f)./s;
if ~nargout
   n=length(x);
   plot(1:n,x,'-',t,zeros(length(t),1),'o');
end

