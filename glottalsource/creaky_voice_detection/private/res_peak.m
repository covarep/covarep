% Function to derive residual peak prominence parameter
%
% Description
% Function to derive residual peak prominence parameter
%
%
% Inputs
%  x        : [samples] [Nx1]  Speech signal
%  fs       : [Hz]      [1x1]  Sampling frequency
%  res      : [samples] [Nx1]  Linear prediction residual of speech
%  Es       : [dB]      [Mx1]  Energy contour
%
% Outputs
%  peak_inter:[samples] [Nx2]  Interpolated peak prominence contour
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

function [peak_inter,peak_prom,peak_t,rep] = res_peak(x,fs,F0mean,res,Es)

% Function to generate residual peak prominence contour which makes up the
% second component of the creak detection algorithm.


% Different settings according to the speakers baseline pitch
if F0mean>100
    maxCreakF0=90;
elseif F0mean < 100 && F0mean>=85
    maxCreakF0=65;
elseif F0mean<85
    maxCreakF0=55;
else maxCreakF0=80;
end

% Set window length based on maximum possible creaky F0
winLen=round(fs/maxCreakF0)*2; 

% Resonator settings
Phi=2*pi*1*F0mean/fs;
Rho=0.8;
rep=filter([1 0 0],[1 -2*Rho*cos(Phi) Rho^2],res);

% Measure residual peak prominence
[peak_prom,peak_t] = get_res_peak_prom(rep,fs,winLen,x,Es);

% Interpolate
if length(peak_prom)>1
    peak_inter=interp1(peak_t,peak_prom,1:length(x));
else peak_inter=zeros(1,length(x));
end

peak_inter=[peak_inter(:) (1:length(x))'];