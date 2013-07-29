% RCVD_reson_GCI.m
% Function to use the resonator used in RCVD (creaky voice detection), 
% applied to the LP-residual signal and give output
%
% Octave compatible
%
% Description
% Function to use the resonator used in RCVD (creaky voice detection), 
% applied to the LP-residual signal and give output
%
% Inputs
%  res             : [samples] [Nx1] Linear Prediction residual
%  fs              : [Hz]      [1x1] sampling frequency
%  F0mean          : [Hz]      [1x1] Mean fundamental frequency
%
% Outputs
%  y               : [samples] [Nx1] Resonator output
%
% Example
%  GCI = SE_VQ(x,fs);
%
% References
%  [1] Kane, J., Gobl, C., (2013) `Evaluation of glottal closure instant 
%       detection in a range of voice qualities', Speech Communication 
%       55(2), pp. 295-314.
%
% Copyright (c) 2013 Trinity College Dublin
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
%
% $Id <info set by the versioning system> $

function y = RCVD_reson_GCI(res,fs,F0mean)

Phi=2*pi*1*F0mean/fs;
Rho=0.9; % Set to a narrow bandwidth
rep=filtfilt([1 0 0],[1 -2*Rho*cos(Phi) Rho^2],res); % Filter forwards and backwards to have zero-phase
y=rep/max(abs(rep));
