% Function to generate mother wavelet function used in peakSlope and MDQ
% functions.
%
% Octave compatible
%
% Description
%  This function generates a cosine modulated gaussian pulse, which is used
%  as a mother wavelet in both peakSlope.m and MDQ.m functions. This is the
%  same mother wavelet as is used in (as is done in d'Alessandro 2011 and
%  Tuan & d'Alessandro 1999).
%
% Inputs
%  i               : [samples] [Nx1] scales to be used, e.g., i=0:6;
%  fs              : [Hz]      [1x1] sampling frequency
%  n               : [sample]  [1x1] index of current scale
%
% Outputs
%  h_i             : [samples] [Mx1] mother wavelet function
%
% Example
%  h_i = daless_MW(i,n,fs);
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
%

function h_i = daless_MW(i,n,fs)

s=2.^i; % convert to octave bands
f_o = fs/2;
tau = 1/(2*f_o);

t=(-1000:1000)./fs;
mother_wavelet = @(s_n) -cos(2*pi*f_o.*(t./s_n)).*exp(-((t./s_n).^2)/(2*(tau^2))); % mother wavelet
h_i = mother_wavelet(s(n));
