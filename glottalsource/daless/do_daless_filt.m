% Function to convolve a given input signal, s, with the wavelet function h_i
%
% Octave compatible
%
% Description
%  Function to convolve a given input signal, s, with the wavelet function h_i
%
% Inputs
%  s               : [samples] [Nx1] input signal (speech signal or glottal source)
%  h_i             : [samples] [Mx1] mother wavelet function
%
% Outputs
%  y_i             : [samples] [Nx1] Decomposed input signal 
%
% Example
%  y_i = do_daless_filt(s,h_i);
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

function y_i = do_daless_filt(s,h_i)

s_h_i=conv(s,h_i);
halfLen = ceil(length(h_i)/2);
y_i = s_h_i(halfLen:halfLen+length(s)-1);
