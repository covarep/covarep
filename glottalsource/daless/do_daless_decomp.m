% Function to do full wavelet based decomposition of an input signal
%
% Octave compatible
%
% Description
%  This function decomposes an input signal into different octaves bands
%  using the mother wavelet function generated in daless_MW.m
%
% Inputs
%  s               : [samples] [Nx1] input signal (speech signal or glottal source)
%  i               : [samples] [Mx1] scales to be used, e.g., i=0:6;
%  fs              : [Hz]      [1x1] sampling frequency
%
% Outputs
%  y               : [samples] [MxN] Decomposed input signal
%  y_norm          : [samples] [MxN] Decomposed input signal, with each row normalised in amplitude. 
%
% Example
%  [y,y_norm] = do_daless_decomp(s,fs,i)
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

function [y,y_norm] = do_daless_decomp(s,fs,i)

if nargin < 3
    i=0:6; % i=0:6 => 8 kHz, 4 kHz, 2 kHz, 1 kHz, 500 Hz, 250 Hz, 125 Hz
end

%% Allocate space
y=zeros(length(i),length(s));
y_norm=zeros(length(i),length(s));

%% Do filtering
for n=1:length(i)
    
    h_i = daless_MW(i,n,fs);
    y_i = do_daless_filt(s,h_i);
    y(n,:) = y_i;
    y_norm(n,:) = y_i/max(abs(y_i));
end
