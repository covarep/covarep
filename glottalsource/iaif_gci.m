% IAIF Glottal Inverse Filtering
%
% Description
%  This function estimates vocal tract linear prediction coefficients and
%  the glottal volume velocity waveform from a speech signal frame using
%  Iterative Adaptive Inverse Filtering (IAIF) method. Analysis is carried
%  out on a GCI-synchronous basis and waveforms are generated using overlap 
%  and add.
%
% Inputs
%  x       : Speech signal frame [samples]
%  fs      : Sampling frequency [Hz]
%  GCI     : Glottal closure instants [seconds] 
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
%
% Example
%  Simplest, type g = iaif_gci(x,fs,GCI) to estimate the glottal flow of a speech
%  signal.
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
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Octave compatible
% 
% Author
%  Tuomo Raitio <tuomo.raitio@aalto.fi>
%  Overlap and add part by John Kane <kanejo@tcd.ie>

function [g,gd,a,ag] = iaif_gci(x,fs,GCI,p_vt,p_gl,d,hpfilt)

%% Initial settings
% Use default settings from Tuomo Raitios iaif.m implementation
if nargin < 4
    hpfilt = 1;
    d = 0.99;
    p_gl = 2*round(fs/4000);
    p_vt = 2*round(fs/2000)+4;
end
GCI=round(GCI*fs);

% Allocate space
N=length(GCI);
gd=zeros(1,length(x));
g=zeros(1,length(x));
a=zeros(p_vt+1,N);
ag=zeros(p_gl+1,N);
hpfilter_out = [];

%% Do processing
for n=1:N
   
    % Get windowed frame
    if n==1
        T0=GCI(n+1)-GCI(n);
    else T0=GCI(n)-GCI(n-1);
    end
    start=GCI(n)-T0;
    stop=GCI(n)+T0;
    
    if start < 1
        start=1;
    end
    if stop > length(x)
        stop=length(x);
    end
    
    x_frame=x(start:stop);
    x_win=x_frame(:).*hanning(length(x_frame));
    
    % Do IAIF
    [g_frame,gd_frame,a_frame,ag_frame,hpfilter_out]=iaif(x_win,fs,p_vt,p_gl,d,hpfilt,hpfilter_out);
    if isempty(g_frame)==0
        a(:,n)=a_frame(:);
        ag(:,n)=ag_frame(:);

        % Overlap and add
        g(start:stop)=g(start:stop)+g_frame(:)';
        gd(start:stop)=gd(start:stop)+gd_frame(:)';
    end
    
end
