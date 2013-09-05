% IAIF Glottal Inverse Filtering
%  [g,dg,a,ag] = iaif_ola(x,fs,winLen,winShift,p_vt,p_gl,d,hpfilt)
%
% Description
%  This function estimates vocal tract linear prediction coefficients and
%  the glottal volume velocity waveform from a speech signal frame using
%  Iterative Adaptive Inverse Filtering (IAIF) method. Analysis is carried
%  out on a fixed frame basis and waveforms are generated using overlap and
%  add.
%
% Inputs
%  x       : Speech signal frame [samples]
%  fs      : Sampling frequency [Hz]
%  winLen  : Window Length [samples] 
%  winShift: Window Shift [samples]
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
% Notes
%  This function does not perform pitch synchronous analysis. This ensures
%  the robustness of the method regardless of the GCI estimation
%  performance.
%
% Example
%  Simplest, type g = iaif(x,fs) to estimate the glottal flow of a speech
%  frame.
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

function [g,gd,a,ag] = iaif_ola(x,fs,winLen,winShift,p_vt,p_gl,d,hpfilt)

    %% Initial settings
    if nargin < 3
        winLen=25/1000*fs;
        winShift=5/1000*fs;
    end

    % Use default settings from Tuomo Raitios iaif.m implementation
    if nargin < 5
        hpfilt = 1;
        d = 0.99;
        p_gl = 2*round(fs/4000);
        p_vt = 2*round(fs/2000)+4;
    end

    % Allocate space
    gd=zeros(1,length(x));
    g=zeros(1,length(x));
    wins=zeros(1,length(x));
    a=zeros(p_vt+1,round((length(x)-winLen)/winShift));
    ag=zeros(p_gl+1,round((length(x)-winLen)/winShift));

    %% Do processing
    start=1;
    stop=start+winLen-1;
    cnt=1;

    win = hanning(winLen);
    while stop <= length(x)
    
        % Get windowed frame
        x_frame=x(start:stop);

        % Do IAIF
        [g_frame,gd_frame,a_frame,ag_frame]=iaif(x_frame,fs,p_vt,p_gl,d,hpfilt);
        a(:,cnt)=a_frame(:);
        ag(:,cnt)=ag_frame(:);
        
        % Overlap and add
        g(start:stop)=g(start:stop)+g_frame(:)'.*win.';
        gd(start:stop)=gd(start:stop)+gd_frame(:)'.*win.';

        wins(start:stop)=wins(start:stop)+win.';

        % Increment
        start=start+winShift-1;
        stop=start+winLen-1;
        cnt=cnt+1;
    end

    idx = find(wins>0);
    g(idx) = g(idx)./wins(idx);
    gd(idx) = gd(idx)./wins(idx);

return
