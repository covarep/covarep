% The Short Time Fan-Chirp Transform (spectrogram-like FChT)
% 
% This function is Octave compatible, but dreadfully slow, because of
%  the complex exponential function in fcht.m.
%
% Description
%  Compute the FChT (see fcht.m) along frames
%  Have a look at fcht.m
%
% Inputs
%  x        : A full signal
%  win      : The window to use
%  noverlap : [samples] Number of overlapping samples from one frame to the next
%  dftlen   : [samples] The number of bins in the discrete spectral representation
%  fs       : [Hz] The sampling frequency
%  f0s      : [Nx2] [s,Hz] A time-data vector specifying an estimated fundamental
%             frequency at given times.
%
% Outputs
%  X       : The FChT transform
%  F       : The frequency values corresponding to each bin
%  T       : Temporal position of the the window centers
%  as      : Estimated relative f0 slope at window centers
%
% Example
%  Please see the HOWTO_spectra
%
% References
%  [1] Kepesi, M., Weruaga, L.: Adaptive Chirp-based time-frequency analysis of
%      speech signals, Speech communication 48(5), 474-492, 2006.
%  [2] Weruaga, L, Kepesi, M: The Fan-Chirp Transform for non-stationary harmonic
%      signals, Signal Processing 87(6), 1504-1522, 2007.
%
% Copyright (c) 2012 University of Crete - Computer Science Department
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
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function [X F T as] = fchtgram(x, win, noverlap, dftlen, fs, f0s)

    winlen = length(win);
    F = fs*(0:dftlen/2)/dftlen;

    % Seems to correspond to Matlab times od spectrogram function
    % TODO faster
    T = [];
    ids = 1;
    ind = 1;
    while ids(end)<=length(x)
        ids = (ind-1)*(winlen-noverlap) + (1:winlen);
        if ids(end)<=length(x)
            T = [T (ids((winlen-1)/2+1)-1)/fs];
        end
        ind = ind+1;
    end

    f0 = interp1(f0s(:,1), f0s(:,2), T, 'linear','extrap');

%      hold off;
%      plot((0:length(x)-1)/fs, x, 'k');
%      hold on;
%      plot(f0s(:,1), hz2oct(f0s(:,2)), 'b');
%      plot(T, hz2oct(f0), 'k');

%      hold on;
%      plot(T, f0, 'b');

    as = ones(length(T),1);
    X = zeros(dftlen/2+1,length(T));
    pb = progressbar(length(T));
    for ind=1:length(T)
        ids = (ind-1)*(winlen-noverlap) + (1:winlen);

        if ind==1 || ind==length(T)
            a = 0;
        else
            % 1st order centered derivative approximation
            a = (1/f0(ind))*(f0(ind+1)-f0(ind-1))/(T(ind+1)-T(ind-1));
        end
        a=a/fs;

        as(ind) = a;

        Xs = fcht(x(ids).*win, a, dftlen);

        X(:,ind) = Xs(1:end/2+1);

        pb = progressbar(pb, ind);
    end

return
