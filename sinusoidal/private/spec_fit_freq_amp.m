% Get the freq. and amp. parameters of a spectral peak by quandratic fitting
%
% Octave compatible
% 
% Inputs
%  S        : The spectrum containing the peak to fit
%  fs       : [Hz] The sampling frequency
%  index    : The peak position where the parabola has to be fit
%  [zp]     : the zero padding factor used for the DFT computation
%             (used for the Abe & Smith corrections)
%  [wintype]: The window type used for the DFT computation
%             wintype=2 => Hann
%
% Outputs
%  freq     : [Hz] The frequency parameter
%  amp      : The amplitude parameter (linear scale)
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

function [freq amp p] = spec_fit_freq_amp(S, fs, index, zp, wintype)

    if index<=1
        freq = index;
        amp = log(abs(S(1)));
        return
    elseif index>=length(S)
        freq = index;
        amp = log(abs(S(end)));
        return
    end

    % fit a parabola in log amplitudes
    y1 = log(abs(S(index-1)));
    y2 = log(abs(S(index)));
    y3 = log(abs(S(index+1)));
    A = (y3-y2)/2 + (y1-y2)/2;
    B = -( y1 - y2 - A );
    C = y1 - A + B;
    p = [A B C];

    % max abscissa
    di = -p(2)/(2*p(1));

    % max amplitude
    logamp = p(1)*di*di + p(2)*di + p(3);

    if nargin>3 && ~isempty(zp)
        % Use Abe & Smith corrections
        if wintype==2;
            c0=0.247560; c1=0.084372; c2=-0.090608; c3=-0.055781;   % Hann
        else
            error('Unknown correction terms for the given window for Abe&Smith corrections');
        end
        delta = di;

        % Frequency correction
        ksi = c0*zp^(-2) + c1*zp^(-4);                  % (3)
        di = delta + ksi*(delta-0.5)*(delta+0.5)*delta; % (1)

        % Amplitude correction
        ita = c2*zp^(-4) + c3*zp^(-6);                  % (4)
        logamp = logamp + ita*di^2;                     % (2)
    end

    freq = index + di;
    freq = (freq-1)*(fs/length(S));

    amp = exp(logamp);
return

