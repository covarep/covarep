% Minimal interface-function between Octave's specgram and Matlab's spectrogram
%
% This function has the only purpose to make the HOWTO_spectra working
% under Octave. Therefore, currently, only the following Matlab's usage of
% specrogram.m is supported:
%  [S,F,T] = spectrogram(x,window,noverlap,F,fs);
% 
% Copyright (c) 2013 University of Crete - Computer Science Department
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

function [S, F, T] = spectrogram(x, window, noverlap, dftlen, fs)

    [S, F, T] = specgram(x, dftlen, fs, window, noverlap);

return
