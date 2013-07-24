% Complete a half spectrum to a full spectrum creating a real signal
%
% A half spectrum contains only positive frequencies, DC and Nyquist frequencies
% The output is supposed to be of even length
%
% Input
%  hS         : a half spectrum: positive frequencies with DC an Nyquist frequency
%               indices S(1:end/2+1)
%
% Output
%  S          : a spectrum as made by the fft function 
%               S(1:end/2+1) = hS
%               (even length)
%
% Example
%  X = fft(x, fftlen);
%  hX = X(1:end/2+1);
%  X-hspec2spec(hX)
%
% Copyright (c) 2008 Ircam/CNRS-UMR9912-STMS
%
% License
%  This file is part of libphoni.
%  libphoni is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  libphoni is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%  You should have received a copy of the GNU Lesser General Public License
%  along with libphoni.  If not, see <http://www.gnu.org/licenses/>.
%  
% Author
%  Gilles Degottex <gilles.degottex@ircam.fr>
%
% $Id$

function S = hspec2spec(hS, has_nyquist)

    hS = hS(:).';

    % full-spectrum has even length 
    if nargin<2 || has_nyquist
        S = [hS, conj(hS(end-1:-1:2))];

    % full-spectrum has odd length
    else
        S = [hS, conj(hS(end:-1:2))];
    end
    
return
