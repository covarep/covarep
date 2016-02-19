% Complete a half spectrum to a full spectrum corresponding to a real signal
%
% Description
%  A half spectrum contains only positive frequencies, DC and Nyquist frequencies.
%  The output is supposed to be of even length
%
% Inputs
%  hS         : a half spectrum: positive frequencies with DC an Nyquist frequency
%               indices S(1:end/2+1). If hS is a vector, hS should have the
%               size [1 length(S(1:end/2+1))]. hS can be matrix of size [n length(S(1:end/2+1))].
%
% Outputs
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
%  Gilles Degottex <gilles.degottex@ircam.fr>
%

function S = hspec2spec(hS, has_nyquist)

    if size(hS,2)==1, hS = hS'; end

    % full-spectrum has even length 
    if nargin<2 || has_nyquist
        S = [hS, conj(hS(:,end-1:-1:2))];

    % full-spectrum has odd length
    else
        S = [hS, conj(hS(:,end:-1:2))];
    end

end
