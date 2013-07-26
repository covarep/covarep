% Drop the negative frequencies of a DFT
%
% Description
%  A half spectrum contains only positive frequencies, DC and Nyquist frequencies.
%  The output is supposed to be of even length
%
% Inputs
%  S          : A spectrum as made by the fft function 
%               S(1:end/2+1) = hS
%               (even length)
%
% Outputs
%  hS         : A half spectrum: positive frequencies with DC an Nyquist frequency
%               indices S(1:end/2+1)
%
% Example
%  X = spec2hspec(fft(x, fftlen));
%  x = hspec2spec(X);
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

function hS = spec2hspec(S)

    S = S(:).';

    if mod(length(S),2)==0; hS=S(1:end/2+1);
    else;                   fhS=S(1:(end-1)/2+1);   end
    
return
