% Gives the spectrum of a derivative effect (i.e. a zero at zero frequency)
%
% Input
%  dftlen : the number of bin in spectrum (full DFT length).
%  fs     : [Hz] The sampling frequency
%
% Output
%  S      : The spectrum of the derivative effect
%
% Copyright (c) 2011 University of Crete - Computer Science Department
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

function S = spec_derivative(dftlen, fs)

    % The half-spectrum of a zero at frequency zero is:
    hS = fs*(2i*pi/dftlen)*(0:dftlen/2).';
    
    if mod(dftlen,2)==1
        % If DFT length is odd
        S = [hS; hS(end:-1:2)];
    else
        % If DFT length is even
        hS(end) = 0;
        S = hspec2spec(hS);
    end 

return
