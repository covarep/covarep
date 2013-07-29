% Create a spectrum from ARMA coefficients (g,b,a)
%
% Inputs
%  g          : gain
%  b          : MA coefficients
%  a          : AR coefficients
%  dftlen     : length of the spectrum
%
% Outputs
%  S          : a spectrum made from transfer function coefficients
%
% Example
%
% Copyright (c) 2007 Ircam-CNRS UMR9912-STMS
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
%  Gilles Degottex <degottex@ircam.fr>
%

function S = gba2spec(g, b, a, dftlen)

    a = a(:)';
    b = b(:)';

    if length(b)==1
        if length(a)==1
            S = g.*(b/a) * ones(dftlen,1);
        else
            S = g.*b ./fft(a, dftlen);
        end
    else
        if length(a)==1
            S = g.*fft(b, dftlen)./a;
        else
            S = g.*fft(b, dftlen)./fft(a, dftlen);
        end
    end

    if mod(dftlen,2)==0; S=S(1:end/2+1);
    else;                S=S(1:(end-1)/2+1);   end

return
