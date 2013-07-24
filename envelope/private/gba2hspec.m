% Create a spectrum from ARMA coefficients (g,b,a)
%
% USAGE
%  S = Fgba2spec(g, b, a, dftlen)
%
% INPUT
%  g          : gain
%  b          : MA coefficients
%  a          : AR coefficients
%  dftlen     : length of the spectrum
%
% OUTPUT
%  S          : a spectrum made from transfer function coefficients
%
% EXAMPLE
%
% REFERENCE
%
% SEE AlSO
%  Fcc2spec
%
% AUTHOR
%  degottex@ircam.fr
%
% COPYRIGHT
%  Copyright (c) 2008 IRCAM/CNRS - Degottex
%
% $Id: gba2spec.m 214 2012-12-14 12:22:05Z degottex $

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
