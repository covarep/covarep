% Create a spectrum from ARMA coefficients (g,b,a)
%
% Inputs
%  g          : gain [linear energy]
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

function S = gba2hspec(g, b, a, dftlen)

    if prod(size(g))>length(g) || prod(size(b))>length(b) || prod(size(a))>length(a)
        Nframes = max([size(g,1), size(b,1), size(a,1)]);
        if mod(dftlen,2)==0; S = zeros(Nframes, dftlen/2+1);
        else;                S = zeros(Nframes, (dftlen-1)/2+1);   end
        for k=1:Nframes
            if length(b)>0
                if length(a)>0
                    S(k,:) = gba2hspec(g(k), b(k,:), a(k,:), dftlen);
                else
                    S(k,:) = gba2hspec(g(k), b(k,:), 1, dftlen);
                end
            else
                if length(a)>0
                    S(k,:) = gba2hspec(g(k), 1, a(k,:), dftlen);
                else
                    S(k,:) = g(k);
                end
            end
        end

    else
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
    end

return
