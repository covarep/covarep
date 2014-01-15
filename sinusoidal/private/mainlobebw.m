% Estimate numerically the bandwidth of a given window
%
% Octave compatible
% 
% Description
%  The half bandwidth is defined as the crossing of the frequency amplitude response
%  by the x axis, thus on a linear scale and using the real part only.
%
% Input
%  win    : The sampled window
%  fs     : Sampling frequency
%
% Output
%  bw     : Estimated window bandwidth
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

function bw = mainlobebw(win, fs)

    winlen = length(win);
    dftlen = 2^(nextpow2(winlen)+2);
%      dftlen = 4*winlen;

    if 0
        W = lin2db(fft(win,dftlen));
        W = W(1:dftlen/2)-W(1);
        wids = max(find(W>-12));
        Bi = interp1(W(wids:wids+1), wids:wids+1, -12, 'linear', 'extrap');
        bw = 2*fs*(Bi-1)/dftlen;
    else
        W = fft(win, dftlen).';
        D = delay2spec((winlen-1)/2, dftlen);
        W = W.*D;
        B1 = min(find(real(W) ~= abs(real(W))))-2;
        B2 = v_findpeaks(-abs(W(1:end/2)));
        B2 = B2(1)-1;
        B = min(B1,B2);

%          hold off;
%          plot((0:length(W)-1), abs(real(W)), 'b');
%          hold on;

        opts.TolX = 0.01/fs;

%          ks = (0:0.01:dftlen/2)';
%          for n=1:length(ks)
%              Bs(n) = fourierconv(ks(n), win,dftlen);
%          end
%  %          plot(ks, Bs, 'k');
%          hold on;
%  keyboard

%          [B,fval,exitflag] = fzero(@(k)fourierconv(k,win,dftlen),B-2,opts);

        [B,fval,exitflag] = fminbnd(@(k)fourierconv(k,win,dftlen), B-1,B+1,opts);

        if exitflag<0; error('Cannot find the bandwidth of this window'); end

        bw = 2*fs*B/dftlen;
    end

return

function e = fourierconv(k, win, N)

    n = (-(length(win)-1)/2:(length(win)-1)/2)';

    e = sum(win.*exp((-1i*2*pi*k/N).*n));

    e = abs(real(e));

%      plot(k, e, 'xr');

return
