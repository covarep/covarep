% Estimate of partials frequency, amplitude and phase given a gross fundamental frequency
%
% Octave compatible
% 
% Description
%  Estimate of partials frequency, amplitude and phase
%  All returned partials frequencies are relative to spectrum size !
%
% INPUT
%  S            : spectrum (FFT representation)
%
%  fs           : sampling frequency
%
%  f0           : undamental frquency in Hertz
%
%  [method]     : 0:'nearest': no refinement, simple spectrum sampling
%                 1:'linear':  refine amplitude linearly
%                              phase linearly
%                 2:'refine':  refine frequency by quadratic interpolation
%                   (default)  amplitude by quadratic interpolation
%                              phase by linear interpolation
%
%  [sampling]   : first estimate of partial sampling
%                 sampling(1) is suposed to be on the first harmonic
%                 eg. sampling=w0*(1:floor(length(S)/2-1)/w0)+1;
%
%  [zp]         : used zero-padding
%
% OUPUT
%  partials     : vector of partials [freq;amp;phase;harm;amp_var]
%
%  f0           : refined fundamental frequency
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

function sins = spec_getsins_f0(S, fs, f0, max_h, varargin)

    if nargin<4; max_h=[]; end

    S = S(:).';

    dftlen = length(S);

    step = dftlen*f0/fs;

    if ~isempty(max_h); max_h = min(max_h,floor((fs/2-f0/2)/f0));
    else                max_h = round((fs/2-f0/2)/f0); end

    sins = zeros(4,max_h);
    for h=1:max_h

        ind = round(h*step)+1;

        [ind ok] = spec_getnearest_max(S, ind, round(0.5*step));

        if ok
            [sins(1,h) sins(2,h)] = spec_fit_freq_amp(S, fs, ind, varargin{:});
        else
            freq = fs*(ind-1)/dftlen;
            amp = abs(S(ind));
        end
    end

    % Estimate the phase from linear interpolation of the phase spectrum
    sins(3,:) = wrap(interp1q(fs*(0:dftlen/2)'/dftlen, unwrap(angle(S(1:end/2+1))).', sins(1,:)'));

    sins(4,:) = 1:max_h;

    % Add the DC
    sins = [[0; abs(S(1)); angle(S(1)); 0], sins];

return

function [i ok] = spec_getnearest_max(S, i, limit)

    if i<=1 || i>=length(S)
        ok = false;
        return
    end

    ok = true;
    % get the nearest summit around i
    if ~(abs(S(i-1))<abs(S(i)) && abs(S(i))>abs(S(i+1)))
        il = i; ir = i;
        ilsum = false; irsum = false;
        imin = i-limit;
        imax = i+limit;
        while ok && ~ilsum && ~irsum
            il = il - 1;
            ok = il>1 && il>=imin;
            if ok; ilsum = abs(S(il-1))<abs(S(il)) && abs(S(il))>abs(S(il+1)); end
            ir = ir + 1;
            ok = ir<length(S) && ir<=imax;
            if ok; irsum = abs(S(ir-1))<abs(S(ir)) && abs(S(ir))>abs(S(ir+1)); end
        end
        if ok;
            if ilsum;       i = il;
            else            i = ir; end
        end
    end
return

