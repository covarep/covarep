% Estimate of harmonic frequencies, amplitude and phase given a gross f0
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
%  [max_h]      : number of harmonic to return (+1 bcs of DC)
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

function sins = spec_getsins_f0(S, fs, f0, max_h)

    if nargin<4 || isempty(max_h); max_h = round((fs/2-f0/2)/f0);
    else;                          max_h = min(max_h,floor((fs/2-f0/2)/f0)); end

    S = S(:).';

    dftlen = length(S);

    step = dftlen*f0/fs; % f0 in number of bins

    sins = zeros(5,max_h); % Allocate the whole matrix first

    % By default, use a simple harmonic sampling ...
    sins(1,:) = step*(1:max_h);
    sins(2,:) = abs(S(1+round(sins(1,:))));

    % ... and replace by the peaks found.
    [k, v] = v_findpeaks(log(abs(S(1:end/2+1))), 'q');
    if length(k)==1
        idx = abs(k-(1+sins(1,:)))<step/2;
        sins(1,idx) = k-1;
        sins(2,idx) = exp(v);
        sins(5,idx) = true;
    elseif length(k)>1
        D = abs(repmat(k,1,max_h)-repmat(1+sins(1,:), length(k),1));
        [mind, mindi] = min(D);
        idx = mind<step/2;
        sins(1,idx) = k(mindi(idx))-1;
        sins(2,idx) = exp(v(mindi(idx)));
        sins(5,idx) = true;
    end
    sins(1,:) = (fs/dftlen)*sins(1,:);

    % Estimate the phase from linear interpolation of the phase spectrum
    sins(3,:) = wrap(interp1q(fs*(0:dftlen/2)'/dftlen, unwrap(angle(S(1:end/2+1))).', sins(1,:)'));

    % The harmonic number
    sins(4,:) = 1:size(sins,2);

    % Add the DC
    sins = [[0; abs(S(1)); angle(S(1)); 0; abs(S(1))>abs(S(2))], sins];

    if 0
        F = fs*(0:dftlen/2)/dftlen;
        plot(F, mag2db(abs(S(1:end/2+1))), 'k');
        hold on;
        plot(sins(1,:), mag2db(sins(2,:)), 'xb');
        idx = find(sins(5,:));
        plot(sins(1,idx), mag2db(sins(2,idx)), 'xr');
        keyboard
    end

return
