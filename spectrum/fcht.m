% The Fan-Chirp Transform (FChT)
% 
% This function is Octave compat., but dreadfully slow, because of the exp.
%
% Description
%  This is an implementation of the Fan-Chirp Transform [1][2].
%
% Inputs
%  x       : A short time signal (usually windowed)
%  ahat    : a/fs where a=f0'/f0
%
% Outputs
%  X       : The FChT transform.
%
% Example
%  Please see the HOWTO_spectra
%
% References
%  [1] Kepesi, M., Weruaga, L.: Adaptive Chirp-based time-frequency analysis of
%      speech signals, Speech communication 48(5), 474–492, 2006.
%  [2] Weruaga, L, Kepesi, M: The Fan-Chirp Transform for non-stationary harmonic
%      signals, Signal Processing 87(6), 1504–1522, 2007.
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

function X = fcht(x, ahat, N, ks)

    x = x(:);
    M = length(x);

    if abs(ahat)>2/M; warning('The focal point of the FChT is inside the window.'); end % (63)

    if nargin<3;
		N=length(x);
	else
		if N>M
			% Nothing to do
		elseif N<M
			warning('The requested spectrum size is smaller than the provided signal. Therefore, the signal will be truncated.');
			x = x(1:N);
		end
	end

	% If ks is not provided, compute the whole FCHT
	% (over all possible discrete frequencies)
	if nargin<4
		if mod(N,2)==0    % if even
			ks = (-N/2+1:N/2)';
			ks = [ks(N/2:end); ks(1:N/2-1)];           % same as in FFT
		else              % else odd
			ks = (-(N-1)/2:(N-1)/2)';
			ks = [ks((N-1)/2+1:end); ks(1:(N-1)/2)];   % same as in FFT
		end
	end

    n = -(M-1)/2:(M-1)/2;

    E = (-2*pi/N).*ks*((1+0.5*ahat*n).*n);
    E = exp(j*E); % This kills Octave
    X = sum(ones(N,1)*(x'.*sqrt(abs((1+ahat*n)))).*E,2);

end
