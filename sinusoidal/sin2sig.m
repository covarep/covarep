% Synthesis of a stationary time signal given sinusoids parameters
%
% Octave compatible
% 
% Inputs
%  sins   : matrix of size(3,N) [freq;amp;phase]
%  fs     : [Hz] Sampling frequency
%  winlen : length of the time signal to synthesize
%  derive : if deriv=1, synthesize the derivative of the signal
%
% Outputs
%  s       : The synthesized signal
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

function [s t] = sin2sig(sins, fs, winlen, deriv)
    if nargin<4; deriv=0; end

%      t = (0:winlen-1)/fs;
%      t = (0:winlen-1)/fs;
    % The time reference used for the phase in sins is in the window center
    t = (-(winlen-1)/2:(winlen-1)/2)/fs;

    if deriv==0
        A = cos(2*pi*t'*sins(1,:) + ones(length(t),1)*sins(3,:));
        A = ones(length(t),1)*sins(2,:).*A;
        s = sum(A,2);
    elseif deriv==1
        A = -sin(2*pi*t'*sins(1,:) + ones(length(t),1)*sins(3,:));
        A = ones(length(t),1)*(sins(2,:).*(2*pi*partials(1,:))).*A;
        s = sum(A,2);
    end

return
