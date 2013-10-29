% STE_WLP Short-time-energy-weighted linear prediction
%   [A,w] = ste_wlp(s,p,m,k)
%
% Description
%   This function fits linear prediction coefficients to the analysis
%   frame using the basic form of weighted linear prediction (WLP), i.e., 
%   by weighting each value of the squared prediction error by the short-time 
%   energy (STE) of the previous samples.
%
% Inputs
%   s      : Speech signal frame [samples]
%   p      : Order of WLP analysis
%   m      : Length of the STE window (default to p)
%   k      : Delay of the STE window (default to 1)
%
% Outputs
%   A      : Linear prediction inverse filter coefficients
%   w      : the STE weighting function
%
% Notes
%   Even though coefficients are solved using a criterion similar to the
%   autocorrelation method of linear prediction, because of the weighting
%   the filter is not guaranteed to be stable. If a stable synthesis filter
%   is required, use STE_SWLP instead.
%
% Example
%   A = ste_wlp(s,p) gives the linear predictive inverse filter
%   coefficients optimized using STE-WLP
%
% References
%  [1] C. Ma, Y. Kamp and L. F. Willems, "Robust signal selection for
%  linear prediction analysis of voiced speech", Speech Communication, vol.
%  12, no. 1, pp. 6981, 1993.
%
% Copyright (c) 2013 Aalto University
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
% This function is part of the Common Speech Processing Repository 
% http://covarep.github.io/covarep/
%
% Octave Compatible
% 
% Author
%  Jouni Pohjalainen jouni.pohjalainen@aalto.fi
%
% $Id <info set by the versioning system> $

function [A,w] = ste_wlp(s,p,m,k)

% handle the special case of all-zero frames
if all(s==0)
  s = eps*randn(size(s));
end

N = length(s);

% default value for the STE window length
if nargin<3
    m = p;
end

% default value for the STE window lag
if nargin<4
    k = 1;
end

% compute STE weighting function
w = stew(s,m,p,k)+eps;

% apply the square root of the weighting function to the delayed versions
% of the signal from lag 0 to lag p
wsr = sqrt(w(1:(N+p)));
Y = zeros(N+p,p+1);
for i1=0:p
    Y(:,i1+1) = [zeros(i1,1);s;zeros(p-i1,1)].*wsr;
end

% compute weighted autocorrelations
R = (Y'*Y)/N;

% solve the p predictor coefficients
A = R(2:end,2:end)\R(2:end,1);

% convert to inverse filter form
A = [1;-A]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = stew(x,m,p,k)
% STE weighting for an order p predictor using window of length m delayed
% by k samples
w = conv([zeros(k,1);x(:)].^2,ones(m,1));
w = [w;zeros(p-m,1)];
