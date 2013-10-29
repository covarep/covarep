% STE_SWLP Stabilized weighted linear prediction using short-time-energy
% weighting
%   [A,w] = ste_swlp(s,p,m,k)
%
% Description
%   This function fits linear prediction coefficients to the analysis
%   frame using the basic form of stabilized weighted linear prediction
%   (SWLP) that weights each value of the squared prediction error by the
%   short-time energy (STE) of the previous samples.
%
% Inputs
%   s      : Speech signal frame [samples]
%   p      : Order of SWLP analysis
%   m      : Length of the STE window (default to p)
%   k      : Delay of the STE window (default to 1)
%
% Outputs
%   A      : Linear prediction inverse filter coefficients
%   w      : the STE weighting function
%
% Notes
%   SWLP generally gives smoother spectral shapes than the corresponding
%   weighted linear prediction (WLP) using the same weighting function.
%
% Example
%   A = ste_swlp(s,p) gives the linear predictive inverse filter
%   coefficients optimized using STE-SWLP
%
% References
%  [1] C. Magi, J. Pohjalainen, T. Bäckström and P. Alku, "Stabilised
%  weighted linear prediction", Speech Communication, vol. 51, no. 5, pp.
%  401411, 2009.
%  [2] J. Pohjalainen, H. Kallasjoki, K. J. Palomäki, M. Kurimo and P.
%  Alku, "Weighted linear prediction for speech analysis in noisy
%  conditions", in Proc. Interspeech, Brighton, UK, pp. 1315-1318,
%  2009.
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

function [A,w] = ste_swlp(s,p,m,k)

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

w = w(1:(N+p));

% initialize partial weights for recursive computation (see [2], Eqs.
% (5)-(7), for the partial-weight formulation of SWLP)
Z = zeros(N+p,p+1);
Z(:,1) = sqrt(w);

% delayed and weighted versions of the signal
Y = zeros(N+p,p+1);
Y(:,1) = Z(:,1).*[s;zeros(p,1)];

% recursion for partial weights and weighting of differently lagged
% versions of the signal
for i1=1:p
    Z((i1+1):(N+p),i1+1) = max([ones(N+p-i1,1) sqrt(w((i1+1):(N+p))./w(i1:(N+p-1)))],[],2) .* Z(i1:(N+p-1),i1);
    Y(:,i1+1) = Z(:,i1+1) .* [zeros(i1,1);s;zeros(p-i1,1)];
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
