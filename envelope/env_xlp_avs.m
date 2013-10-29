% AVS_XLP Extended weighted linear prediction using absolute-value-sum
% weighting
%   [A,Z] = avs_xlp(s,p,m)
%
% Description
%   This function fits linear prediction coefficients to the analysis
%   frame using the initial version of extended weighted linear prediction
%   (XLP) which weights each lagged sample at each prediction instant by
%   the absolute value sum (AVS) scheme, i.e., a low-pass-filtered sum of
%   absolute values of the predicted and lagged sample.
%
% Inputs
%   s      : Speech signal frame [samples]
%   p      : Order of XLP analysis
%   m      : Comparable fixed window for first-order IIR lowpass (whose
%   memory coefficient is obtained as (m-1)/m) which is used to filter
%   the lag weights (defaults to p)
%
% Outputs
%   A      : Linear prediction inverse filter coefficients
%   Z      : the partial-weight matrix of size (p+1) times (N+p), where N
%   is the length of s
%
% Notes
%   Extended weighted linear prediction (XLP) generalizes weighted linear
%   prediction (WLP). The particular weighting scheme used in this version
%   of XLP has shown improved robustness of feature extraction in
%   some studies on speaker verification and automatic speech recognition
%   [1][2]. It is not guaranteed to produce a stable synthesis filter but
%   appears to produce a stable filter more often than the basic version of
%   WLP (STE-WLP).
%
% Example
%   A = avs_xlp(s,p) gives the linear predictive inverse filter
%   coefficients optimized using AVS-XLP
%
% References
%  [1] J. Pohjalainen, R. Saeidi, T. Kinnunen and P. Alku, "Extended
%  weighted linear prediction (XLP) analysis of speech and its application
%  to speaker verification in adverse conditions", in Proc. Interspeech,
%  Makuhari, Japan, 2010.
%  [2] S. Keronen, J. Pohjalainen, P. Alku and M. Kurimo, "Noise robust
%  feature extraction based on extended weighted linear prediction in
%  LVCSR", in Proc. Interspeech, Florence, Italy, 2011.
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
% 
% 
% Author
%  Jouni Pohjalainen jouni.pohjalainen@aalto.fi
%
% $Id <info set by the versioning system> $

function [A,Z] = avs_xlp(s,p,m,varargin)

% handle the special case of all-zero frames
if all(s==0)
  s = eps*randn(size(s));
end

N = length(s);

% default value for the parameter determining the memory coefficient of
% first-order IIR low-pass filtering used in computing the weights
if nargin<3
  m = p;
end
 
% partial-weight matrix initialized with zeros
if nargout>1
  Z = zeros(p+1,N+p);
end

% p+1 lagged signal samples for each time instant
% use the buffer function from Matlab's Signal processing toolbox
D = flipud(buffer([s;zeros(p,1)],p+1,p));

% AVS weighting
ac = zeros(p+1,1);
% memory coefficient for lowpass filtering
mem = (m-1)./m;
for i1=1:N+p
    ac = mem.*ac + (1-mem).*(abs(D(1,i1))+abs(D(:,i1)));
    D(:,i1) = D(:,i1) .* ac;
    % if the weighting matrix is required as output
    if nargout>1
        Z(:,i1) = ac;
    end
end

% compute weighted autocorrelations
R = (D*D')/N;

% solve the p predictor coefficients
A = R(2:end,2:end)\R(2:end,1);

% convert to inverse filter form
A = [1;-A]';
