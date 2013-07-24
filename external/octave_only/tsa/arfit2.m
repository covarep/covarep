function [w, MAR, C, sbc, fpe, th] = arfit2(Y, pmin, pmax, selector, no_const)
% ARFIT2 estimates multivariate autoregressive parameters
% of the MVAR process Y
%
%   Y(t,:)' = w' + A1*Y(t-1,:)' + ... + Ap*Y(t-p,:)' + x(t,:)'
%
% ARFIT2 uses the Nutall-Strand method (multivariate Burg algorithm) 
% which provides better estimates the ARFIT [1], and uses the 
% same arguments. Moreover, ARFIT2 is faster and can deal with 
% missing values encoded as NaNs. 
%
% [w, A, C, sbc, fpe] = arfit2(v, pmin, pmax, selector, no_const)
%
% INPUT: 
%  v		data - each channel in a column
%  pmin, pmax 	minimum and maximum model order
%  selector	'fpe' or 'sbc' [default] 
%  no_const	'zero' indicates no bias/offset need to be estimated 
%		in this case is w = [0, 0, ..., 0]'; 
%
% OUTPUT: 
%  w		mean of innovation noise
%  A		[A1,A2,...,Ap] MVAR estimates	
%  C		covariance matrix of innovation noise
%  sbc, fpe	criteria for model order selection 
%
% see also: ARFIT, MVAR
%
% REFERENCES:
%  [1] A. Schloegl, 2006, Comparison of Multivariate Autoregressive Estimators.
%       Signal processing, p. 2426-9.
%  [2] T. Schneider and A. Neumaier, 2001. 
%	Algorithm 808: ARFIT-a Matlab package for the estimation of parameters and eigenmodes 
%	of multivariate autoregressive models. ACM-Transactions on Mathematical Software. 27, (Mar.), 58-65.

%       $Id: arfit2.m 9609 2012-02-10 10:18:00Z schloegl $
%	Copyright (C) 1996-2005,2008,2012 by Alois Schloegl <alois.schloegl@ist.ac.at>	
%       This is part of the TSA-toolbox. See also 
%       http://pub.ist.ac.at/~schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


%%%%% checking of the input arguments was done the same way as ARFIT
if (pmin ~= round(pmin) || pmax ~= round(pmax))
        error('Order must be integer.');
end
if (pmax < pmin)
        error('PMAX must be greater than or equal to PMIN.')
end

% set defaults and check for optional arguments
if (nargin == 3)              	% no optional arguments => set default values
        mcor       = 1;         % fit intercept vector
        selector   = 'sbc';	% use SBC as order selection criterion
elseif (nargin == 4)          	% one optional argument
        if strcmp(selector, 'zero')
                mcor     = 0;   % no intercept vector to be fitted
                selector = 'sbc';	% default order selection 
        else
                mcor     = 1;	% fit intercept vector
        end
elseif (nargin == 5)		% two optional arguments
        if strcmp(no_const, 'zero')
                mcor     = 0;   	% no intercept vector to be fitted
        else
                error(['Bad argument. Usage: ', '[w,A,C,SBC,FPE,th]=AR(v,pmin,pmax,SELECTOR,''zero'')'])
        end
end



%%%%% Implementation of the MVAR estimation 
[N,M]=size(Y);
    
if mcor,
        [m,N] = sumskipnan(Y,1);                    % calculate mean 
        m = m./N;
	Y = Y - repmat(m,size(Y)./size(m));    % remove mean    
end;

[MAR,RCF,PE] = mvar(Y, pmax, 2);   % estimate MVAR(pmax) model

N = min(N);

%if 1;nargout>3;
ne = N-mcor-(pmin:pmax);
for p=pmin:pmax, 
        % Get downdated logarithm of determinant
        logdp(p-pmin+1) = log(det(PE(:,p*M+(1:M))*(N-p-mcor))); 
end;

% Schwarz's Bayesian Criterion
sbc = logdp/M - log(ne) .* (1-(M*(pmin:pmax)+mcor)./ne);

% logarithm of Akaike's Final Prediction Error
fpe = logdp/M - log(ne.*(ne-M*(pmin:pmax)-mcor)./(ne+M*(pmin:pmax)+mcor));

% Modified Schwarz criterion (MSC):
% msc(i) = logdp(i)/m - (log(ne) - 2.5) * (1 - 2.5*np(i)/(ne-np(i)));

% get index iopt of order that minimizes the order selection 
% criterion specified by the variable selector
if strcmpi(selector,'fpe'); 
    [val, iopt]  = min(fpe); 
else %if strcmpi(selector,'sbc'); 
    [val, iopt]  = min(sbc); 
end; 

% select order of model
popt = pmin + iopt-1; % estimated optimum order 

if popt<pmax, 
        [MAR, RCF, PE] = mvar(Y, popt, 2);
end;
%end

C = PE(:,size(PE,2)+(1-M:0));

if mcor,
        I = eye(M);        
        for k = 1:popt,
                I = I - MAR(:,k*M+(1-M:0));
        end;
        w = -I*m';
else
        w = zeros(M,1);
end;

if nargout>5,	th=[];	fprintf(2,'Warning ARFIT2: output TH not defined\n'); end;
