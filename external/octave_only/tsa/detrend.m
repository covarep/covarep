function [X,T]=detrend(t,X,p)
% DETREND removes the trend from data, NaN's are considered as missing values
% 
% DETREND is fully compatible to previous Matlab and Octave DETREND with the following features added:
% - handles NaN's by assuming that these are missing values
% - handles unequally spaced data
% - second output parameter gives the trend of the data
% - compatible to Matlab and Octave 
%
% [...]=detrend([t,] X [,p])
%	removes trend for unequally spaced data
%	t represents the time points
%	X(i) is the value at time t(i)
%	p must be a scalar
%
% [...]=detrend(X,0)
% [...]=detrend(X,'constant')
%	removes the mean
%
% [...]=detrend(X,p)
%	removes polynomial of order p (default p=1)
%
% [...]=detrend(X,1) - default
% [...]=detrend(X,'linear')
%	removes linear trend 
%
% [X,T]=detrend(...) 
%
% X is the detrended data
% T is the removed trend
% 
% see also: SUMSKIPNAN, ZSCORE		

% Copyright (C) 1995, 1996 Kurt Hornik <Kurt.Hornik@ci.tuwien.ac.at>
%       $Id: detrend.m 5090 2008-06-05 08:12:04Z schloegl $
%       Copyright (C) 2001,2007 by Alois Schloegl <a.schloegl@ieee.org>	
%       This function is part of the TSA-toolbox
%       http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/
%	Copyright (C) 1997, 1998, 2008 by Alois Schloegl <a.schloegl@ieee.org>
%
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


if (nargin == 1)
    	p = 1;
       	X = t;
	t = [];
elseif (nargin == 2)
	if strcmpi(X,'constant'),
		p = 0; 
		X = t; 
		t = [];
	elseif strcmpi(X,'linear'),
		p = 1; 
		X = t; 
		t = [];
	elseif ischar(X)
		error('unknown 2nd input argument');	
        elseif all(size(X)==1), 
                p = X;
                X = t;
                t = [];
        else
                p = 1;
    	end;            
elseif (nargin == 3)
        if ischar(X),
		warning('input arguments are not supported');	
        end; 
        
elseif (nargin > 3)
    	fprintf (1,'usage: detrend (x [, p])\n');
end;

% check data, must be in culomn order
[m, n] = size (X);
if (m == 1)
        X = X';
        r=n;
else
        r=m;
end
% check time scale
if isempty(t),
	t = (1:r).'; % make time scale 
elseif ~all(size(t)==size(X)) 
        t = t(:);
end;
% check dimension of t and X
if ~all(size(X,1)==size(t,1))
        fprintf (2,'detrend: size(t,1) must same as size(x,1) \n');
end;
% check the order of the polynomial 
if (~(all(size(p)==1) & (p == round (p)) & (p >= 0)))
	fprintf (2,'detrend:  p must be a nonnegative integer\n');
end

if (nargout>1)  , % needs more memory
        T = zeros(size(X))+nan; 
        %T=repmat(nan,size(X)); % not supported by Octave 2.0.16
        
        
        if (size(t,2)>1),	% for multiple time scales
                for k=1:size(X,2),
	                idx=find(~isnan(X(:,k)));
                        b = (t(idx,k) * ones (1, p + 1)) .^ (ones (length(idx),1) * (0 : p));
		        T(idx,k) = b * (b \ X(idx,k));
		end;
        	        
        else			% if only one time scale is used
		b = (t * ones (1, p + 1)) .^ (ones (length(t),1) * (0 : p));
                for k=1:size(X,2),
	                idx=find(~isnan(X(:,k)));
		        T(idx,k) = b(idx,:) * (b(idx,:) \ X(idx,k));
	        	%X(idx,k) = X(idx,k) - T(idx,k); % 1st alternative implementation
                        %X(:,k) = X(:,k) - T(:,k); % 2nd alternative 
                end;
	end;
        X = X-T;  % 3nd alternative 
        
        if (m == 1)
	        X = X';
	        T = T';
	end
else % needs less memory
        if (size(t,2)>1),	% for multiple time scales
                for k = 1:size(X,2),
	                idx = find(~isnan(X(:,k)));
                        b = (t(idx,k) * ones (1, p + 1)) .^ (ones (length(idx),1) * (0 : p));
		        X(idx,k) = X(idx,k) -  b * (b \ X(idx,k));
		end;
        else			% if only one time scale is used
		b = (t * ones (1, p + 1)) .^ (ones (length(t),1) * (0 : p));
		for k = 1:size(X,2),
		        idx = find(~isnan(X(:,k)));
	        	X(idx,k) = X(idx,k) - b(idx,:) * (b(idx,:) \ X(idx,k));
		end;
	end;

        if (m == 1)
	        X = X';
	end
end;



