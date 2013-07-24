function [ACF,NN] = acovf(Z,KMAX,Mode,Mode2);
% ACOVF estimates autocovariance function (not normalized)
% NaN's are interpreted as missing values. 
%
% [ACF,NN] = acovf(Z,MAXLAG,Mode);
%
% Input:
%  Z    Signal (one channel per row);
%  MAXLAG  maximum lag
%  Mode	'biased'  : normalizes with N [default]
%	'unbiased': normalizes with N-lag
%	'coeff'	  : normalizes such that lag 0 is 1	
%        others	  : no normalization
%
% Output:
%  ACF autocovariance function
%  NN  number of valid elements 
%
% REFERENCES:
%  A.V. Oppenheim and R.W. Schafer, Digital Signal Processing, Prentice-Hall, 1975.
%  S. Haykin "Adaptive Filter Theory" 3ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.
%  J.S. Bendat and A.G.Persol "Random Data: Analysis and Measurement procedures", Wiley, 1986.

%	$Id: acovf.m 7687 2010-09-08 18:39:23Z schloegl $
%	Copyright (C) 1998-2003,2008,2010 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the TSA-toolbox. See also 
%       http://biosig-consulting.com/matlab/tsa/
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



if nargin<3, Mode='biased'; end;

[lr,lc] = size(Z);

MISSES = sum(isnan(Z)')';
if any(MISSES); % missing values
	M = real(~isnan(Z));
	Z(isnan(Z))=0;
end;

if (nargin == 1) 
	KMAX = lc-1; 
elseif (KMAX >= lc-1) 
	KMAX = lc-1;
end;

ACF = zeros(lr,KMAX+1);

if nargin>3,		% for testing, use arg4 for comparing the methods,

elseif 	(KMAX*KMAX > lc*log2(lc)) % & isempty(MISSES);	
	Mode2 = 1;
elseif 	(10*KMAX > lc);
	Mode2 = 3;
else
	Mode2 = 4;
end;


%%%%% ESTIMATION of non-normalized ACF %%%%%

% the following algorithms gve equivalent results, however, the computational effort is different,
% depending on lr,lc and KMAX, a different algorithm is most efficient.
if Mode2==1; % KMAX*KMAX > lc*log(lc);        % O(n.logn)+O(K²)
        tmp = fft(Z',2^nextpow2(size(Z,2))*2);
        tmp = ifft(tmp.*conj(tmp));
        ACF = tmp(1:KMAX+1,:)'; 
        if ~any(any(imag(Z))), ACF=real(ACF); end; % should not be neccessary, unfortunately it is.
elseif Mode2==3; % (10*KMAX > lc)   % O(n*K)     % use fast Built-in filter function
        for L = 1:lr,
                acf = filter(Z(L,lc:-1:1),1,Z(L,:));
                ACF(L,:)= acf(lc:-1:lc-KMAX);
        end;    
else Mode2==4; % O(n*K)
        for L = 1:lr,
                for K = 0:KMAX, 
                        ACF(L,K+1) = Z(L,1:lc-K) * Z(L,1+K:lc)';
                end;
        end;    
end;


%%%%% GET number of elements used for estimating ACF - is needed for normalizing ACF %%%%%

if any(MISSES),
    % the following algorithms gve equivalent results, however, the computational effort is different,
    % depending on lr,lc and KMAX, a different algorithm is most efficient.
    if Mode2==1; % KMAX*KMAX > lc*log(lc);        % O(n.logn)+O(K²)
        tmp = fft(M',2^nextpow2(size(M,2))*2);
        tmp = ifft(tmp.*conj(tmp));
        NN = tmp(1:KMAX+1,:)'; 
        if ~any(any(imag(M))), NN=real(NN); end; % should not be neccessary, unfortunately it is.
    elseif Mode2==3; % (10*KMAX > lc)   % O(n*K)     % use fast Built-in filter function
        for L = 1:lr,
                acf = filter(M(L,lc:-1:1),1,M(L,:));
                NN(L,:)= acf(lc:-1:lc-KMAX);
        end;    
    else Mode2==4; % O(n*K)
        for L = 1:lr,
                for K = 0:KMAX, 
                        NN(L,K+1) = M(L,1:lc-K) * M(L,1+K:lc)';
                end;
        end;    
    end;
else
    NN = (ones(lr,1)*(lc:-1:lc-KMAX));
end;


if strcmp(Mode,'biased')
	if ~any(MISSES),
	        ACF=ACF/lc;
	else
	        %ACF=ACF./((lc-MISSES)*ones(1,KMAX+1));
	        ACF=ACF./max(NN + ones(lr,1)*(0:KMAX),0);
	end;

elseif strcmp(Mode,'unbiased')
        ACF=ACF./NN; 
	%if ~any(MISSES),
	%       ACF=ACF./(ones(lr,1)*(lc:-1:lc-KMAX));
	%else
	%	ACF=ACF./((lc-MISSES)*ones(1,KMAX+1) - ones(lr,1)*(0:KMAX));
	%end;

elseif strcmp(Mode,'coeff')
        %ACF = ACF ./ ACF(:,ones(1,KMAX+1)) .* ((lc-MISSES)*ones(1,KMAX+1));
        ACF = ACF./NN; 
	ACF = ACF./(ACF(:,1)*ones(1,size(ACF,2)));
else 

end;
