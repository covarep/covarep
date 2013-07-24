function [AutoCov,AutoCorr,MX,E,NC]=invest0(Y,Pmax,Mode);
% First Investigation of a signal (time series) - automated part
% [AutoCov,AutoCorr,ARPMX,E,ACFsd,NC]=invest0(Y,Pmax);
%
% [AutoCov,AutoCorr,ARPMX,E,ACFsd,NC]=invest0(AutoCov,Pmax,Mode);
% 
%
% Y	time series
% Pmax	maximal order (optional)
%
% AutoCov	Autocorrelation 
% AutoCorr	normalized Autocorrelation
% PartACF	Partial Autocorrelation
% ARPMX     Autoregressive Parameter for order Pmax-1
% E	        Error function E(p)
% NC            Number of values (length-missing values)
%
% REFERENCES:
%  P.J. Brockwell and R.A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  M.S. Grewal and A.P. Andrews "Kalman Filtering" Prentice Hall, 1993. 
%  S. Haykin "Adaptive Filter Theory" 3ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.

%	$Id: invest0.m 7687 2010-09-08 18:39:23Z schloegl $
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

if nargin<3
   Mode=0;
else
   Mode=1;
end;	

[nr,nc]=size(Y);
NC = sumskipnan(real(~isnan(Y)),2);             % number of valid components (data points)

if Mode==0
	if nargin<2, Pmax = min([100 nc/3]); end;
	M = min(Pmax,nc-1);
        AutoCov = acovf(Y,M);
else
        AutoCov=Y;
        M  = min(Pmax,nc-1);
        nc = Pmax;
        %M=size(AutoCov,2)-1;
end;

AutoCorr = AutoCov(:,2:M+1)./AutoCov(:,ones(M,1));

if 1,Pmax<100; % this needs change of sinvest1
	K=M-1;
	[MX,E]=lattice(Y,Pmax);
	%[MX,E]=ulsar(Y,Pmax);
	%[MX,E]=durlev(AutoCov);
%    	ARP=MX(:,K*(K-1)/2+(1:K));
%    	PartACF =MX(:,(1:K).*(2:K+1)/2);
else %if nargout > 2
	[ARP,RC,E]=lattice(Y,Pmax);
        %[ARP,PartACF,E]=durlev(AutoCov);
	%MX=[ARP,RC];        
end;



                    
                    
                    
                    
                    
                    
