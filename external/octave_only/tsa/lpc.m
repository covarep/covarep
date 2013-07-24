function [A] = lpc(Y,P,mode);
% LPC Linear prediction coefficients 
% The Burg-method is used to estimate the prediction coefficients
%
% A = lpc(Y [,P]) finds the coefficients  A=[ 1 A(2) ... A(N+1) ],
%     	of an Pth order forward linear predictor
%     
% 	 Xp(n) = -A(2)*X(n-1) - A(3)*X(n-2) - ... - A(N+1)*X(n-P)
%	    
% 	such that the sum of the squares of the errors
%		
%       err(n) = X(n) - Xp(n)
%	       
%	is minimized.  X can be a vector or a matrix.  If X is a matrix
%       containing a separate signal in each column, LPC returns a model
%	estimate for each column in the rows of A. N specifies the order
%	of the polynomial A(z).
%				       
%	If you do not specify a value for P, LPC uses a default P = length(X)-1.
%
%
% see also ACOVF ACORF AR2POLY RC2AR DURLEV SUMSKIPNAN LATTICE 
% 

% REFERENCE(S):
%  J.P. Burg, "Maximum Entropy Spectral Analysis" Proc. 37th Meeting of the Society of Exp. Geophysiscists, Oklahoma City, OK 1967
%  J.P. Burg, "Maximum Entropy Spectral Analysis" PhD-thesis, Dept. of Geophysics, Stanford University, Stanford, CA. 1975.
%  P.J. Brockwell and R. A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  S.   Haykin "Adaptive Filter Theory" 3rd ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.

%	$Id: lpc.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>
%       This is part of the TSA-toolbox. See also 
%       http://hci.tugraz.at/schloegl/matlab/tsa/
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


[yr,yc] = size(Y);
if yr < yc,
	fprintf(2,'Warning LCP: data vector Y must be a column not a row vector\n');
end;

if nargin < 2,
        P = yr-1;
end;

% you can use any of the following routines. 
% the lattice methods are preferable for stochastic time series. 
% but can fail for deterministic signals see: 
% http://sourceforge.net/mailarchive/message.php?msg_name=20080516115110.GB20642%40localhost

% [AR,RC,PE] = lattice(Y.',P);		% Burg method
% [AR,RC,PE] = lattice(Y.',P,'GEOL');	% geometric lattice
[AR,RC,PE] = durlev(acovf(Y.',P));  	% Yule-Walker

A = ar2poly(AR);

