function  [AR,RC] = arcext(MX,P);
% ARCEXT extracts AR and RC of order P from Matrix MX
% function  [AR,RC] = arcext(MX,P);
%
%  INPUT:
% MX 	AR and RC matrix calculated by durlev 
% P 	model order (default maximum possible)
%
%  OUTPUT
% AR    autoregressive model parameter	
% RC    reflection coefficients (= -PARCOR coefficients)
%
% All input and output parameters are organized in rows, one row 
% corresponds to the parameters of one channel
%
% see also ACOVF ACORF DURLEV 
% 
% REFERENCES:
%  P.J. Brockwell and R. A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  S. Haykin "Adaptive Filter Theory" 3rd ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.

%  $Id: arcext.m 9609 2012-02-10 10:18:00Z schloegl $
%  Copyright (C) 1998-2003,2008,2012 by Alois Schloegl <alois.schloegl@ist.ac.at>	
%       This is part of the TSA-toolbox. See also 
%       http://pub.ist.ac.at/~schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/
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


[lr,lc]=size(MX);

if ~mod(sqrt(8*lc+1)-1,2)	% full number of elements
	K = (sqrt(8*lc+1)-1)/2;
elseif ~mod(lc,2),   % compressed form of MX
	K = lc/2;
else		% invalid number of elements
	fprintf(2,'Warning ARCEXT: Number of elements is different than a triangular matrix\n');
end;

if (K~=P) && (lc~=K*(K+1)/2),
	[AR,RC,PE]=rc2ar(MX(:,(1:P).*(2:P+1)/2));
else
	AR = MX(:,P*(P-1)/2+(1:P));
	RC = MX(:,(1:P).*(2:P+1)/2);
end;
