function Y=flix(D,x)
% floating point index - interpolates data in case of non-integer indices
%
% Y=flix(D,x)
%   FLIX returns Y=D(x) if x is an integer 
%   otherwise D(x) is interpolated from the neighbors D(ceil(x)) and D(floor(x)) 
% 
% Applications: 
% (1)  discrete Dataseries can be upsampled to higher sampling rate   
% (2)  transformation of non-equidistant samples to equidistant samples
% (3)  [Q]=flix(sort(D),q*(length(D)+1)) calculates the q-quantile of data series D   
%
% FLIX(D,x) is the same as INTERP1(D,X,'linear'); Therefore, FLIX might
% become obsolete in future. 
%
% see also: HIST2RES, Y2RES, PLOTCDF, INTERP1

%	$Id: flix.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) by 2001-2005,2008 Alois Schloegl <a.schloegl@ieee.org>	
%	This is part of the TSA-toolbox see also: 
% 	   http://www.dpmi.tu-graz.ac.at/schloegl/matlab/tsa/
%	   http://octave.sf.net/
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

D  = D(:);
Y  = x;

k1 = ((x >= 1) & (x <= size(D,1)));	
Y(~k1) = NaN;

k  = x - floor(x);	% distance to next sample	 

ix = ~k & k1; 	        % find integer indices
Y(ix) = D(x(ix)); 	% put integer indices

ix = k & k1;     	% find non-integer indices

Y(ix) = D(floor(x(ix))).*(1-k(ix)) + D(ceil(x(ix))).*k(ix);  

