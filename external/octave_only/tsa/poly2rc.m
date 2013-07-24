function [RC,r0] = poly2rc(a,efinal);
% converts AR-polynomial into reflection coefficients
% [RC,R0] = poly2rc(A [,Efinal])
%
%  INPUT:
% A     AR polynomial, each row represents one polynomial
% Efinal    is the final prediction error variance (default value 1)
%
%  OUTPUT
% RC    reflection coefficients
% R0    is the variance (autocovariance at lag=0) based on the 
%	prediction error
%
%
% see also ACOVF ACORF AR2RC RC2AR DURLEV AC2POLY, POLY2RC, RC2POLY, RC2AC, AC2RC, POLY2AC

%       $Id: poly2rc.m 5090 2008-06-05 08:12:04Z schloegl $
%       Copyright (C) 1998-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>
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

if all(size(a))>1,
        fprintf(2,'Error poly2rc: "a" must be a vector\n');
        return;
end;
a=a(:).';

mfilename='POLY2RC';
if ~exist('ar2rc','file')
        fprintf(2,'Error %s: AR2RC.M not found. \n Download TSA toolbox from http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/\n',mfilename);
        return;
end;

if nargin<2, efinal=1; end;

[AR,RC,PE] = ar2rc(poly2ar(a));
if nargout>1,
	r0=efinal.*PE(:,1)./PE(:,size(PE,2));
end;
