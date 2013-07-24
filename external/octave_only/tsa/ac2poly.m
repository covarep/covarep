function [A,E] = ac2poly(acf);
% converts the autocorrelation sequence into an AR polynomial
% [A,Efinal] = ac2poly(r)
%
% see also ACOVF ACORF AR2RC RC2AR DURLEV AC2POLY, POLY2RC, RC2POLY, RC2AC, AC2RC, POLY2AC
% 

%       $Id: ac2poly.m 5090 2008-06-05 08:12:04Z schloegl $
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


if all(size(acf))>1,
        fprintf(2,'Error AC2POLY: "r" must be a vector\n');
        return;
end;


mfilename='AC2POLY';
if ~exist('durlev','file')
        fprintf(2,'Error %s: DURLEV.M not found. \n Download TSA toolbox from http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/\n',mfilename);
        return;
end;

[AR,RC,PE] = durlev(acf);

A = [ones(size(AR,1),1),-AR];
E = PE(:,size(PE,2));
