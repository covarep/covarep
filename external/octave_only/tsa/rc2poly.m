function [a,efinal] = rc2poly(RC,E);
% converts reflection coefficients into an AR-polynomial
% [a,efinal] = rc2poly(K)
%
% see also ACOVF ACORF AR2RC RC2AR DURLEV AC2POLY, POLY2RC, RC2POLY, RC2AC, AC2RC, POLY2AC
% 

%       $Id: rc2poly.m 9609 2012-02-10 10:18:00Z schloegl $
%       Copyright (C) 1998-2002,2008,2012 by Alois Schloegl <alois.schloegl@ist.ac.at>
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


if nargout>1 && nargin<2
   fprintf('Zero-lag autocorrelation, R0 not specified\n')
   return; 
end;
   
mfilename='RC2POLY';
if all(size(RC))>1,
        fprintf(2,'Error %s: "K" must be a vector\n',mfilename);
        return;
end;

if ~exist('rc2ar','file')
        fprintf(2,'Error %s: RC2AR.M not found. \n Download TSA toolbox from http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/\n',mfilename);
        return;
end;

[AR,RC,PE] = rc2ar(RC(:).');

a=[1,-AR];
efinal=PE(length(PE))/PE(1);
