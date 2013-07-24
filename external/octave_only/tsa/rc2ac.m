function ACF=rc2ac(RC,R0)
% converts reflection coefficients to autocorrelation sequence
% [R] = rc2ac(RC,R0);
%
% see also ACOVF ACORF AR2RC RC2AR DURLEV AC2POLY, POLY2RC, RC2POLY, RC2AC, AC2RC, POLY2AC
% 

%       $Id: rc2ac.m 5090 2008-06-05 08:12:04Z schloegl $
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


fprintf(2,'ERROR: RC2AC does not work yet. Sorry\n');
return;

if all(size(RC)>1),
        fprintf(2,'Error RC2AC: "K" must be a vector\n');
        return;
end;

mfilename='RC2AC';
if ~exist('rc2ar','file')
        fprintf(2,'Error %s: RC2AR.M not found. \n Download TSA toolbox from http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/\n',mfilename);
        return;
end;

[AR,RC,PE] = rc2ar(RC(:).');

ACF=cumprod(ones(size(RC))-RC.^2,2);
ACF = ACF./ACF(:,ones(1,size(ACF,2)));

if nargin>1
        ACF=ACF*R0;
end;
