function ACF=poly2ac(a,efinal)
% converts an AR polynomial into an autocorrelation sequence
% [R] = poly2ac(a [,efinal] );
%
% see also ACOVF ACORF AR2RC RC2AR DURLEV AC2POLY, POLY2RC, RC2POLY, RC2AC, AC2RC, POLY2AC
% 

%       $Id: poly2ac.m 5090 2008-06-05 08:12:04Z schloegl $
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


fprintf(2,'ERROR: POLY2AC does not work yet. Sorry\n');
return;

mfilename='POLY2AC';
if ~exist('rc2ar','file')
        fprintf(2,'Error %s: RC2AR.M not found. \n Download TSA toolbox from http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/\n',mfilename);
        return;
end;

[AR,RC,PE] = ar2rc(poly2ar(a));
%[AR,RC,PE,ACF] = rc2ar(RC);

if nargin<2, efinal=PE(:,size(PE,2)); end;

ACF=zeros(size(a));
ACF(:,1) = 1;
for k=2:size(a,2),
        ACF(:,k)=sum(AR(:,1:k-1).*ACF(:,k-(1:k-1)),2);
end;
R0=(sum(AR(:,1:k-1).*ACF(:,2:size(ACF,2)),2)+1).*efinal; %PE(:,size(PE,2));

ACF = ACF.*R0(:,ones(1,size(ACF,2)));

%ACF=ACF*efinal*PE(1)/PE(length(PE));
