function [A] = poly2ar(A);
% Converts AR polymials into autoregressive parameters. 
% Multiple polynomials can be converted. 
%
% function  [AR] = poly2ar(A);
%
%  INPUT:
% A     AR polynomial, each row represents one polynomial
%
%  OUTPUT
% AR    autoregressive model parameter	
%
% see also ACOVF ACORF DURLEV RC2AR AR2POLY

%       $Id: poly2ar.m 5090 2008-06-05 08:12:04Z schloegl $
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

% Inititialization
[lr,lc]=size(A);

if ~all(A(:,1)==1)
	fprintf(2,'Warning POLY2AR: input argument might not be an AR-polynom');
end;	

A = -A(:,2:size(A,2))./A(:,ones(1,size(A,2)-1));
