function R = histo2(Y, W)
% HISTO2 calculates histogram for multiple columns with separate bin values 
%    for each data column.
%
% R = HISTO2(Y)
% R = HISTO2(Y, W)
%	Y	data
%	W	weight vector containing weights of each sample, 
%		number of rows of Y and W must match.
%		default W=[] indicates that each sample is weighted with 1. 
%
% R = HISTO(...)            
% 	R is 	a struct with th fields 
%       R.X  	the bin-values, bin-values are computed separately for each 
%		data column, thus R.X is a matrix, each column contains the 
%		the bin values of for each data column, unused elements are indicated with NaN.
%		In order to have common bin values, use HISTO3.  
%       R.H  is the frequency of occurence of value X 
%  	R.N  are the number of valid (not NaN) samples (i.e. sum of weights)
%
% more histogram-based results can be obtained by HIST2RES2  
%
% see also: HISTO, HISTO2, HISTO3, HISTO4
%
% REFERENCE(S):
%  C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).

%	$Id: histo2.m 8383 2011-07-16 20:06:59Z schloegl $
%	Copyright (C) 1996-2002,2008,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%    	This is part of the TSA-toolbox 
%	http://hci.tugraz.at/~schloegl/matlab/tsa/
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


%%%%% check input arguments %%%%%
[yr,yc] = size(Y);
if nargin < 2, 
	W = []; 
end; 
if ~isempty(W) && (yr ~= numel(W)),
	error('number of rows of Y does not match number of elements in W');
end; 

%%%%% identify all possible X's and generate overall Histogram %%%%%
N  = sum(~isnan(Y), 1);
NN = N; 
if isempty(W)
	sY  = sort(Y,1);
else
	[sY, idx] = sort(Y,1);
	W = cumsum(W(idx));	%% W becomes cumulative sum 
end; 
[ix,iy] = find( diff(sY, [], 1) > 0);
nn0 = 0;

for k = 1:yc,
	tmp = [ix(iy==k); N(k)];
        nn1 = length(tmp);
        
        if isempty(W)
		H(1:nn1,k) = [tmp(1); diff(tmp)];
	else
		%%% Note that W is the cumulative sum
		H(1:nn1,k) = [W(tmp(1),k); diff(W(tmp,k))];
		NN(k) = W(N(k), k);	
	end;
	X(1:nn1, k) = sY(tmp, k);

        if k==1;
		nn0 = nn1;
	elseif nn1 < nn0,
		H (1+nn1:nn0, k) = NaN;
		X (1+nn1:nn0, k) = NaN;
        elseif nn1 > nn0,
		H (1+nn0:nn1, 1:k-1) = NaN;
		X (1+nn0:nn1, 1:k-1) = NaN;
		nn0 = nn1;
        end;
end;

R.datatype = 'HISTOGRAM';
R.H = H;
R.X = X;
R.N = NN;

