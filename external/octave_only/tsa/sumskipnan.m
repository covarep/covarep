function [o,count,SSQ] = sumskipnan(x, DIM, W)
% SUMSKIPNAN adds all non-NaN values. 
%
% All NaN's are skipped; NaN's are considered as missing values. 
% SUMSKIPNAN of NaN's only  gives O; and the number of valid elements is return. 
% SUMSKIPNAN is also the elementary function for calculating 
% various statistics (e.g. MEAN, STD, VAR, RMS, MEANSQ, SKEWNESS, 
% KURTOSIS, MOMENT, STATISTIC etc.) from data with missing values.  
% SUMSKIPNAN implements the DIMENSION-argument for data with missing values.
% Also the second output argument return the number of valid elements (not NaNs) 
% 
% Y = sumskipnan(x [,DIM])
% [Y,N,SSQ] = sumskipnan(x [,DIM])
% [...] = sumskipnan(x, DIM, W)
% 
% x	input data 	
% DIM	dimension (default: [])
%	empty DIM sets DIM to first non singleton dimension	
% W	weight vector for weighted sum, numel(W) must fit size(x,DIM)
% Y	resulting sum
% N	number of valid (not missing) elements
% SSQ	sum of squares
%
% the function FLAG_NANS_OCCURED() returns whether any value in x
%  is a not-a-number (NaN)
%
% features:
% - can deal with NaN's (missing values)
% - implements dimension argument. 
% - computes weighted sum 
% - compatible with Matlab and Octave
%
% see also: FLAG_NANS_OCCURED, SUM, NANSUM, MEAN, STD, VAR, RMS, MEANSQ, 
%      SSQ, MOMENT, SKEWNESS, KURTOSIS, SEM


%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.

%	$Id: sumskipnan.m 9033 2011-11-08 20:58:07Z schloegl $
%    	Copyright (C) 2000-2005,2009,2011 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


global FLAG_NANS_OCCURED;

if nargin<2,
        DIM = [];
end;
if nargin<3,
        W = [];
end;

% an efficient implementation in C of the following lines 
% could significantly increase performance 
% only one loop and only one check for isnan is needed
% An MEX-Implementation is available in sumskipnan.cpp
%
% Outline of the algorithm: 
% for { k=1,o=0,count=0; k++; k<N} 
% 	if ~isnan(i(k)) 
% 	{ 	o     += x(k);
%               count += 1;
%		tmp    = x(k)*x(k)
%		o2    += tmp;
%		o3    += tmp.*tmp;
%       }; 

if isempty(DIM),
        DIM = find(size(x)>1,1);
        if isempty(DIM), DIM = 1; end;
end
if (DIM<1), DIM = 1; end; %% Hack, because min([])=0 for FreeMat v3.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-float data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  (isempty(W) && (~(isa(x,'float') || isa(x,'double')))) || ~flag_implicit_skip_nan(), %%% skip always NaN's
	if ~isempty(W)
		error('SUMSKIPNAN: weighted sum of integers not supported, yet');
	end; 
	x = double(x); 
	o = sum(x,DIM);
	if nargout>1
		sz = size(x);
		N  = sz(DIM); 
		sz(DIM) = 1; 	
		count = repmat(N,sz);
		if nargout>2
			x = x.*x; 
			SSQ = sum(x,DIM);
		end; 
	end; 	
	return; 
end; 	

if (length(size(x))<DIM)
	error('SUMSKIPNAN: DIM argument larger than number of dimensions of x');
elseif ~isempty(W) && (size(x,DIM)~=numel(W))
	error('SUMSKIPNAN: size of weight vector does not match size(x,DIM)');
end; 

%% mex and oct files expect double
x = double(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use Matlab-MEX function when available  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if 1,	
try
	
	%% using sumskipnan_mex.mex

	%% !!! hack: FLAG_NANS_OCCURED is an output argument, reserve memory !!!
	if isempty(FLAG_NANS_OCCURED),
		FLAG_NANS_OCCURED = logical(0);  % default value 
	end;

	if (nargout<2),
		o = sumskipnan_mex(real(x),DIM,FLAG_NANS_OCCURED,W);
		if (~isreal(x))
			io = sumskipnan_mex(imag(x),DIM,FLAG_NANS_OCCURED,W);
			o  = o + i*io;
		end; 
		return; 
	elseif (nargout==2),
		[o,count] = sumskipnan_mex(real(x),DIM,FLAG_NANS_OCCURED,W);
		if (~isreal(x))
			[io,icount] = sumskipnan_mex(imag(x),DIM,FLAG_NANS_OCCURED,W);
			if any(count(:)-icount(:))
				error('Number of NaNs differ for REAL and IMAG part');
			else
				o  = o+i*io;
			end; 
		end; 
		return; 
	elseif (nargout>=3),
		[o,count,SSQ] = sumskipnan_mex(real(x),DIM,FLAG_NANS_OCCURED,W);
		if (~isreal(x))
			[io,icount,iSSQ] = sumskipnan_mex(imag(x),DIM,FLAG_NANS_OCCURED,W);
			if any(count(:)-icount(:))
				error('Number of NaNs differ for REAL and IMAG part');
			else
				o  = o+i*io;
				SSQ = SSQ+iSSQ;
			end; 
		end; 
		return; 
	end; 	
end; 

if ~isempty(W) 
	error('weighted sumskipnan requires sumskipnan_mex');
end; 	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% count non-NaN's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>1,
        count = sum(x==x,DIM); 
	FLAG_NANS_OCCURED = any(count(:)<size(x,DIM));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace NaN's with zero 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(x~=x) = 0;	
o = sum(x,DIM);

if nargout>2,
        x = real(x).^2 + imag(x).^2;
        SSQ = sum(x,DIM);
end;

%!assert(sumskipnan([1,2],1),[1,2])
%!assert(sumskipnan([1,NaN],2),1)
%!assert(sumskipnan([1,NaN],2),1)
%!assert(sumskipnan([nan,1,4,5]),10)
%!assert(sumskipnan([nan,1,4,5]',1,[3;2;1;0]),6)



