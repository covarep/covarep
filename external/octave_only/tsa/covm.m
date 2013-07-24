function [CC,NN] = covm(X,Y,Mode,W)
% COVM generates covariance matrix
% X and Y can contain missing values encoded with NaN.
% NaN's are skipped, NaN do not result in a NaN output. 
% The output gives NaN only if there are insufficient input data
%
% COVM(X,Mode);
%      calculates the (auto-)correlation matrix of X
% COVM(X,Y,Mode);
%      calculates the crosscorrelation between X and Y
% COVM(...,W);
%	weighted crosscorrelation 
%
% Mode = 'M' minimum or standard mode [default]
% 	C = X'*X; or X'*Y correlation matrix
%
% Mode = 'E' extended mode
% 	C = [1 X]'*[1 X]; % l is a matching column of 1's
% 	C is additive, i.e. it can be applied to subsequent blocks and summed up afterwards
% 	the mean (or sum) is stored on the 1st row and column of C
%
% Mode = 'D' or 'D0' detrended mode
%	the mean of X (and Y) is removed. If combined with extended mode (Mode='DE'), 
% 	the mean (or sum) is stored in the 1st row and column of C. 
% 	The default scaling is factor (N-1). 
% Mode = 'D1' is the same as 'D' but uses N for scaling. 
%
% C = covm(...); 
% 	C is the scaled by N in Mode M and by (N-1) in mode D.
% [C,N] = covm(...);
%	C is not scaled, provides the scaling factor N  
%	C./N gives the scaled version. 
%
% see also: DECOVM, XCOVF

%	$Id: covm.m 9032 2011-11-08 20:25:36Z schloegl $
%	Copyright (C) 2000-2005,2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

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


global FLAG_NANS_OCCURED;

if nargin<3,
	W = []; 
        if nargin==2,
		if isnumeric(Y),
			Mode='M';
		else
			Mode=Y;
			Y=[];
		end;
        elseif nargin==1,
                Mode = 'M';
                Y = [];
        elseif nargin==0,
                error('Missing argument(s)');
        end;

elseif (nargin==3) && isnumeric(Y) && ~isnumeric(Mode);
	W = [];

elseif (nargin==3) && ~isnumeric(Y) && isnumeric(Mode);
	W = Mode; 
	Mode = Y;
	Y = [];

elseif (nargin==4) && ~isnumeric(Mode) && isnumeric(Y);
	; %% ok 
else 
	error('invalid input arguments');
end;

Mode = upper(Mode);

[r1,c1]=size(X);
if ~isempty(Y)
        [r2,c2]=size(Y);
        if r1~=r2,
                error('X and Y must have the same number of observations (rows).');
        end;
else
        [r2,c2]=size(X);
end;

persistent mexFLAG2; 
persistent mexFLAG; 
if isempty(mexFLAG2) 
	mexFLAG2 = exist('covm_mex','file');	
end; 
if isempty(mexFLAG) 
	mexFLAG = exist('sumskipnan_mex','file');	
end; 


if ~isempty(W)
	W = W(:);
	if (r1~=numel(W))
		error('Error COVM: size of weight vector does not fit number of rows');
	end;
	%w = spdiags(W(:),0,numel(W),numel(W));
	%nn = sum(W(:)); 
	nn = sum(W);
else
	nn = r1;
end; 


if mexFLAG2 && mexFLAG && ~isempty(W),
	%% the mex-functions here are much slower than the m-scripts below 
	%% however, the mex-functions support weighting of samples. 
	if isempty(FLAG_NANS_OCCURED),
		%% mex-files require that FLAG_NANS_OCCURED is not empty, 
		%% otherwise, the status of NAN occurence can not be returned. 
		FLAG_NANS_OCCURED = logical(0);  % default value 
	end;

	if any(Mode=='D') || any(Mode=='E'),
		[S1,N1] = sumskipnan(X,1,W);
		if ~isempty(Y)
	               	[S2,N2] = sumskipnan(Y,1,W);
	        else
	        	S2 = S1; N2 = N1;
		end;
                if any(Mode=='D'), % detrending mode
       			X  = X - ones(r1,1)*(S1./N1);
                        if ~isempty(Y)
                                Y  = Y - ones(r1,1)*(S2./N2);
                        end;
                end;
	end;

	[CC,NN] = covm_mex(real(X), real(Y), FLAG_NANS_OCCURED, W);
	%% complex matrices 
	if ~isreal(X) && ~isreal(Y)
		[iCC,inn] = covm_mex(imag(X), imag(Y), FLAG_NANS_OCCURED, W);
		CC = CC + iCC;
	end; 
	if isempty(Y) Y = X; end;
	if ~isreal(X)
		[iCC,inn] = covm_mex(imag(X), real(Y), FLAG_NANS_OCCURED, W);
		CC = CC - i*iCC;
	end;
	if ~isreal(Y)
		[iCC,inn] = covm_mex(real(X), imag(Y), FLAG_NANS_OCCURED, W);
		CC = CC + i*iCC;
	end;
	
        if any(Mode=='D') && ~any(Mode=='1'),  %  'D1'
                NN = max(NN-1,0);
        end;
        if any(Mode=='E'), % extended mode
                NN = [nn, N2; N1', NN];
                CC = [nn, S2; S1', CC];
        end;


elseif ~isempty(W),

	error('Error COVM: weighted COVM requires sumskipnan_mex and covm_mex but it is not available');

	%% weighted covm without mex-file support
	%% this part is not working.

elseif ~isempty(Y),
	if (~any(Mode=='D') && ~any(Mode=='E')), % if Mode == M
        	NN = real(X==X)'*real(Y==Y);
		FLAG_NANS_OCCURED = any(NN(:)<nn);
	        X(X~=X) = 0; % skip NaN's
	        Y(Y~=Y) = 0; % skip NaN's
        	CC = X'*Y;

        else  % if any(Mode=='D') | any(Mode=='E'), 
        	[S1,N1] = sumskipnan(X,1);
               	[S2,N2] = sumskipnan(Y,1);
               	NN = real(X==X)'*real(Y==Y);

                if any(Mode=='D'), % detrending mode
			X  = X - ones(r1,1)*(S1./N1);
			Y  = Y - ones(r1,1)*(S2./N2);
			if any(Mode=='1'),  %  'D1'
				NN = NN;
			else	%  'D0'       
				NN = max(NN-1,0);
			end;
                end;
		X(X~=X) = 0; % skip NaN's
		Y(Y~=Y) = 0; % skip NaN's
               	CC = X'*Y;

                if any(Mode=='E'), % extended mode
                        NN = [nn, N2; N1', NN];
                        CC = [nn, S2; S1', CC];
                end;
	end;

else
	if (~any(Mode=='D') && ~any(Mode=='E')), % if Mode == M
		tmp = real(X==X);
		NN  = tmp'*tmp;
		X(X~=X) = 0; % skip NaN's
	        CC = X'*X;
		FLAG_NANS_OCCURED = any(NN(:)<nn);

        else  % if any(Mode=='D') | any(Mode=='E'), 
	        [S,N] = sumskipnan(X,1);
       		tmp = real(X==X);
               	NN  = tmp'*tmp;
       	        if any(Mode=='D'), % detrending mode
        	        X  = X - ones(r1,1)*(S./N);
                       	if any(Mode=='1'),  %  'D1'
                               	NN = NN;
                        else  %  'D0'      
       	                        NN = max(NN-1,0);
               	        end;
                end;
                
       	        X(X~=X) = 0; % skip NaN's
		CC = X'*X;
                if any(Mode=='E'), % extended mode
                        NN = [nn, N; N', NN];
                        CC = [nn, S; S', CC];
                end;
	end

end;


if nargout<2
        CC = CC./NN; % unbiased
end;

%!assert(covm([1;NaN;2],'D'),0.5)
%!assert(covm([1;NaN;2],'M'),2.5)
%!assert(covm([1;NaN;2],'E'),[1,1.5;1.5,2.5])
