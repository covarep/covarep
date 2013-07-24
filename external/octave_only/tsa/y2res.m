function [R]=y2res(Y)
% Y2RES evaluates basic statistics of a data series
% 
% R = y2res(y)
%	several statistics are estimated from each column of y
% 
% OUTPUT:
%   R.N     number of samples, NaNs are not counted 
%   R.SUM   sum of samples
%   R.MEAN  mean
%   R.STD   standard deviation 
%   R.VAR   variance
%   R.Max   Maximum
%   R.Min   Minimum 
%   ...   and many more including:  
%	MEDIAN, Quartiles, Variance, standard error of the mean (SEM), 
%	Coefficient of Variation, Quantization (QUANT), TRIMEAN, SKEWNESS, 
%	KURTOSIS, Root-Mean-Square (RMS), ENTROPY 
% 

%	$Id: y2res.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2005,2008 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the TSA-toolbox 
%       http://octave.sourceforge.net/
%	http://www.dpmi.tugraz.at/~schloegl/matlab/tsa/
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


[R.SUM, R.N, R.SSQ] = sumskipnan(Y,1);
%R.S3P = sumskipnan(Y.^3,1);
R.S4P = sumskipnan(Y.^4,1);
%R.S5P = sumskipnan(Y.^5,1);

R.MEAN	= R.SUM./R.N;
R.MSQ   = R.SSQ./R.N;
R.RMS	= sqrt(R.MSQ);
R.SSQ0  = R.SSQ-R.SUM.*R.MEAN;		% sum square of mean removed

if 1,%flag_implicit_unbiased_estim,
    n1 	= max(R.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and STE are INF
else
    n1	= R.N;
end;

R.VAR  	= R.SSQ0./n1;	     		% variance (unbiased) 
R.STD  	= sqrt(R.VAR);		     	% standard deviation
R.SEM  	= sqrt(R.SSQ0./(R.N.*n1)); 	% standard error of the mean
R.SEV	= sqrt(n1.*(n1.*R.S4P./R.N+(R.N.^2-2*R.N+3).*(R.SSQ./R.N).^2)./(R.N.^3)); % standard error of the variance
R.Coefficient_of_variation = R.STD./R.MEAN;

R.CM2	= R.SSQ0./n1;

R.Max   = max(Y,[],1);
R.Min   = min(Y,[],1);

%R.NormEntropy = log2(sqrt(2*pi*exp(1)))+log2(R.STD);

Q0500=repmat(nan,1,size(Y,2));
Q0250=Q0500;
Q0750=Q0500;
%MODE=Q0500;
for k = 1:size(Y,2),
        tmp = sort(Y(:,k));
        Q0250(k) = flix(tmp,R.N(k)/4   + 0.75);
        Q0500(k) = flix(tmp,R.N(k)/2   + 0.50);
        Q0750(k) = flix(tmp,R.N(k)*3/4 + 0.25);
        tmp = diff(tmp);

	pdf   = diff([0; find(tmp>0); R.N(k)])/R.N(k); % empirical probability distribution 
	R.ENTROPY(k) = -sumskipnan(pdf.*log(pdf));

        tmp = tmp(find(tmp));
        q   = min(tmp);
	qerror = 0; 
        if isempty(q),
                q = NaN;
        else
                tmp = tmp/q; 
		qerror = max(abs(tmp-round(tmp)));
        end;
        R.QUANT(k) = q; 
	R.Qerror(k) = qerror; 
end;
if any(R.Qerror*1e6>R.QUANT)
	warning('(Y2RES) Quantization might not be equidistant')
end;	

R.MEDIAN 	= Q0500;
R.Quartiles   	= [Q0250; Q0750];
% R.IQR = H_spread   	= [Q0750 - Q0250];
R.TRIMEAN	= [Q0250 + 2*Q0500 + Q0750]/4;

Y       = Y - repmat(R.MEAN,size(Y)./size(R.MEAN));
R.CM3 	= sumskipnan(Y.^3,1)./n1;
R.CM4 	= sumskipnan(Y.^4,1)./n1;
%R.CM5 	= sumskipnan(Y.^5,1)./n1;

R.SKEWNESS = R.CM3./(R.STD.^3);
R.KURTOSIS = R.CM4./(R.VAR.^2)-3;

%R.Skewness.Fisher = (R.CM3)./(R.STD.^3);	%%% same as R.SKEWNESS

%R.Skewness.Pearson_Mode   = (R.MEAN-R.MODE)./R.STD;
%R.Skewness.Pearson_coeff1 = (3*R.MEAN-R.MODE)./R.STD;
R.Skewness.Pearson_coeff2 = (3*R.MEAN-R.MEDIAN)./R.STD;
R.Skewness.Bowley = (Q0750+Q0250 - 2*Q0500)./(Q0750-Q0250); % quartile skewness coefficient

R.datatype = 'STAT Level 4';

