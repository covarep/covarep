function X = selmo2(y,Pmax); 
% SELMO2 - model order selection for univariate and multivariate 
%   autoregressive models 
% 
%  X = selmo(y,Pmax); 
%  
%  y 	data series
%  Pmax maximum model order 
%  X.A, X.B, X.C parameters of AR model 
%  X.OPT... various optimization criteria
% 
% see also: SELMO, MVAR, 

%	$Id: selmo2.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 2007 by Alois Schloegl <a.schloegl@ieee.org>		
%       This is part of the TSA-toolbox. See also 
%       http://hci.tugraz.at/schloegl/matlab/tsa/
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


[M,N]=size(y); 
if M>N,
	y=y';
end; 
[M,N]=size(y); 

% Univariate AR 
[AutoCov, AutoCorr, ARPMX, E, NC] = invest0(y,Pmax);
[FPE,AIC,BIC,SBC,MDL,CATcrit,PHI,optFPE,optAIC,optBIC,optSBC,optMDL,optCAT,optPHI,s,C] = selmo(E,NC);


% AIC for MVAR model 
[AR,RC,PE] = mvar(y',Pmax); 
E2 = repmat(NaN,1,Pmax); 
for k = 0:Pmax,
	S = PE(:,k*M+(1:M));
	E2(k+1) = trace(S); 

	%%%% FIX ME %%%%% this does not seem right because it depends on the scaling of y
	AIC_MV(k+1) = 2*log(det(S))+2*k*M*M/N; % Ding et al. 2000 refering to Akaike 1974
end; 	


X.A = [eye(M),-AR];
X.B = eye(M);
X.C = S;
X.datatype = 'MVAR'; 
X.MV.AIC = AIC_MV; 
X.UV.FPE = FPE; 
X.UV.AIC = AIC; 
X.UV.BIC = BIC; 
X.UV.MDL = MDL; 
X.UV.CAT = CATcrit; 
X.UV.PHI = PHI; 

X.OPT.FPE = optFPE; 
X.OPT.AIC = optAIC; 
X.OPT.BIC = optBIC; 
X.OPT.MDL = optMDL; 
X.OPT.CAT = optCAT; 
X.OPT.PHI = optPHI; 

[tmp,X.OPT.MVAIC] = min(AIC_MV); 
