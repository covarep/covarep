function [FPE,AIC,BIC,SBC,MDL,CATcrit,PHI,optFPE,optAIC,optBIC,optSBC,optMDL,optCAT,optPHI,p,C]=selmo(e,NC);
% Model order selection of an autoregrssive model
% [FPE,AIC,BIC,SBC,MDL,CAT,PHI,optFPE,optAIC,optBIC,optSBC,optMDL,optCAT,optPHI]=selmo(E,N);
%
% E	Error function E(p)
% N	length of the data set, that was used for calculating E(p)
% show  optional; if given the parameters are shown
%
% FPE	Final Prediction Error (Kay 1987, Wei 1990, Priestley 1981  -> Akaike 1969)
% AIC	Akaike Information Criterion (Marple 1987, Wei 1990, Priestley 1981 -> Akaike 1974)
% BIC	Bayesian Akaike Information Criterion (Wei 1990, Priestley 1981 -> Akaike 1978,1979)
% CAT	Parzen's CAT Criterion (Wei 1994 -> Parzen 1974)
% MDL	Minimal Description length Criterion (Marple 1987 -> Rissanen 1978,83)
% SBC	Schwartz's Bayesian Criterion (Wei 1994; Schwartz 1978)
% PHI	Phi criterion (Pukkila et al. 1988, Hannan 1980 -> Hannan & Quinn, 1979)
% HAR	Haring G. (1975)
% JEW	Jenkins and Watts (1968)
%
% optFPE 	order where FPE is minimal
% optAIC 	order where AIC is minimal
% optBIC 	order where BIC is minimal
% optSBC 	order where SBC is minimal
% optMDL 	order where MDL is minimal
% optCAT 	order where CAT is minimal
% optPHI 	order where PHI is minimal
%
% usually is 
% AIC > FPE > *MDL* > PHI > SBC > CAT ~ BIC
%
% REFERENCES:
%  P.J. Brockwell and R.A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  S. Haykin "Adaptive Filter Theory" 3ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.
%  Jenkins G.M. Watts D.G "Spectral Analysis and its applications", Holden-Day, 1968.
%  G. Haring  "Über die Wahl der optimalen Modellordnung bei der Darstellung von stationären Zeitreihen mittels Autoregressivmodell als Basis der Analyse von EEG - Biosignalen mit Hilfe eines Digitalrechners", Habilitationschrift - Technische Universität Graz, Austria, 1975.
%                  (1)"About selecting the optimal model at the representation of stationary time series by means of an autoregressive model as basis of the analysis of EEG - biosignals by means of a digital computer)"
%

%       $Id: selmo.m 9609 2012-02-10 10:18:00Z schloegl $
%       Copyright (C) 1997-2002,2008,2012 by Alois Schloegl <alois.schloegl@ist.ac.at>
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

[lr,lc]=size(e);
if (lr>1) && (lc>1), 
        p=zeros(lr+1,9)+NaN;
else
        p=zeros(1,9)+NaN;
end;

if nargin<2 
        NC=lc*ones(lr,1); 
	NC=(lc-sum(isnan(e)')')*(NC<lc) + NC.*(NC>=lc); % first part 
%end;% Pmax=min([100 N/3]); end;
	%if NC<lc N=lc; end; 
        %NC=(lc-sum(isnan(e)')')*(NC<lc) + NC.*(NC>=lc); % first part 
else
        % NC=NC;
end;

M=lc-1;
m=0:M;

e = e./e(:,ones(1,lc));

for k=0:lr,
        if k>0, % 
                E=e(k,:);
                N=NC(k);
        elseif lr>1
                tmp = e;%(NC>0,:);
                tmp(isnan(tmp)) = 0;
                E = sum(tmp.*(NC*ones(1,lc)))/sum(NC); % weighted average, weigths correspond to number of valid (not missing) values 
                N = sum(NC)./sum(NC>0); % corresponding number of values, 
        else
                E = e;
                N = NC;
        end;
FPE = E.*(N+m)./(N-m);	%OK
		optFPE=find(FPE==min(FPE))-1;	%optimal order
        if isempty(optFPE), optFPE=NaN; end;
AIC = N*log(E)+2*m; 	%OK
	optAIC=find(AIC==min(AIC))-1;	%optimal order
        if isempty(optAIC), optAIC=NaN; end;
AIC4=N*log(E)+4*m;	%OK
	optAIC4=find(AIC4==min(AIC4))-1;	%optimal order
        if isempty(optAIC4), optAIC4=NaN; end;

m=1:M;
BIC=[ N*log(E(1)) N*log(E(m+1)) - (N-m).*log(1-m/N) + m*log(N) + m.*log(((E(1)./E(m+1))-1)./m)];
%BIC=[ N*log(E(1)) N*log(E(m+1)) - m + m*log(N) + m.*log(((E(1)./E(m+1))-1)./m)];
%m=0:M; BIC=N*log(E)+m*log(N);          % Hannan, 1980 -> Akaike, 1977 and Rissanen 1978
        optBIC=find(BIC==min(BIC))-1;	%optimal order
        if isempty(optBIC), optBIC=NaN; end;
        
HAR(2:lc)=-(N-m).*log((N-m).*E(m+1)./(N-m+1)./E(m));         
        HAR(1)=HAR(2);
	optHAR=min(find(HAR<=(min(HAR)+0.2)))-1;	%optimal order
%	optHAR=find(HAR==min(HAR))-1;	%optimal order
        if isempty(optHAR), optHAR=NaN; end;
        
m=0:M;
SBC = N*log(E)+m*log(N);
	optSBC=find(SBC==min(SBC))-1;	%optimal order
        if isempty(optSBC), optSBC=NaN; end;
MDL = N*log(E)+log(N)*m;
	optMDL=find(MDL==min(MDL))-1;	%optimal order
        if isempty(optMDL), optMDL=NaN; end;
        
m=0:M;
%CATcrit= (cumsum(1./E(m+1))/N-1./E(m+1));
E1=N*E./(N-m);
CATcrit= (cumsum(1./E1(m+1))/N-1./E1(m+1));	
	optCAT=find(CATcrit==min(CATcrit))-1;	%optimal order
        if isempty(optCAT), optCAT=NaN; end;

PHI = N*log(E)+2*log(log(N))*m;
	optPHI=find(PHI==min(PHI))-1;	%optimal order
        if isempty(optPHI), optPHI=NaN; end;
        
JEW = E.*(N-m)./(N-2*m-1);	% Jenkins-Watt
	optJEW=find(JEW==min(JEW))-1;	%optimal order
        if isempty(optJEW), optJEW=NaN; end;
        
% in case more than 1 minimum is found, the smaller model order is returned;
p(k+1,:) = [optFPE(1), optAIC(1), optBIC(1), optSBC(1), optCAT(1), optMDL(1), optPHI(1), optJEW(1), optHAR(1)];

end;
C=[FPE;AIC;BIC;SBC;MDL;CATcrit;PHI;JEW;HAR(:)']';
