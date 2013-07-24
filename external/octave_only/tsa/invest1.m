function [AutoCov,AutoCorr,ARPMX,E,C,s]=invest1(Y,Pmax,D);
% First Investigation of a signal (time series) - interactive
% [AutoCov,AutoCorr,ARPMX,E,CRITERIA,MOPS]=invest1(Y,Pmax,show);
%
% Y	time series
% Pmax	maximal order (optional)
% show  optional; if given the parameters are shown
%
% AutoCov	Autocorrelation 
% AutoCorr	normalized Autocorrelation
% PartACF	Partial Autocorrelation
% E	Error function E(p)
% CRITERIA curves of the various (see below) criteria, 
% MOPS=[optFPE optAIC optBIC optSBC optMDL optCAT optPHI];
%      optimal model order according to various criteria
%
% FPE	Final Prediction Error (Kay, 1987)
% AIC	Akaike Information Criterion (Marple, 1987)
% BIC	Bayesian Akaike Information Criterion (Wei, 1994)
% SBC	Schwartz's Bayesian Criterion (Wei, 1994)
% MDL	Minimal Description length Criterion (Marple, 1987)
% CAT	Parzen's CAT Criterion (Wei, 1994)
% PHI	Phi criterion (Pukkila et al. 1988)
% minE		order where E is minimal
%
% REFERENCES:
%  P.J. Brockwell and R.A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  S.   Haykin "Adaptive Filter Theory" 3ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.

% optFPE 	order where FPE is minimal
% optAIC 	order where AIC is minimal
% optBIC 	order where BIC is minimal
% optSBC 	order where SBC is minimal
% optMDL 	order where MDL is minimal
% optCAT 	order where CAT is minimal
% optPHI 	order where PHI is minimal
% optRC2        max reflection coefficient larger than std-error

%	$Id: invest1.m 7687 2010-09-08 18:39:23Z schloegl $
%	Copyright (C) 1998-2002,2008,2010 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the TSA-toolbox. See also 
%       http://biosig-consulting.com/matlab/tsa/
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

N=length(Y); 
[nr,nc]=size(Y);
if nc==1 Y=transpose(Y); nc=nr; nr=1; end;
    
if nargin<2 Pmax=min([100 nc/3]); end;

if exist('OCTAVE_VERSION'),
        fprintf(2,'Warning INVEST1: DIFF-based optimization not possible\n');
        %%% missing DIM-argument in DIFF.M    
else
        %tmp=Y-mean(Y,2)*ones(1,nc);
        RMS(:,1) = mean(Y.^2,2);
        Dmax = min(Pmax,5);
        for k = 1:Dmax,
                RMS(:,k+1) = mean(diff(Y,k,2).^2,2);
        end;
        [tmp, orderDIFF] = min(RMS,[],2);
        
        % show a nice histogram
        h = histo3(orderDIFF-1);
        X = 0:Dmax; H = zeros(1,Dmax+1); for k=1:length(h.X), H(find(X==h.X(k)))=h.H(k); end;
        %X = 0:Dmax; H = zeros(1,Dmax+1); for k=1:length(x), H(find(X==x(k)))=h(k); end;
        bar(X,H);
        drawnow;

	if nargin>2
		oD=0;
	else	        
    		oD=input('Which order should be used for differentiating [default=0] ?: ');
        end;
	if oD>0
                Y=diff(Y,oD,2);
        end;
end;

[AutoCov, AutoCorr, ARPMX, E, NC] = invest0(Y,Pmax);

[FPE,AIC,BIC,SBC,MDL,CATcrit,PHI,optFPE,optAIC,optBIC,optSBC,optMDL,optCAT,optPHI,s,C] = selmo(E,NC);

if 0,
optRC2=zeros(nr+1,1);
for k=0:nr,
        if k>0
                optRC2(k+1)=max(find(abs(ARPMX(k,(1:Pmax).*(2:Pmax+1)/2))*sqrt(size(Y,2))>1));
        else 
                optRC2(k+1)=max(find(mean(abs(ARPMX(:,(1:Pmax).*(2:Pmax+1)/2))*sqrt(size(Y,2)),2)>1));
        end;
end;
%GERSCH=min(find(rc.^2<(0.05/1.05)));
s=[s optRC2];
end;

%CRITERIA=([FPE;AIC;BIC;SBC;MDL;CATcrit;PHI])';
MOPS = s(1:size(s,1),:); %[optFPE optAIC optBIC optSBC optMDL optCAT optPHI];

if nargin==3,  
        if size(ARPMX,2)==2*Pmax,
                %invest1(eeg8s,30,'s');
                AR=ARPMX(:,1:Pmax);
                RC=ARPMX(:,Pmax+1:2*Pmax);
        else
                AR=ARPMX(:,Pmax/2*(Pmax-1)+(1:Pmax));
                RC=ARPMX(:,(1:Pmax).*(2:Pmax+1)/2);
        end;
        oo=optBIC; 
        sinvest1;
end;
