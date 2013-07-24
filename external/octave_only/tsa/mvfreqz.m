function [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF, pCOH2, PDCF, coh,GGC,Af,GPDC,GGC2]=mvfreqz(B,A,C,N,Fs)
% MVFREQZ multivariate frequency response
% [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF,pCOH2,PDCF,coh,GGC,Af,GPDC] = mvfreqz(B,A,C,f,Fs)
% [...]  = mvfreqz(B,A,C,N,Fs)
%  
% INPUT: 
% ======= 
% A, B	multivariate polynomials defining the transfer function
%
%    a0*Y(n) = b0*X(n) + b1*X(n-1) + ... + bq*X(n-q)
%                          - a1*Y(n-1) - ... - ap*Y(:,n-p)
%
%  A=[a0,a1,a2,...,ap] and B=[b0,b1,b2,...,bq] must be matrices of
%  size  Mx((p+1)*M) and Mx((q+1)*M), respectively. 
%
%  C is the covariance of the input noise X (i.e. D'*D if D is the mixing matrix)
%  N if scalar, N is the number of frequencies 
%    if N is a vector, N are the designated frequencies. 
%  Fs sampling rate [default 2*pi]
% 
%  A,B,C and D can by obtained from a multivariate time series 
%       through the following commands: 
%  [AR,RC,PE] = mvar(Y,P);
%       M = size(AR,1); % number of channels       
%       A = [eye(M),-AR];
%       B = eye(M); 
%       C = PE(:,M*P+1:M*(P+1)); 
%
% Fs 	sampling rate in [Hz]
% (N 	number of frequencies for computing the spectrum, this will become OBSOLETE),  
% f	vector of frequencies (in [Hz])  
%
%
% OUTPUT: 
% ======= 
% S   	power spectrum
% h	transfer functions, abs(h.^2) is the non-normalized DTF [11]
% PDC 	partial directed coherence [2]
% DC  	directed coupling	
% COH 	coherency (complex coherence) [5]
% DTF 	directed transfer function
% pCOH 	partial coherence
% dDTF 	direct Directed Transfer function
% ffDTF full frequency Directed Transfer Function 
% pCOH2 partial coherence - alternative method 
% GGC	a modified version of Geweke's Granger Causality [Geweke 1982]
%	   !!! it uses a Multivariate AR model, and computes the bivariate GGC as in [Bressler et al 2007]. 
%	   This is not the same as using bivariate AR models and GGC as in [Bressler et al 2007]
% Af	Frequency transform of A(z), abs(Af.^2) is the non-normalized PDC [11]
% PDCF 	Partial Directed Coherence Factor [2]
% GPDC 	Generalized Partial Directed Coherence [9,10]
%
% see also: FREQZ, MVFILTER, MVAR
% 
% REFERENCE(S):
% [1] H. Liang et al. Neurocomputing, 32-33, pp.891-896, 2000. 
% [2] L.A. Baccala and K. Samashima, Biol. Cybern. 84,463-474, 2001. 
% [3] A. Korzeniewska, et al. Journal of Neuroscience Methods, 125, 195-207, 2003. 
% [4] Piotr J. Franaszczuk, Ph.D. and Gregory K. Bergey, M.D.
% 	Fast Algorithm for Computation of Partial Coherences From Vector Autoregressive Model Coefficients
%	World Congress 2000, Chicago. 
% [5] Nolte G, Bai O, Wheaton L, Mari Z, Vorbach S, Hallett M.
%	Identifying true brain interaction from EEG data using the imaginary part of coherency.
%	Clin Neurophysiol. 2004 Oct;115(10):2292-307. 
% [6] Schlogl A., Supp G.
%       Analyzing event-related EEG data with multivariate autoregressive parameters.
%       (Eds.) C. Neuper and W. Klimesch, 
%       Progress in Brain Research: Event-related Dynamics of Brain Oscillations. 
%       Analysis of dynamics of brain oscillations: methodological advances. Elsevier. 
%       Progress in Brain Research 159, 2006, p. 135 - 147
% [7] Bressler S.L., Richter C.G., Chen Y., Ding M. (2007)
%	Cortical fuctional network organization from autoregressive modelling of loal field potential oscillations.
%	Statistics in Medicine, doi: 10.1002/sim.2935 
% [8] Geweke J., 1982	
%	J.Am.Stat.Assoc., 77, 304-313.
% [9] L.A. Baccala, D.Y. Takahashi, K. Sameshima. (2006) 
% 	Generalized Partial Directed Coherence. 
%	Submitted to XVI Congresso Brasileiro de Automatica, Salvador, Bahia.  
% [10] L.A. Baccala, D.Y. Takahashi, K. Sameshima. 
% 	Computer Intensive Testing for the Influence Between Time Series, 
%	Eds. B. Schelter, M. Winterhalder, J. Timmer: 
%	Handbook of Time Series Analysis - Recent Theoretical Developments and Applications
%	Wiley, p.413, 2006.
% [11] M. Eichler
%	On the evaluation of informatino flow in multivariate systems by the directed transfer function
%	Biol. Cybern. 94: 469-482, 2006. 	

%	$Id: mvfreqz.m 8141 2011-03-02 08:01:58Z schloegl $
%	Copyright (C) 1996-2008 by Alois Schloegl <a.schloegl@ieee.org>	
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

[K1,K2] = size(A);
p = K2/K1-1;
%a=ones(1,p+1);
[K1,K2] = size(B);
q = K2/K1-1;
%b=ones(1,q+1);
if nargin<3
        C = eye(K1,K1);
end;
if nargin<5,
        Fs= 1;        
end;
if nargin<4,
        N = 512;
        f = (0:N-1)*(Fs/(2*N));
end;
if all(size(N)==1),	
	fprintf(1,'Warning MVFREQZ: The forth input argument N is a scalar, this is ambigous.\n');
	fprintf(1,'   In the past, N was used to indicate the number of spectral lines. This might change.\n');
	fprintf(1,'   In future versions, it will indicate the spectral line.\n');
        f = (0:N-1)*(Fs/(2*N));
else
        f = N;
end;
N = length(f);
s = exp(i*2*pi*f/Fs);
z = i*2*pi/Fs;

h=zeros(K1,K1,N);
Af=zeros(K1,K1,N);
g=zeros(K1,K1,N);
S=zeros(K1,K1,N);
S1=zeros(K1,K1,N);
DTF=zeros(K1,K1,N);
COH=zeros(K1,K1,N);
%COH2=zeros(K1,K1,N);
PDC=zeros(K1,K1,N);
%PDC3=zeros(K1,K1,N);
PDCF = zeros(K1,K1,N);
pCOH = zeros(K1,K1,N);
GGC=zeros(K1,K1,N);
GGC2=zeros(K1,K1,N);
invC=inv(C);
tmp1=zeros(1,K1);
tmp2=zeros(1,K1);

M = zeros(K1,K1,N);
detG = zeros(N,1);

%D = sqrtm(C);
%iD= inv(D);
ddc2 = diag(diag(C).^(-1/2)); 
for n=1:N,
        atmp = zeros(K1);
        for k = 1:p+1,
                atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;   
        
        % compensation of instantaneous correlation 
        % atmp = iD*atmp*D;
        
        btmp = zeros(K1);
        for k = 1:q+1,
                btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;        
        h(:,:,n)  = atmp\btmp;        
        Af(:,:,n)  = atmp/btmp;        
        S(:,:,n)  = h(:,:,n)*C*h(:,:,n)'/Fs;        
        S1(:,:,n) = h(:,:,n)*h(:,:,n)';        
        ctmp = ddc2*atmp;	%% used for GPDC 
        for k1 = 1:K1,
                tmp = squeeze(atmp(:,k1));
                tmp1(k1) = sqrt(tmp'*tmp);
                tmp2(k1) = sqrt(tmp'*invC*tmp);

                %tmp = squeeze(atmp(k1,:)');
                %tmp3(k1) = sqrt(tmp'*tmp);

                tmp = squeeze(ctmp(:,k1));
                tmp3(k1) = sqrt(tmp'*tmp);
        end;
        
        PDCF(:,:,n)  = abs(atmp)./tmp2(ones(1,K1),:);
        PDC(:,:,n)   = abs(atmp)./tmp1(ones(1,K1),:);
        GPDC(:,:,n)  = abs(ctmp)./tmp3(ones(1,K1),:);
        %PDC3(:,:,n) = abs(atmp)./tmp3(:,ones(1,K1));
        
        g = atmp/btmp;        
        G(:,:,n) = g'*invC*g;
        detG(n) = det(G(:,:,n));        
end;

if nargout<4, return; end;

%%%%% directed transfer function
for k1=1:K1;
        DEN=sum(abs(h(k1,:,:)).^2,2);	        
        for k2=1:K2;
                %COH2(k1,k2,:) = abs(S(k1,k2,:).^2)./(abs(S(k1,k1,:).*S(k2,k2,:)));
                COH(k1,k2,:)   = (S(k1,k2,:))./sqrt(abs(S(k1,k1,:).*S(k2,k2,:)));
                coh(k1,k2,:)   = (S1(k1,k2,:))./sqrt(abs(S1(k1,k1,:).*S1(k2,k2,:)));
                %DTF(k1,k2,:)  = sqrt(abs(h(k1,k2,:).^2))./DEN;	        
                DTF(k1,k2,:)   = abs(h(k1,k2,:))./sqrt(DEN);
                ffDTF(k1,k2,:) = abs(h(k1,k2,:))./sqrt(sum(DEN,3));
                pCOH2(k1,k2,:) = abs(G(k1,k2,:).^2)./(G(k1,k1,:).*G(k2,k2,:));
                
                %M(k2,k1,:) = ((-1)^(k1+k2))*squeeze(G(k1,k2,:))./detG; % oder ist M = G?
        end;
end;

dDTF = pCOH2.*ffDTF; 

if nargout<6, return; end;

DC = zeros(K1);
for k = 1:p,
        DC = DC + A(:,k*K1+(1:K1)).^2;
end;        


if nargout<13, return; end;

for k1=1:K1;
        for k2=1:K2;
		% Bivariate Granger Causality (similar to Bressler et al. 2007. )
                GGC(k1,k2,:) = ((C(k1,k1)*C(k2,k2)-C(k1,k2)^2)/C(k2,k2))*real(h(k1,k2,:).*conj(h(k1,k2,:)))./abs(S(k2,k2,:));
                %GGC2(k1,k2,:) = -log(1-((C(k1,k1)*C(k2,k2)-C(k1,k2)^2)/C(k2,k2))*real(h(k1,k2,:).*conj(h(k1,k2,:)))./S(k2,k2,:));
        end;
end;

return;

if nargout<7, return; end;

for k1=1:K1;
        for k2=1:K2;
                M(k2,k1,:) = ((-1)^(k1+k2))*squeeze(G(k1,k2,:))./detG; % oder ist M = G?
        end;
end;


for k1=1:K1;
        for k2=1:K2;
                pCOH(k1,k2,:) = abs(M(k1,k2,:).^2)./(M(k1,k1,:).*M(k2,k2,:));
        end;
end;

