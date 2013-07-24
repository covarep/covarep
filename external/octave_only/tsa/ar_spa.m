function  [w,A,B,R,P,F,ip] = ar_spa(ARP,nhz,E);
% AR_SPA decomposes an AR-spectrum into its compontents 
% [w,A,B,R,P,F,ip] = ar_spa(AR,fs,E);
%
%  INPUT:
% AR   autoregressive parameters
% fs    sampling rate, provide w and B in [Hz], if not given the result is in radians 
% E     noise level (mean square),  gives A and F in units of E, if not given as relative amplitude
%
%  OUTPUT
% w	center frequency
% A     Amplitude
% B     bandwidth
%       - less important output parameters - 
% R	residual
% P	poles
% ip	number of complex conjugate poles
% real(F)     	power, absolute values are obtained by multiplying with noise variance E(p+1) 
% imag(F)	assymetry, - " -
%
% All input and output parameters are organized in rows, one row 
% corresponds to the parameters of one channel
%
% see also ACOVF ACORF DURLEV IDURLEV PARCOR YUWA 
% 
% REFERENCES:
% [1] Zetterberg L.H. (1969) Estimation of parameter for linear difference equation with application to EEG analysis. Math. Biosci., 5, 227-275. 
% [2] Isaksson A. and Wennberg, A. (1975) Visual evaluation and computer analysis of the EEG - A comparison. Electroenceph. clin. Neurophysiol., 38: 79-86.
% [3] G. Florian and G. Pfurtscheller (1994) Autoregressive model based spectral analysis with application to EEG. IIG - Report Series, University of Technolgy Graz, Austria.

% 	$Id: ar_spa.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2003 by Alois Schloegl <a.schloegl@ieee.org>
%	This is part of the TSA-toolbox see also: 
% 	   http://hci.tugraz.at/schloegl/matlab/tsa/
%	   http://octave.sf.net/
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


[NTR,pp]=size(ARP);

R=zeros(size(ARP));
P=zeros(size(ARP));
w=zeros(size(ARP));
A=zeros(size(ARP));
B=zeros(size(ARP));
F=zeros(size(ARP));
F1 = F;
for k = 1:NTR, %if ~mod(k,100),k, end;
	[r,p,tmp] = residue(1,[1 -ARP(k,:)]);
	[tmp,idx] = sort(-abs(r));   
	R(k,:) = r(idx)';		% Residual, 
   	P(k,:) = p(idx)';		% Poles
   	%r(k,:)=roots([1 -ARP(k,:)])';
   	w(k,:) = angle(p(idx)');	% center frequency (in [radians])
   	A(k,:) = 1./abs(polyval([1 -ARP(k,:)],exp(i*w(k,:))));	% Amplitude 
   	%A(k,:) = freqz(1,[1 -ARP(k,:)],w(k,:));	% Amplitude 
   	%A2(k,:) = abs(r)'./abs(exp(i*w(k,:))-r');   % Amplitude
   	B(k,:) = -log(abs(p(idx)'));  % Bandwidth
           
        if nargout < 6,  

	elseif 0,	
	        F(k,:) = (1+sign(imag(r(idx)')))./(polyval([-ARP(k,pp-1:-1:1).*(1:pp-1) pp],1./p(idx).').*polyval([-ARP(k,pp:-1:1) 1],p(idx).'));        

        elseif 1;
	        a3 = polyval([-ARP(k,pp-1:-1:1).*(1:pp-1), pp],1./p(idx).');
	        a  = polyval([-ARP(k,pp:-1:1) 1],p(idx).');
		%F(k,:) = (1+(imag(P(k,:))~=0))./(a.*a3); 
		F(k,:) = (1+sign(imag(P(k,:))))./(a.*a3); 
        end;	
end;

A = A.*sqrt(E(:,ones(1,pp))/(2*pi*nhz));
if nargin>1,
        if size(nhz,1)==1,
                nhz = nhz(ones(NTR,1),:);
        end;
        w = w.*nhz(:,ones(1,pp))/(2*pi);
        B = B.*nhz(:,ones(1,pp))/(2*pi);
end;
if nargin>2,
        F = F.*E(:,ones(1,pp));
        F1 = F1.*E(:,ones(1,pp));
end;

ip = sum(imag(P)~=0,2)/2; 
return;

np(:,1) = sum(imag(P')==0)';	% number of real poles
np(:,2) = pp-np(:,1);		% number of imaginary poles


