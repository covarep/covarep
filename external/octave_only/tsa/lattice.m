 function [MX,PE,arg3] = lattice(Y,lc,Mode);
% Estimates AR(p) model parameter with lattice algorithm (Burg 1968) 
% for multiple channels. 
% If you have the NaN-tools, LATTICE.M can handle missing values (NaN), 
%
% [...] = lattice(y [,Pmax [,Mode]]);
%
% [AR,RC,PE] = lattice(...);
% [MX,PE] = lattice(...);
%
%  INPUT:
% y	signal (one per row), can contain missing values (encoded as NaN)
% Pmax	max. model order (default size(y,2)-1))
% Mode  'BURG' (default) Burg algorithm
%	'GEOL' geometric lattice
%
%  OUTPUT
% AR    autoregressive model parameter	
% RC    reflection coefficients (= -PARCOR coefficients)
% PE    remaining error variance
% MX    transformation matrix between ARP and RC (Attention: needs O(p^2) memory)
%        AR(:,K) = MX(:, K*(K-1)/2+(1:K)); = MX(:,sum(1:K-1)+(1:K)); 
%        RC(:,K) = MX(:,cumsum(1:K));      = MX(:,(1:K).*(2:K+1)/2);
%
% All input and output parameters are organized in rows, one row 
% corresponds to the parameters of one channel
%
% see also ACOVF ACORF AR2RC RC2AR DURLEV SUMSKIPNAN 
% 
% REFERENCE(S):
%  J.P. Burg, "Maximum Entropy Spectral Analysis" Proc. 37th Meeting of the Society of Exp. Geophysiscists, Oklahoma City, OK 1967
%  J.P. Burg, "Maximum Entropy Spectral Analysis" PhD-thesis, Dept. of Geophysics, Stanford University, Stanford, CA. 1975.
%  P.J. Brockwell and R. A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  S.   Haykin "Adaptive Filter Theory" 3rd ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.

%	$Id: lattice.m 7687 2010-09-08 18:39:23Z schloegl $ 
%	Copyright (C) 1996-2002,2008,2010 by Alois Schloegl <a.schloegl@ieee.org>
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


if nargin<3, Mode='BURG'; 
else Mode=upper(Mode(1:4));end;
BURG=~strcmp(Mode,'GEOL');

% Inititialization
[lr,N]=size(Y);
if nargin<2, lc=N-1; end;
F=Y;
B=Y;
[DEN,nn] = sumskipnan((Y.*Y),2);
PE = [DEN./nn,zeros(lr,lc)];

if nargout<3         % needs O(p^2) memory 
        MX = zeros(lr,lc*(lc+1)/2);   
        idx= 0;
        
        % Durbin-Levinson Algorithm
        for K=1:lc,
                [TMP,nn] = sumskipnan(F(:,K+1:N).*B(:,1:N-K),2);
                MX(:,idx+K) = TMP./DEN; %Burg
                if K>1,   %for compatibility with OCTAVE 2.0.13
                        MX(:,idx+(1:K-1))=MX(:,(K-2)*(K-1)/2+(1:K-1))-MX(:,(idx+K)*ones(K-1,1)).*MX(:,(K-2)*(K-1)/2+(K-1:-1:1));
                end;   
                
                tmp = F(:,K+1:N) - MX(:,(idx+K)*ones(1,N-K)).*B(:,1:N-K);
                B(:,1:N-K) = B(:,1:N-K) - MX(:,(idx+K)*ones(1,N-K)).*F(:,K+1:N);
                F(:,K+1:N) = tmp;
                
                [PE(:,K+1),nn] = sumskipnan([F(:,K+1:N).^2,B(:,1:N-K).^2],2);        
                if ~BURG,
                        [f,nf] = sumskipnan(F(:,K+1:N).^2,2);
                        [b,nb] = sumskipnan(B(:,1:N-K).^2,2); 
                        DEN = sqrt(b.*f); 
                else
                        DEN = PE(:,K+1);
                end;
                idx=idx+K;
		PE(:,K+1) = PE(:,K+1)./nn; 	% estimate of covariance
        end;
else            % needs O(p) memory 
        arp=zeros(lr,lc-1);
        rc=zeros(lr,lc-1);
        % Durbin-Levinson Algorithm
        for K=1:lc,
                [TMP,nn] = sumskipnan(F(:,K+1:N).*B(:,1:N-K),2);
                arp(:,K) = TMP./DEN; %Burg
                rc(:,K)  = arp(:,K);
                if K>1,	% for compatibility with OCTAVE 2.0.13
                        arp(:,1:K-1) = arp(:,1:K-1) - arp(:,K*ones(K-1,1)).*arp(:,K-1:-1:1);
                end;
                
                tmp = F(:,K+1:N) - rc(:,K*ones(1,N-K)).*B(:,1:N-K);
                B(:,1:N-K) = B(:,1:N-K) - rc(:,K*ones(1,N-K)).*F(:,K+1:N);
                F(:,K+1:N) = tmp;
                
                [PE(:,K+1),nn] = sumskipnan([F(:,K+1:N).^2,B(:,1:N-K).^2],2);        
                if ~BURG,
                        [f,nf] = sumskipnan(F(:,K+1:N).^2,2);
                        [b,nb] = sumskipnan(B(:,1:N-K).^2,2); 
                        DEN = sqrt(b.*f); 
                else
                        DEN = PE(:,K+1);
                end;
		PE(:,K+1) = PE(:,K+1)./nn; 	% estimate of covariance
        end;
% assign output arguments
	arg3=PE;
        PE=rc;
        MX=arp;
end; %if

