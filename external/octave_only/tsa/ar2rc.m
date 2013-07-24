function [MX,res,arg3] = ar2rc(ar);
% converts autoregressive parameters into reflection coefficients 
% with the Durbin-Levinson recursion for multiple channels
% function  [AR,RC,PE] = ar2rc(AR);
% function  [MX,PE] = ar2rc(AR);
%
%  INPUT:
% AR    autoregressive model parameter	
%
%  OUTPUT
% AR    autoregressive model parameter	
% RC    reflection coefficients (= -PARCOR coefficients)
% PE    remaining error variance (relative to PE(1)=1)
% MX    transformation matrix between ARP and RC (Attention: needs O(p^2) memory)
%        AR = MX(:,K*(K-1)/2+(1:K));
%        RC = MX(:,(1:K).*(2:K+1)/2);
%
% All input and output parameters are organized in rows, one row 
% corresponds to the parameters of one channel
%
% see also ACOVF ACORF DURLEV RC2AR 
% 
% REFERENCES:
%  P.J. Brockwell and R. A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  S. Haykin "Adaptive Filter Theory" 3rd ed. Prentice Hall, 1996.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.

%       $Id: ar2rc.m 5090 2008-06-05 08:12:04Z schloegl $
%       Copyright (C) 1998-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>
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

% Inititialization
[lr,lc]=size(ar);
res=[ones(lr,1) zeros(lr,lc)];

if nargout<3         % needs O(p^2) memory 
        MX=zeros(lr,lc*(lc+1)/2);   
	MX(:,lc*(lc-1)/2+(1:lc))=ar;
	
        % Durbin-Levinson Algorithm
        idx=lc*(lc-1)/2;
        for K=lc:-1:2; 
                %idx=K*(K-1)/2;  %see below
                MX(:,(K-2)*(K-1)/2+(1:K-1)) = (MX(:,idx+(1:K-1)) + MX(:,(idx+K)*ones(K-1,1)).*MX(:,idx+(K-1:-1:1)))./((ones(lr,1)-abs(MX(:,idx+K)).^2)*ones(1,K-1));
		idx=idx-K+1;
        end;
	for K=1:lc
		idx=K*(K-1)/2;  %see below
                res(:,K+1) = res(:,K).*(1-abs(MX(:,idx+K)).^2);
        end;
	    
        %arp=MX(:,K*(K-1)/2+(1:K));
        %rc =MX(:,(1:K).*(2:K+1)/2);

else            % needs O(p) memory 
        
        %ar=zeros(lr,lc);
        rc=zeros(lr,lc);
	rc(:,lc)=ar(:,lc);
	MX=ar; % assign output
	
	% Durbin-Levinson Algorithm
        for K=lc-1:-1:1,
                ar(:,1:K)=(ar(:,1:K)+ar(:,(K+1)*ones(K,1)).*ar(:,K:-1:1))./((ones(lr,1)-abs(ar(:,K+1)).^2)*ones(1,K));
                rc(:,K)=ar(:,K);
	end;
	
	for K=1:lc,
                res(:,K+1) = res(:,K) .* (1-abs(ar(:,K)).^2);
        end;
        
        % assign output arguments
        arg3=res;
        res=rc;
        %MX=ar;
end; %if
