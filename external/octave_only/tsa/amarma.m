function [z,e,REV,ESU,V,Z,SPUR] = amarma(y, Mode, MOP, UC, z0, Z0, V0, W); 
% Adaptive Mean-AutoRegressive-Moving-Average model estimation
% [z,E,ESU,REV,V,Z,SPUR] = amarma(y, mode, MOP, UC, z0, Z0, V0, W); 
% Estimates AAR parameters with Kalman filter algorithm
% 	y(t) = sum_i(a(i,t)*y(t-i)) + mu(t) + E(t)
%
% State space model:
%	z(t)=G*z(t-1) + w(t)      w(t)=N(0,W) 
%	y(t)=H*z(t)   + v(t)	  v(t)=N(0,V)	
%
% G = I, 
% z = [µ(t)/(1-sum_i(a(i,t))),a_1(t-1),..,a_p(t-p),b_1(t-1),...,b_q(t-q)];
% H = [1,y(t-1),..,y(t-p),e(t-1),...,e(t-q)];
% W = E{(z(t)-G*z(t-1))*(z(t)-G*z(t-1))'}
% V = E{(y(t)-H*z(t-1))*(y(t)-H*z(t-1))'}
%
% Input:
%       y	Signal (AR-Process)
%       Mode
%	    [0,0] uses V0 and W  
%
%       MOP     Model order [m,p,q], default [0,10,0] 
%			m=1 includes the mean term, m=0 does not. 
%			p and q must be positive integers
%			it is recommended to set q=0. 
%	UC	Update Coefficient, default 0
%	z0	Initial state vector
%	Z0	Initial Covariance matrix
%      
% Output:
%	z	AR-Parameter
%	E	error process (Adaptively filtered process)
%       REV     relative error variance MSE/MSY
%
%
% see also: AAR
%
% REFERENCE(S): 
% [1] A. Schloegl (2000), The electroencephalogram and the adaptive autoregressive model: theory and applications. 
%     ISBN 3-8265-7640-3 Shaker Verlag, Aachen, Germany. 
% [2] Schlögl A, Lee FY, Bischof H, Pfurtscheller G
%     Characterization of Four-Class Motor Imagery EEG Data for the BCI-Competition 2005.
%     Journal of neural engineering 2 (2005) 4, S. L14-L22
%
% More references can be found at 
%     http://www.dpmi.tu-graz.ac.at/~schloegl/publications/

%	$Id: amarma.m 5376 2008-10-13 15:53:47Z schloegl $
%       Copyright (C) 1998-2002,2005,2006,2007,2008 by  Alois Schloegl <a.schloegl@ieee.org>
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


[nc,nr]=size(y);

if nargin<2 Mode=0; 
elseif isnan(Mode) return; end;
if nargin<3, MOP=[0,10,0]; end;
if nargin<8, W  = nan ; end;
if length(MOP)==0,      m=0;p=10; q=0; MOP=p;
elseif length(MOP)==1,  m=0;p=MOP(1); q=0; MOP=p;
elseif length(MOP)==2,  fprintf(1,'Error AMARMA: MOP is ambiguos\n');
elseif length(MOP)>2,   m=MOP(1); p=MOP(2); q=MOP(3);MOP=m+p+q;
end;
       
if prod(size(Mode))>1
        aMode=Mode(1);
        eMode=Mode(2);
end;
%fprintf(1,['a' int2str(aMode) 'e' int2str(eMode) ' ']);

       
e = zeros(nc,1);
V = zeros(nc,1);V(1)=V0;
T = zeros(nc,1);
ESU = zeros(nc,1)+nan;
SPUR = zeros(nc,1)+nan;
z = z0(ones(nc,1),:);

dW = UC/MOP*eye(MOP);                % Schloegl

%------------------------------------------------
%	First Iteration
%------------------------------------------------

H = zeros(MOP,1); 
if m, 
    %M0   = z0(1)/(1-sum(z0(2:p+1))); %transformierter Mittelwert
    H(1) = 1;%M0; 
    %z0(1)= 1;
end; 

Z = Z0;
zt= z0;

A1 = zeros(MOP); A2 = A1;

%------------------------------------------------
%	Update Equations
%------------------------------------------------
        
for t=1:nc,
        %H=[y(t-1); H(1:p-1); E ; H(p+1:MOP-1)]
        
 
        if t<=p, H(m+(1:t-1)) = y(t-1:-1:1);     %H(p)=mu0;          % Autoregressive 
        else     H(m+(1:p)) = y(t-1:-1:t-p); %mu0]; 
	end;
        
        if t<=q, H(m+p+(1:t-1)) = e(t-1:-1:1);       % Moving Average
        else     H(m+p+(1:q)) = e(t-1:-1:t-q); 
	end;
        
        % Prediction Error 
        E = y(t) - zt*H;
        
        e(t) = E;
        
        if ~isnan(E),
                E2 = E*E;
	        AY = Z*H; 
                
%                [zt, t, y(t), E,ESU(t),V(t),H,Z],pause,
                
                ESU(t) = H'*AY;
  
	        if eMode==0
	    		V(t) = V0;        
	        elseif eMode==1
			V0 = V(t-1);
	    		V(t) = V0*(1-UC)+UC*E2;        
	        elseif eMode==2
			V0 = 1;
	        	V(t) = V0; %V(t-1)*(1-UC)+UC*E2;        
	        elseif eMode==3
			V0 = 1-UC;
	        	V(t) = V0; %(t-1)*(1-UC)+UC*E2;        
	        elseif eMode==4
		        V0 = V0*(1-UC)+UC*E2;        
			V(t) = V0;
	        elseif eMode==5
			V(t)=V0;
			%V0 = V0;
	        elseif eMode==6
        	        if E2>ESU(t) 
                	        V0=(1-UC)*V0+UC*(E2-ESU(t));
                        end;
                        V(t)=V0;
	        elseif eMode==7
                        V0=V(t); 
        	        if E2>ESU(t) 
                	        V(t) = (1-UC)*V0+UC*(E2-ESU(t));
	                else 
        	                V(t) = V0;
                        end;
	        elseif eMode==8
			V0=0;
	        	V(t) = V0; % (t-1)*(1-UC)+UC*E2;        
	        end;

%[t,size(H),size(Z)]
                  
                k = AY / (ESU(t) + V0);		% Kalman Gain
                zt = zt + k'*E;
                %z(t,:) = zt;
		
                if aMode==0
                        %W = W; %nop                  % Schloegl et al. 2003
                elseif aMode==2
                        T(t)=(1-UC)*T(t-1)+UC*(E2-Q(t))/(H'*H);   % Roberts I 1998
                        Z=Z*V(t-1)/Q(t);  
                        if T(t)>0 W=T(t)*eye(MOP); else W=zeros(MOP);end;          
                elseif aMode==5
                        Q_wo = (H'*C*H + V(t-1));                 % Roberts II 1998
                        T(t)=(1-UC)*T(t-1)+UC*(E2-Q_wo)/(H'*H);      
                        if T(t)>0 W=T(t)*eye(MOP); else W=zeros(MOP); end;          
                elseif aMode==6
                        T(t)=(1-UC)*T(t-1)+UC*(E2-Q(t))/(H'*H);      
                        Z=Z*V(t)/Q(t);  
                        if T(t)>0 W=T(t)*eye(MOP); else W=zeros(MOP); end;          
                elseif aMode==11
                        %Z = Z - k*AY';
                        W = sum(diag(Z))*dW;
                elseif aMode==12
                        W = UC*UC*eye(MOP);
                elseif aMode==13
                        W = UC*diag(diag(Z));
                elseif aMode==14
                        W = (UC*UC)*diag(diag(Z));
                elseif aMode==15
                        W = sum(diag(Z))*dW;
                elseif aMode==16
                        W = UC*eye(MOP);               % Schloegl 1998
                elseif aMode==17
      			Z = 0.5*(Z+Z');
       			W = UC*Z;
                elseif aMode==18
       			W = 0.5*UC*(Z+Z');
			%W=W;
                end;

	        Z = Z - k*AY';               % Schloegl 1998
    	else

		V(t) = V0;

    	end;     
        if any(any(isnan(W))), W=UC*Z; end;
            
	z(t,:) = zt;
	Z   = Z + W;               % Schloegl 1998
	SPUR(t)=trace(Z);
end;


if 0,m,
    z(:,1)=M0*z(:,1)./(1-sum(z(:,2:p),2));
end;

REV = mean(e.*e)/mean(y.*y);
if any(~isfinite(Z(:))), REV=inf; end;

