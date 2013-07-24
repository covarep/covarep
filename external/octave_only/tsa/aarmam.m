function [z,e,REV,ESU,V,Z,SPUR] = aarmam(y, Mode, MOP, UC, z0, Z0, V0, W); 
% Estimating Adaptive AutoRegressive-Moving-Average-and-mean model (includes mean term) 
%
% !! This function is obsolete and is replaced by AMARMA
%
% [z,E,REV,ESU,V,Z,SPUR] = aarmam(y, mode, MOP, UC, z0, Z0, V0, W); 
% Estimates AAR parameters with Kalman filter algorithm
% 	y(t) = sum_i(a_i(t)*y(t-i)) + m(t) + e(t) + sum_i(b_i(t)*e(t-i))
%
% State space model
%	z(t) = G*z(t-1) + w(t)    w(t)=N(0,W) 
%	y(t) = H*z(t)   + v(t)	  v(t)=N(0,V)	
%
% G = I, 
% z = [m(t),a_1(t-1),..,a_p(t-p),b_1(t-1),...,b_q(t-q)];
% H = [1,y(t-1),..,y(t-p),e(t-1),...,e(t-q)];
% W = E{(z(t)-G*z(t-1))*(z(t)-G*z(t-1))'}
% V = E{(y(t)-H*z(t-1))*(y(t)-H*z(t-1))'}
%
%
% Input:
%       y	Signal (AR-Process)
%       Mode	determines the type of algorithm
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
% REFERENCE(S): 
% [1] A. Schloegl (2000), The electroencephalogram and the adaptive autoregressive model: theory and applications. 
%     ISBN 3-8265-7640-3 Shaker Verlag, Aachen, Germany. 
%
% More references can be found at 
%     http://pub.ist.ac.at/~schloegl/publications/

%       $Id: aarmam.m 9609 2012-02-10 10:18:00Z schloegl $
%       Copyright (C) 1998-2002,2008,2012 by Alois Schloegl <alois.schloegl@ist.ac.at>
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

%#realonly 
%#inbounds

warning('AARMAM is obsolete. Use AMARMA instead!')

[nc,nr]=size(y);

if nargin<2 Mode=0; 
elseif ischar(Mode) Mode=bin2dec(Mode); 
elseif isnan(Mode) return; end;
if nargin<3 MOP=[0,10,0]; end;
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


e=zeros(nc,1);
V=zeros(nc,1);V(1)=V0;
T=zeros(nc,1);
ESU=zeros(nc,1)+nan;
SPUR=zeros(nc,1)+nan;
z=z0(ones(nc,1),:);

arc=poly((1-UC*2)*[1;1]);b0=sum(arc); % Whale forgetting factor for Mode=258,(Bianci et al. 1997)

dW=UC/MOP*eye(MOP);                % Schloegl


%------------------------------------------------
%	First Iteration
%------------------------------------------------

H = zeros(MOP,1); 
if m, 
        H(1) = 1;%M0; 
        if m~=1,
                fprintf(2,'Warning AARMAM: m must be 0 or 1\n');
		return;        
        end;
end; 
if (p<0) || (q<0) || (round(p)~=p) || (round(q)~=q),
        fprintf(2,'Error AARMAM: p and q must be positive integers\n');
	return;        
end;

E = 0;
Z = Z0;
zt= z0;

A1 = zeros(MOP); A2 = A1;

y_1=0;

%------------------------------------------------
%	Update Equations
%------------------------------------------------

for t=1:nc,
        
        % make measurement matrix
        if 0,
                if t>1, 
                        y_1 = y(t-1);
                end;
                H=[1; y_1; H(m+(1:p-1)'); E(1:min(1,q-1)) ; H(p+m+(1:q-1)')];  % shift y and e
                
        else    % this seem to be slightly faster 
                if t<=p, H(m+(1:t-1)) = y(t-1:-1:1);    % Autoregressive 
                else     H(m+(1:p))   = y(t-1:-1:t-p); 
                end;
                
                if t<=q, H(m+p+(1:t-1)) = e(t-1:-1:1);  % Moving Average
                else     H(m+p+(1:q)) = e(t-1:-1:t-q); 
                end;
        end;
        
        % Prediction Error 
        E = y(t) - zt*H;
        e(t) = E;
        
        if ~isnan(E),
                E2 = E*E;
                AY = Z*H; 
                
                ESU(t) = H'*AY;
                
                if eMode==1
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
                        if E2>ESU(t), 
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
                
                k = AY / (ESU(t) + V0);		% Kalman Gain
                zt = zt + k'*E;
                %z(t,:) = zt;
                
                if aMode==2
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
                        %elseif aMode==17
                        %W=W;
                end;
                
                Z = Z - k*AY';               % Schloegl 1998
        else
                
                V(t) = V0;
                
        end;     
        
        if any(any(isnan(W))), W=UC*Z;end;
        
        z(t,:) = zt;
        Z   = Z + W;               % Schloegl 1998
        SPUR(t) = trace(Z);
end;

REV = mean(e.*e)/mean(y.*y);
if any(~isfinite(Z(:))), REV=inf; end;

