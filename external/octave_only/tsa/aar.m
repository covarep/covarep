function [a,e,REV,TOC,CPUTIME,ESU] = aar(y, Mode, arg3, arg4, arg5, arg6, arg7, arg8, arg9); 
% Calculates adaptive autoregressive (AAR) and adaptive autoregressive moving average estimates (AARMA)
% of real-valued data series using Kalman filter algorithm.
% [a,e,REV] = aar(y, mode, MOP, UC, a0, A, W, V); 
%
% The AAR process is described as following  
%       y(k) - a(k,1)*y(t-1) -...- a(k,p)*y(t-p) = e(k);
% The AARMA process is described as following  
%       y(k) - a(k,1)*y(t-1) -...- a(k,p)*y(t-p) = e(k) + b(k,1)*e(t-1) + ... + b(k,q)*e(t-q);
%
% Input:
%       y       Signal (AR-Process)
%       Mode    is a two-element vector [aMode, vMode], 
%               aMode determines 1 (out of 12) methods for updating the co-variance matrix (see also [1])
%               vMode determines 1 (out of 7) methods for estimating the innovation variance (see also [1])
%               aMode=1, vmode=2 is the RLS algorithm as used in [2]
%               aMode=-1, LMS algorithm (signal normalized)
%               aMode=-2, LMS algorithm with adaptive normalization  
%                                     
%       MOP     model order, default [10,0] 
%               MOP=[p]         AAR(p) model. p AR parameters
%               MOP=[p,q]       AARMA(p,q) model, p AR parameters and q MA coefficients
%       UC      Update Coefficient, default 0
%       a0      Initial AAR parameters [a(0,1), a(0,2), ..., a(0,p),b(0,1),b(0,2), ..., b(0,q)]
%                (row vector with p+q elements, default zeros(1,p) )
%       A       Initial Covariance matrix (positive definite pxp-matrix, default eye(p))
%	W	system noise (required for aMode==0)
%	V	observation noise (required for vMode==0)
%      
% Output:
%       a       AAR(MA) estimates [a(k,1), a(k,2), ..., a(k,p),b(k,1),b(k,2), ..., b(k,q]
%       e       error process (Adaptively filtered process)
%       REV     relative error variance MSE/MSY
%
%
% Hint:
% The mean square (prediction) error of different variants is useful for determining the free parameters (Mode, MOP, UC) 
%
% REFERENCE(S): 
% [1] A. Schloegl (2000), The electroencephalogram and the adaptive autoregressive model: theory and applications. 
%     ISBN 3-8265-7640-3 Shaker Verlag, Aachen, Germany. 
%
% More references can be found at 
%     http://www.dpmi.tu-graz.ac.at/~schloegl/publications/

%
%	$Id: aar.m 8383 2011-07-16 20:06:59Z schloegl $
%       Copyright (C) 1998-2003 by Alois Schloegl <a.schloegl@ieee.org>
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
%if nc<nr y=y'; end; tmp=nr;nc=nr; nr=tmp;end;

if nargin<2 Mode=0; end;
% check Mode (argument2)
if prod(size(Mode))==2
        aMode=Mode(1);
        vMode=Mode(2);
end;
if any(aMode==(0:14)) && any(vMode==(0:7)), 
        fprintf(1,['a' int2str(aMode) 'e' int2str(vMode) ' ']);
else
        fprintf(2,'Error AAR.M: invalid Mode argument\n');
        return;
end;

% check model order (argument3)
if nargin<3 MOP=[10,0]; else MOP= arg3; end;
if length(MOP)==0 p=10; q=0; MOP=p;
elseif length(MOP)==1 p=MOP(1); q=0; MOP=p;
elseif length(MOP)>=2 p=MOP(1); q=MOP(2); MOP=p+q;
end;

if nargin<4 UC=0; else UC= arg4; end;

a0=zeros(1,MOP); 
A0=eye(MOP);
if nargin>4, 
	if all(size(arg5)==([1,1]*(MOP+1))); 	% extended covariance matrix of AAR parameters 
		a0 = arg5(1,2:size(arg5,2));
		A0 = arg5(2:size(arg5,1),2:size(arg5,2)) - a0'*a0;
	else
		a0 = arg5;  
		if nargin>5 
			A0 = arg6;  
		end;
	end;
end;

if nargin<7, W  = []; else W  = arg7; end;
        
if all(size(W)==MOP), 
        if aMode ~= 0, 
                fprintf(1,'aMode should be 0, because W is given.\n');
        end;
elseif isempty(W),
        if aMode == 0, 
                fprintf(1,'aMode must be non-zero, because W is not given.\n');
        end;
elseif any(size(W)~=MOP), 
        fprintf(1,'size of W does not fit. It must be %i x %i.\n',MOP,MOP);
        return;
end;

if nargin<8, V0 = []; else V0 = arg8; end;
if all(size(V0)==nr), 
        if vMode ~= 0, 
                fprintf(1,'vMode should be 0, because V is given.\n');
        end;
elseif isempty(V0),
        if aMode == 0, 
                fprintf(1,'vMode must be non-zero, because V is not given.\n');
        end;
else 
        fprintf(1,'size of V does not fit. It must be 1x1.\n');
        return;
end;

% if nargin<7 TH=3; else TH = arg7;  end;
%       TH=TH*var(y);
%       TH=TH*mean(detrend(y,0).^2);
MSY=mean(detrend(y,0).^2);

e=zeros(nc,1);
Q=zeros(nc,1);
V=zeros(nc,1);
T=zeros(nc,1);
%DET=zeros(nc,1);
SPUR=zeros(nc,1);
ESU=zeros(nc,1);
a=a0(ones(nc,1),:);
%a=zeros(nc,MOP);
%b=zeros(nc,q);

mu=1-UC; % Patomaeki 1995
lambda=(1-UC); % Schloegl 1996
arc=poly((1-UC*2)*[1;1]);b0=sum(arc); % Whale forgettting factor for Mode=258,(Bianci et al. 1997)

dW=UC/MOP*eye(MOP);                % Schloegl


%------------------------------------------------
%       First Iteration
%------------------------------------------------
Y=zeros(MOP,1);
C=zeros(MOP);
%X=zeros(q,1);
at=a0;
A=A0;
E=y(1);
e(1)=E;
if ~isempty(V0)
        V(1) = V0;
else
        V(1) = (1-UC) + UC*E*E;
end;
ESU(1) = 1; %Y'*A*Y;

A1=zeros(MOP);A2=A1;
tic;CPUTIME=cputime;
%------------------------------------------------
%       Update Equations
%------------------------------------------------
T0=2;

for t=T0:nc,
        
        %Y=[y(t-1); Y(1:p-1); E ; Y(p+1:MOP-1)]
        
        if t<=p Y(1:t-1)=y(t-1:-1:1);           % Autoregressive 
        else    Y(1:p)=y(t-1:-1:t-p); 
        end;
        
        if t<=q Y(p+(1:t-1))=e(t-1:-1:1);       % Moving Average
        else    Y(p+1:MOP)=e(t-1:-1:t-q); 
        end;
        
        % Prediction Error 
        E = y(t) - a(t-1,:)*Y;
        e(t) = E;
        E2=E*E;
        
        AY=A*Y; 
        esu=Y'*AY;
        ESU(t)=esu;
        
        if isnan(E),
                a(t,:)=a(t-1,:);
        else
                V(t) = V(t-1)*(1-UC)+UC*E2;        
                if aMode == -1, % LMS 
                        %       V(t) = V(t-1)*(1-UC)+UC*E2;        
                        a(t,:)=a(t-1,:) + (UC/MSY)*E*Y';
                elseif aMode == -2, % LMS with adaptive estimation of the variance 
                        a(t,:)=a(t-1,:) + UC/V(t)*E*Y';
                        
                else    % Kalman filtering (including RLS) 
                        if vMode==0,            %eMode==4
                                Q(t) = (esu + V0);      
                        elseif vMode==1,            %eMode==4
                                Q(t) = (esu + V(t));      
                        elseif vMode==2,        %eMode==2
                                Q(t) = (esu + 1);          
                        elseif vMode==3,        %eMode==3
                                Q(t) = (esu + lambda);     
                        elseif vMode==4,        %eMode==1
                                Q(t) = (esu + V(t-1));           
                        elseif vMode==5,        %eMode==6
                                if E2>esu 
                                        V(t)=(1-UC)*V(t-1)+UC*(E2-esu);
                                else 
                                        V(t)=V(t-1);
                                end;
                                Q(t) = (esu + V(t));           
                        elseif vMode==6,        %eMode==7
                                if E2>esu 
                                        V(t)=(1-UC)*V(t-1)+UC*(E2-esu);
                                else 
                                        V(t)=V(t-1);
                                end;
                                Q(t) = (esu + V(t-1));           
                        elseif vMode==7,        %eMode==8
                                Q(t) = esu;
                        end;
                        
                        k = AY / Q(t);          % Kalman Gain
                        a(t,:) = a(t-1,:) + k'*E;
                        
                        if aMode==0,                    %AMode=0
                                A = A - k*AY' + W;                   % Schloegl et al. 2003
                        elseif aMode==1,                    %AMode=1
                                A = (1+UC)*(A - k*AY');                   % Schloegl et al. 1997
                        elseif aMode==2,                %AMode=11
                                A = A - k*AY';
                                A = A + sum(diag(A))*dW;
                        elseif aMode==3,                %AMode=5
                                A = A - k*AY' + sum(diag(A))*dW;
                        elseif aMode==4,                %AMode=6
                                A = A - k*AY' + UC*eye(MOP);               % Schloegl 1998
                        elseif aMode==5,                %AMode=2
                                A = A - k*AY' + UC*UC*eye(MOP);
                        elseif aMode==6,                %AMode=2
                                T(t)=(1-UC)*T(t-1)+UC*(E2-Q(t))/(Y'*Y);  
                                A=A*V(t-1)/Q(t);  
                                if T(t)>0 A=A+T(t)*eye(MOP); end;          
                        elseif aMode==7,                %AMode=6
                                T(t)=(1-UC)*T(t-1)+UC*(E2-Q(t))/(Y'*Y);      
                                A=A*V(t)/Q(t);  
                                if T(t)>0 A=A+T(t)*eye(MOP); end;          
                        elseif aMode==8,                %AMode=5
                                Q_wo = (Y'*C*Y + V(t-1));                
                                C=A-k*AY';
                                T(t)=(1-UC)*T(t-1)+UC*(E2-Q_wo)/(Y'*Y);      
                                if T(t)>0 A=C+T(t)*eye(MOP); else A=C; end;          
                        elseif aMode==9,                %AMode=3
                                A = A - (1+UC)*k*AY';
                                A = A + sum(diag(A))*dW;
                        elseif aMode==10,               %AMode=7
                                A = A - (1+UC)*k*AY' + sum(diag(A))*dW;
                        elseif aMode==11,               %AMode=8
                                
                                A = A - (1+UC)*k*AY' + UC*eye(MOP);        % Schloegl 1998
                        elseif aMode==12,               %AMode=4
                                A = A - (1+UC)*k*AY' + UC*UC*eye(MOP);
                        elseif aMode==13
                                A = A - k*AY' + UC*diag(diag(A));
                        elseif aMode==14
                                A = A - k*AY' + (UC*UC)*diag(diag(A));
                        end;
                end;
        end;
end;

%a=a(end,:);
TOC = toc;
CPUTIME = cputime - CPUTIME;
%REV = (e'*e)/(y'*y);

REV = mean(e.*e)./mean(y.*y);
