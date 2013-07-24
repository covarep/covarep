function [x,e,Kalman,Q2] = mvaar(y,p,UC,mode,Kalman)
% Multivariate (Vector) adaptive AR estimation base on a multidimensional
% Kalman filer algorithm. A standard VAR model (A0=I) is implemented. The 
% state vector is defined as X=(A1|A2...|Ap) and x=vec(X')
%
% [x,e,Kalman,Q2] = mvaar(y,p,UC,mode,Kalman)
%
% The standard MVAR model is defined as:
%
%		y(n)-A1(n)*y(n-1)-...-Ap(n)*y(n-p)=e(n)
%
%	The dimension of y(n) equals s 
%	
%	Input Parameters:
%
% 		y			Observed data or signal 
% 		p			prescribed maximum model order (default 1)
%		UC			update coefficient	(default 0.001)
%		mode	 	update method of the process noise covariance matrix 0...4 ^
%					correspond to S0...S4 (default 0)
%
%	Output Parameters
%
%		e			prediction error of dimension s
%		x			state vector of dimension s*s*p
%		Q2			measurement noise covariance matrix of dimension s x s
%

%       $Id: mvaar.m 5090 2008-06-05 08:12:04Z schloegl $
% 	Copyright (C) 2001-2002 Christian Kasess  
% 	Copyright (C) 2003, 2008 Alois Schloegl    
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


if nargin<4,
        mode=0;
end;
if nargin<3,
        UC=0.001   
end;
if nargin<2,
        p=1;
end
if nargin<1,
        fprintf(2,'No arguments supplied\n');
        return   
end;

if ~any(mode==(0:4))
        fprintf(2,'Invalid mode (0...4)\n');
        return   
end;


[M,LEN] = size(y');		%number of channels, total signal length
L = M*M*p;

if LEN<(p+1),
        fprintf(2,'Not enough observed data supplied for given model order\n');
        return   
end

ye = zeros(size(y));	%prediction of y

if nargout>1,
        x=zeros(L,LEN);
end;
if nargout>3,  
        Q2=zeros(M,M,LEN);
end



if nargin<5,
        %Kalman Filter initialsiation (Kp (K predicted or a-priori) equals K(n+1,n) )
        Kalman=struct('F',eye(L),'H',zeros(M,L),'G',zeros(L,M),'x',zeros(L,1),'Kp',eye(L),'Q1',eye(L)*UC,'Q2',eye(M),'ye',zeros(M,1));
        end;

upd = eye(L)/L*UC;		%diagonal matrix containing UC


if(mode==3)
        Block=kron(eye(M),ones(M*p));
        
elseif(mode==4)
        index=[];
        Block1=[];
        Block0=[];
        for i=1:M,
                index=[index ((i-1)*M*p+i:M:i*M*p)];
                mone=eye(M);
                mone(i,i)=0;
                mzero=eye(M)-mone;
                Block1=Blkdiag(Block1,kron(eye(p),mone));
                Block0=Blkdiag(Block0,kron(eye(p),mzero));
        end;
end;


for n = 2:LEN,
        if(n<=p)
                Yr=[y(n-1:-1:1,:)' zeros(M,p-n+1)];	%vector of past observations
                Yr=Yr(:)';
        else
                Yr=y(n-1:-1:n-p,:)';						%vector of past observations
                Yr=Yr(:)';
        end
        
        %Update of measurement matrix
        Kalman.H=kron(eye(M),Yr);
        
        
        %calculate prediction error
        ye(n,:)=(Kalman.H*Kalman.x)';
        err=y(n,:)-ye(n,:);
        
        if ~any(isnan(err(:))),
                %update of Q2 using the prediction error of the previous step
                Kalman.Q2=(1-UC)*Kalman.Q2+UC*err'*err;
                
                
                KpH=Kalman.Kp*Kalman.H';
                HKp=Kalman.H*Kalman.Kp;
                
                %Kalman gain
                Kalman.G=KpH*inv(Kalman.H*KpH+Kalman.Q2);
                
                %calculation of the a-posteriori state error covariance matrix
                %K=Kalman.Kp-Kalman.G*KpH'; Althouh PK is supposed to be symmetric, this operation makes the filter unstable
                K=Kalman.Kp-Kalman.G*HKp; 
                
                %mode==0 no update of Q1
                %update of Q1 using the predicted state error cov matrix
                if(mode==1)      
                        Kalman.Q1=diag(diag(K)).*UC;
                elseif(mode==2)
                        Kalman.Q1=upd*trace(K);
                elseif(mode==3)
                        Kalman.Q1=diag(sum((Block*diag(diag(K)))'))/(p*M)*UC;
                elseif(mode==4)
                        avg=trace(K(index,index))/(p*M)*UC;
                        Kalman.Q1=Block1*UC+Block0*avg;
                end
                
                %a-priori state error covariance matrix for the next time step
                Kalman.Kp=K+Kalman.Q1;
                
                %current estimation of state x
                Kalman.x=Kalman.x+Kalman.G*(err)';
        end; % isnan>(err)   
        
        if nargout>1,
                x(:,n) = Kalman.x;
        end;
        if nargout>3,
                Q2(:,:,n)=Kalman.Q2;   
        end;
end;

e = y - ye;
x = x';

