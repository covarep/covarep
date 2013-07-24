function [a,VAR,S,a_aux,b_aux,e_aux,MLE,pos] = rmle(arg1,arg2); 
% RMLE estimates AR Parameters using the Recursive Maximum Likelihood 
% Estimator according to [1]
% 
% Use: [a,VAR]=rmle(x,p)
    % Input: 
                % x is a column vector of data
                % p is the model order
    % Output:
                % a is a vector with the AR parameters of the recursive MLE
                % VAR is the excitation white noise variance estimate
%
% Reference(s):
% [1] Kay S.M., Modern Spectral Analysis - Theory and Applications. 
%       Prentice Hall, p. 232-233, 1988. 
%
                
%       $Id: rmle.m 9609 2012-02-10 10:18:00Z schloegl $
%	Copyright (C) 2004 by  Jose Luis Gutierrez <jlg@gmx.at>
%	Grupo GENESIS - UTN - Argentina
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


x=arg1*1e-6;
p=arg2;

N=length(x);
S=zeros(p+1,p+1);
a_aux=zeros(p+1,p);, a_aux(1,:)=1;
b_aux=ones(p+1,p);
e_aux=zeros(p,1);, p_aux=zeros(p,1);
MLE=zeros(3,1);
pos=1;

for i=0:p
    for j=0:p      
        for n=0:N-1-i-j
            S(i+1,j+1)=S(i+1,j+1)+x(n+1+i)*x(n+1+j);
        end
    end
end

e0=S(1,1);
c1=S(1,2);
d1=S(2,2);
coef3=1;
coef2=((N-2)*c1)/((N-1)*d1);
coef1=-(e0+N*d1)/((N-1)*d1);
ti=-(N*c1)/((N-1)*d1);
raices=roots([coef3 coef2 coef1 ti]);
for o=1:3
    if raices(o)>-1 && raices(o)<1
        a_aux(2,1)=raices(o);    
        b_aux(p+1,1)=raices(o);
    end
end
e_aux(1,1)=S(1,1)+2*a_aux(2,1)*S(1,2)+(a_aux(2,1)^2)*S(2,2);
p_aux(1,1)=e_aux(1,1)/N;

for k=2:p
    Ck=S(1:k,2:k+1);
    Dk=S(2:k+1,2:k+1);
    ck=a_aux(1:k,k-1)'*Ck*b_aux(p+1:-1:p+2-k,k-1);
    dk=b_aux(p+1:-1:p+2-k,k-1)'*Dk*b_aux(p+1:-1:p+2-k,k-1);
    coef3re=1;
    coef2re=((N-2*k)*ck)/((N-k)*dk);
    coef1re=-(k*e_aux(k-1,1)+N*dk)/((N-k)*dk);
    tire=-(N*ck)/((N-k)*dk);
    raices=roots([coef3re coef2re coef1re tire]);
    for o=1:3
        if raices(o,1)>-1 && raices(o,1)<1
            MLE(o,1)=((1-raices(o)^2)^(k/2))/(((e_aux(k-1)+2*ck*raices(o)+dk*(raices(o)^2))/N)^(N/2));
        end
    end
    [C,I]=max(MLE);
    k_max=raices(I);
    for i=1:k-1
        a_aux(i+1,k)=a_aux(i+1,k-1)+k_max*a_aux(k-i+1,k-1);
    end
    a_aux(k+1,k)=k_max;
    b_aux(p+1-k:p+1,k)=a_aux(1:k+1,k);
    e_aux(k,1)=e_aux(k-1,1)+2*ck*k_max+dk*k_max^2;
    p_aux(k,1)=e_aux(k,1)/N;
end

a=a_aux(:,p)';
VAR=p_aux(p)*1e12;
