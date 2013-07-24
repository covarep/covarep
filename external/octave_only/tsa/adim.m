function [IR, CC, D] = adim(U, UC, IR, CC, arg5); 
% ADIM adaptive information matrix. Estimates the inverse
%   correlation matrix in an adaptive way. 
%
% [IR, CC] = adim(U, UC [, IR0 [, CC0]]); 
%   U 	input data  
%   UC 	update coefficient 0 < UC << 1
%   IR0	initial information matrix
%   CC0 initial correlation matrix
%   IR	information matrix (inverse correlation matrix)
%   CC  correlation matrix
% 	
%  The algorithm uses the Matrix Inversion Lemma, also known as 
%     "Woodbury's identity", to obtain a recursive algorithm.  
%     IR*CC - UC*I should be approx. zero. 
%
% Reference(s):
% [1] S. Haykin. Adaptive Filter Theory, Prentice Hall, Upper Saddle River, NJ, USA 
%     pp. 565-567, Equ. (13.16), 1996.


%       $Id: adim.m 5090 2008-06-05 08:12:04Z schloegl $
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



[ur,p] = size(U);

Mode_E = 1;
%if nargin < 4,
%        m  = zeros(size(U)+[0,Mode_E]);
%end;

if nargin<2,
        fprintf(2,'Error ADIM: missing update coefficient\n');
        return;
else	
        if ~((UC > 0) & (UC <1)),
                fprintf(2,'Error ADIM: update coefficient not within range [0,1]\n');
                return;
        end;
        if UC > 1/p,
                fprintf(2,'Warning ADIM: update coefficient should be smaller than 1/number_of_dimensions\n');
        end;
end;

if nargin<3,
        IR = [];
end;
if nargin<4,
        CC = [];
end;
if nargin<5,
        arg5 = 6;
end;
if isempty(IR),
        IR = eye(p+Mode_E);
end;
if isempty(CC),
        CC = eye(p+Mode_E);
end;

D = zeros(ur,(p+Mode_E)^2);
%D = zeros(ur,1);
W  = eye(p+Mode_E)*UC/(p+Mode_E);
W2 = eye(p+Mode_E)*UC*UC/(p+Mode_E);

for k = 1:ur,
        if ~Mode_E,
                % w/o mean term 
                d  = U(k,:);	
        else
                % w/ mean term 
                d  = [1,U(k,:)];
        end;
        
        if ~any(isnan(d)),
                CC = (1-UC)*CC + UC*(d'*d);
                
                %K = (1+UC)*IR*d'/(1+(1+UC)*d*IR*d');
                v  = IR*d';
                %K = v/(1-UC+d*v);
                
                if arg5==0;  		% original version according to [1], can become unstable 
                        IR = (1+UC)*IR - (1+UC)/(1-UC+d*v)*v*v';
                        
                elseif arg5==1;  	% this provides IR*CC == I, this seems to be stable 
                        IR = IR - (1+UC)/(1-UC+d*v)*v*v' + sum(diag(IR))*W;
                        
                elseif arg5==6;  	% DEFAULT: 
                        IR = ((1+UC)/2)*(IR+IR') - ((1+UC)/(1-UC+d*v))*v*v';
                        
                        % more variants just for testing, do not use them. 
                elseif arg5==2;
                        IR = IR - (1+UC)/(1-UC+d*v)*v*v' + sum(diag(IR))*W2;
                        
                elseif arg5==3;
                        IR = IR - (1+UC)/(1-UC+d*v)*v*v' + eps*eye(p+Mode_E);
                        
                elseif arg5==4;
                        IR = (1+UC)*IR - (1+UC)/(1-UC+d*v)*v*v' + eps*eye(p+Mode_E);
                        
                elseif arg5==5;
                        IR = IR - (1+UC)/(1-UC+d*v)*v*v' + eps*eye(p+Mode_E);
                        
                end;
        end;
        
        %D(k) = det(IR);
        D(k,:) = IR(:)';
end;

IR = IR / UC;
if Mode_E,
        IR(1,1) = IR(1,1) - 1;
end;


