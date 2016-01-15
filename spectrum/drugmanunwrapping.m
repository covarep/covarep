% This script achieves accurate and fast phase unwrapping.
%
% Octave compatible
%
% Description
%  Achieves phase unwrapping based on two mathematical properties: 1) the
%  value of the unwrapped phase at Nyquist frequency, 2) the modified
%  Schur-Cohn algorithm to calculate the number of zeros outside the unit
%  circle. The method has been shown to be both accurate (exact unwrapping)
%  and fast.
%
% Inputs
%  x                : [samples] [Nx1] input signal frame
%  Nfft             : number of samples on which the phase is computed. It
%  is also used as initial resolution for the algorithm. If missing,
%  default value is the next power of 2 of N.
%
% Outputs
%  ph               : [samples] [Nfftx1] unwrapped phase
%
% References
% [1] T. Drugman, Y. Stylianou, Fast and Accurate Phase Unwrapping,
% Interspeech, Dresden Germany, , 2015 
%
%  Publications available at the following link:
%  http://http://tcts.fpms.ac.be/~drugman/files/IS15_PhaseUnwrap.pdf
%
% Copyright (c) 2014 Toshiba
%
% License
%  This code is a part of the GLOAT toolbox with the following
%  licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be


function [ph] = DrugmanUnwrapping(x,Nfft)

if nargin<2
    Nfft=2^(nextpow2(length(x)));
end

if(~isreal(x))
   error('$1 must be real');
end

if size(x,1)<size(x,2)
    x=x';
end

[Nref] = ModifiedSchurCohn(x);

N_OUC=-1000;
if mod(Nfft,2)~=0
    Nfft=Nfft+1;
    display('Nfft not a multiple of 2: corrected')
end

Ind=1;
Nfft_init=Nfft;
while (abs(N_OUC-Nref)>1)
    
   ph = unwrap(angle(fft(x,Nfft)));
   
   N_OUC = -round(ph(length(ph)/2 + 1)/pi);    
   Nfft = 2*Nfft;   
      
   Ind=Ind+1;
   
   if Ind==12

        display('Problem: convergence not reached after 12 iterations')
        pause(0.0000001)       
        break
   end
   
end

Nfft_final=Nfft/2;
ph=ph-ph(1);
ph = ph(1:(Nfft_final/Nfft_init):end);



function [N_OUC,Nder] = ModifiedSchurCohn(P);


while P(1)==0
    P(1)=[];
end

n=length(P)-1;
Pj=P;

k=[];
Nder=0;

for j=n:-1:1
    Success=0;
    Ntry=0;
        
    while Success==0        
                
        if mod(j,50)==0
            Pj=Pj/max(abs(Pj));
        end

        Pj_=flipud(Pj);                
        kj=-Pj(end)/Pj_(end);                       
        Q=Pj+kj*Pj_;            
        
        Q(end)=[];
        itmp=0;
        if length(Q)>1
            while Q(end)==0
                if length(Q)>1
                    Q(end)=[];
                    Nder=Nder+1;
                    itmp=itmp+1;
                else
                    Q=0;
                    break;
                end
            end
        end
        
        Success=1;
        if abs(kj)<1
            Pj=Q;
        elseif abs(kj)>1
            Pj=flipud(Q);
        else
            if Q==0             
                Pj=polyder(Pj);
            else
                if Ntry==0
                    
                    Ctmp=rand(1);
                    for tt=1:itmp+1
                        Q=[0 Q']';
                    end
                    Nder=Nder-itmp;                    
                    
                    Pj=Q+(1/Ctmp*Pj+kj*Ctmp*Pj_);
                    
                    Success=0;
                    Ntry=Ntry+1;
                    
                else
                    Pj=polyder(Pj);
                    Success=0;
                end
                
            end
        end
        
    end
    
    
    k(n-j+1)=kj;
    
    if length(Pj)==1
        break;
    end
end

pos=find(abs(k)>1);

N_OUC=length(pos);
