% MaxP-LP achieves a sparse LP analysis of speech.
%
% Description
%  MaxP-LP is a method allowing maximum-phase modeling for the Linear
%  Prediction of speech. It provides sparser residual signals compared to
%  the conventional approach.
%
% Inputs
%  s             : [samples] [Nx1] input signal (speech signal)
%  fs               : [Hz]      [1x1] sampling frequency
%  order           : [1x1] the total prediction order used for LP analysis
%
% Outputs
%  res              : [samples] [Nx1] the residual excitation signal
%                     obtained by inverse filtering
%
% References
%  [1] T.Drugman, "Maximum Phase Modeling for Sparse Linear Prediction of
%      Speech", IEEE Signal Processing Letters, vol. 21(2), pp. 185-189, 2014.
%
% Copyright (c) 2013 University of Mons, FNRS
%
% License
%  This code will be a part of the GLOAT toolbox with the following
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
% This function is part of the GLOAT toolbox:
% http://tcts.fpms.ac.be/~drugman/Toolbox/
%
% and of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be

function [res] = maxlpresidual(s, fs, order)

L=round(25/1000*fs);          % L is the length of analysis window=25 msec
shift=round(10/1000*fs);       % shift is 5 ms by default

start=1;
stop=start+L;

HannWin = hanning(L+1);
res=zeros(1,length(s));

%% Two possible values for the preemphasis coefficient: -1 and -0.7
Alphas=[-1 -0.7];

while stop<=length(s)
    
    Ind=1;
    Res=[];
    
    for k=1:length(Alphas)
    
        Alpha=Alphas(k);
        
        %% The first LP analysis (of order K-2) estimates the LP coefficients on the pre-emphasized speech signal
        segment=filter([1 Alpha],1,s(start:stop));
        segment=segment.*HannWin;
        [A]=lpc(segment,order-2); 
                    
        segment=s(start:stop);
        segment=segment.*HannWin;        
        segment_orig=segment;
        
        %% The speech signal is inverse-filtered using the coeff. of the first LP analysis
        Res1=filter(A,1,segment);       
        %% The resulting is called ra(n) in the paper
        
        %% First operation of causality inversion
        AntiCausal=flipud(Res1);
        
        %% The second LP analysis (of order 2) and inverse filtering is achieved on ra(-n)
        [A]=lpc(AntiCausal,2);                 
        Res2=filter(A,1,AntiCausal);
        
        %% Final operation of causilty inversion: the time axis is in its original direction
        ResFinal=flipud(Res2);
                
        ResFinal=ResFinal*sqrt(sum(segment_orig.^2)/sum(ResFinal.^2));      
        
        Val(Ind)=GiniMeasure(ResFinal);                
        Res(Ind,:)=ResFinal;        
        Ind=Ind+1;
        
    end

    %% Choose now the value of Alpha (the preemphasis coeff) maximizing the Gini index (here used as a sparsity measure)
    [~,posi]=max(Val);    
    SparseRes=Res(posi,:); 
    
    if Alphas(posi)<-0.7
        %to compensate the polarity inversion due to preemphasis
        SparseRes=-SparseRes;
    end

    res(start:stop)=res(start:stop)+SparseRes;
    
    start=start+shift;
    stop=stop+shift;
end

res=res/max(abs(res));


function [val] = GiniMeasure(x)
x=abs(x);
x=sort(x);

S=0;
N=length(x);

Ab=sum((x));
for k=1:N
    S=S+x(k)/Ab*((N-k+1/2)/N);
end

val=1-2*S;

return
