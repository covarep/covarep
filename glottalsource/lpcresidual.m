% Function to derive the Linear Prediction residual signal
%
% Octave compatible
%
% Description
%  Function to derive the Linear Prediction residual signal
%
% Inputs
%  x               : [samples] [Nx1] Input signal
%  L               : [samples] [1x1] window length (e.g., 25ms =>  25/1000*fs)
%  shift           : [samples] [1x1] window shift (e.g., 5ms => 5/1000*fs)
%  order           : [samples] [1x1] Order of Linear Prediction
%
% Outputs
%  res             : [samples] [Mx1] Linear Prediction residual
%  LPCcoeff        : [samples] [order+1xM] Linear Prediction coefficients
%
% Example
%  [res,LPCcoeff] = lpcresidual(x,L,shift,order);
%
% Copyright (c) 2011 University of Mons, FNRS
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
%  Thomas Drugman, TCTS Lab.
%
% $Id <info set by the versioning system> $

function [res,LPCcoeff] = lpcresidual(x,L,shift,order)


%% Initial settings
x=x(:);
shift = round(shift);
order = round(order);

start=1;
L = round(L);
stop=start+L;    

% Allocate space
res=zeros(1,length(x));
LPCcoeff=zeros(order+1,round(length(x)/shift));

%% Do processing
n=1;
win = hanning(stop-start+1);
while stop<length(x)
     segment=x(start:stop);
     segment=segment.*win;
        
     A=lpc(segment,order);
     LPCcoeff(:,n)=A(:);
        
     inv=filter(A,1,segment);
       
     inv=inv*sqrt(sum(segment.^2)/sum(inv.^2));

     res(start:stop)=res(start:stop)+inv'; % Overlap and add

     % Increment
     start=start+shift;
     stop=stop+shift;
     n=n+1;
end

res=res/max(abs(res)); % Normalise amplitude

