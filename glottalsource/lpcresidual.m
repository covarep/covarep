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
% Copyright (c) TODO
%
% License
%  TODO
%
% This function is part of the Common Speech Processing Repository (TODO)
% TODO URL
% 
% Author 
%  Thomas Drugman, TCTS Lab.
%
% $Id <info set by the versioning system> $

function [res,LPCcoeff] = lpcresidual(x,L,shift,order)


%% Initial settings
x=x(:);

start=1;
stop=start+L;    

% Allocate space
res=zeros(1,length(x));
LPCcoeff=zeros(order+1,round(length(x)/shift));


%% Do processing
n=1;
while stop<length(x)

     segment=x(start:stop);
     segment=segment.*hanning(length(segment));
        
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

