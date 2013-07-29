% get_MBS.m
% Function to derive mean-based signal (MBS) as is derived in Thomas
% Drugman's SEDREAMS algorithm
%
% Octave compatible
%
% Description
% Function to derive mean-based signal (MBS) as is derived in Thomas
% Drugman's SEDREAMS algorithm
%
% Inputs
%  x               : [samples] [Nx1] Speech signal
%  fs              : [Hz]      [1x1] sampling frequency
%  T0mean          : [samples] [1x1] Average glottal period
%
% Outputs
%  MBS             : [samples] [Nx1] Mean-based signal
%
% Example
%  MBS = get_MBS(x,fs,T0mean)
%
% References
%  [1] Kane, J., Gobl, C., (2013) `Evaluation of glottal closure instant 
%       detection in a range of voice qualities', Speech Communication 
%       55(2), pp. 295-314.
%
%
% Copyright (c) 2013 Trinity College Dublin
%
% License
%  This code is a part of the Voice Analysis Toolkit with the following
%  licence:
%  The software product (and any modifications thereof) is distributed under 
%  a dual licence: an open source license for individual, non-commercial 
%  purposes, and a commercial license. The opensource licence under which 
%  the product is distributed is GNU GPL v2. For individual users, this 
%  licence suits their use as these are not distributing proprietary 
%  modifications, additions to, or derivatives of the product and don't 
%  require legal protection of a commercial licence. For commercial users, 
%  where open source does not meet their requirements, we offer commercial 
%  licensing of the product. A commercial license permits customers to 
%  modify, add or produce derivative products without the obligation of 
%  making the subsequent code open source. For more information regarding 
%  our commercial licence, please contact john.whelan@tcd.ie
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
%
% $Id <info set by the versioning system> $

function MBS = get_MBS(x,fs,T0mean)

%% Initial settings
MBS=zeros(1,length(x));
halfL=round((1.6*T0mean(1))/2);

StepExp=3;
Step=2^StepExp;
HP_p=70; % High-pass pass
HP_s=10; % High-pass stop
sm_len=7;

%% Calculate mean-based signal
for m=halfL+1:Step:length(x)-halfL
    
    if length(T0mean)==1
        halfL=round((1.7*T0mean(1))/2);
    else halfL=round((1.7*T0mean(m))/2);
    end
    BlackWin=blackman(2*halfL+1);
    
    start=round(m-halfL);
    stop=round(m+halfL);
    if stop > length(x)
        break
    end

    if start > 0
        vec=x(start:stop);

        vec=vec.*BlackWin;
        MBS(m)=mean(vec);
    end
end

t=find(MBS~=0);
MBS=interp1(t,MBS(t),1:length(x));
MBS(isnan(MBS))=0;
MBS=zeroPhaseHPFilt(MBS,fs,HP_p,HP_s,0);
MBS=MBS/max(MBS);
MBS=smooth(MBS,sm_len);

return


function y = zeroPhaseHPFilt(x,fs,f_p,f_s,plots)

% Function to use a forward and backwards butterworth high pass filter to
% ensure zero phase 
%
% USAGE:    
%       Input:
%             x  : input signal
%             fs : sampling frequency (in Hz)
%             f_p : pass-band (in Hz)
%             f_s : stop band (in Hz)
%             plots : input 1 for a plot of the filters frequency response.
%
%       Output:
%             y  : filtered signal
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rp = 0.5;
Rs = 40;
Wp=f_p/(fs/2);
Ws=f_s/(fs/2);

[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[b,a]=butter(n,Wn,'high');

y = filtfilt(b,a,x);

if nargin > 4
    if plots==1
        freqz(b,a)
    end
end

return
