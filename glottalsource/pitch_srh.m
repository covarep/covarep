% SRH is a robust pitch tracker.
%
% Octave compatible
%
% Description
%  The Summation of the Residual Harmonics (SRH) method is described in [1].
%  This algorithm exploits a criterion taking into the strength of the
%  harmonics and subharmonics of the residual excitation signal in order to
%  determine both voicing decision and F0 estimates. It is shown in [1] to
%  be particularly interesting in adverse conditions (low SNRs with various
%  types of additive noises).
%
%
% Inputs
%  wave            : [samples] [Nx1] input signal (speech signal)
%  fs              : [Hz]      [1x1] sampling frequency
%  f0min           : [Hz] [1x1] minimum possible F0 value
%  f0max           : [Hz]  [1x1] maximum possible F0 value
%  hopsize         : [ms]  [1x1] time interval between two consecutive
%                    frames (i.e. defines the rate of feature extraction).
%
% Outputs
%  F0s             : vector containing the F0 estimates (values are
%                    provided even in unvoiced parts).
%  VUVDecisions    : vector containing the binary voicing decisions.
%  SRHVal          : vector containing the SRH values (according the
%                    harmonic criterion - voicing decision are derived from
%                    these values by simple thresholding).
%  time            : [s] Analysis instants of the features described above.
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] T.Drugman, A.Alwan, "Joint Robust Voicing Detection and Pitch Estimation
%      Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011.
%      Publication available at the following link:
%      http://tcts.fpms.ac.be/~drugman/files/IS11-Pitch.pdf
%
% Copyright (c) 2011 University of Mons, FNRS
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
% This function is part of the Common Speech Processing Repository (TODO)
% TODO URL
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be

function [F0s,VUVDecisions,SRHVal,time] = pitch_srh(wave,fs,f0min,f0max,hopsize)


if f0max<=f0min
    display('You look funny! Your f0min should be lower than f0max!!')
    pause(0.00001)
end
    

if fs>16000
    wave=resample(wave,16000,fs);
    fs=16000;
end

if nargin < 5
    hopsize=10;
end

LPCorder=round(3/4*fs/1000);
Niter=2;

[res] = lpcresidual(wave,round(25/1000*fs),round(5/1000*fs),LPCorder);


%% Estimate the pitch track in 2 iterations
for Iter=1:Niter   

    [F0s,SRHVal,time] = SRH_EstimatePitch(res',fs,f0min,f0max,hopsize);     
    
    posiTmp=find(SRHVal>0.1);   
    
    if length(posiTmp)>1
        
        F0meanEst=median(F0s(posiTmp));        
               
        f0min=round(0.5*F0meanEst);
        f0max=round(2*F0meanEst);        
    end
    
end

time = (time-1)/fs;

% Voiced-Unvoiced decisions are derived from the value of SRH (Summation of
% Residual Harmonics)
VUVDecisions=zeros(1,length(F0s));
Threshold=0.07;

if std(SRHVal)>0.05
    Threshold=0.085;
end

pos= SRHVal>Threshold;
VUVDecisions(pos)=1;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F0s,SRHVal,time] = SRH_EstimatePitch(sig,Fs,F0min,F0max,Hopsize)

time=[];
start=1;
stop=round(100/1000*Fs);
shift=round(Hopsize/1000*Fs);

Nframes=floor((length(sig)-stop)/shift) + 1;
SRHTot=zeros(Nframes,F0max);
F0s = zeros(1,Nframes);
SRHVal = zeros(1,Nframes);

BlackWin=blackman(stop-start+1);

index=1;
while stop<=length(sig)
%      time(index) = start;         % original
    time(index) = 0.5*(start+stop); % degottex
    seg=sig(start:stop);
    seg=seg.*BlackWin;    
    seg=seg-mean(seg);

    Spec=fft(seg,Fs);
    Spec=abs(Spec(1:Fs/2));    
    Spec=Spec/sqrt(sum(Spec.^2));
        
    SRHs=zeros(1,F0max);    
    
    % SRH spectral criterion
    for freq=F0min:F0max
        SRHs(freq)=(Spec(freq)+Spec(2*freq)+Spec(3*freq)+Spec(4*freq)+Spec(5*freq))-(Spec(round(1.5*freq))+Spec(round(2.5*freq))+Spec(round(3.5*freq))+Spec(round(4.5*freq)));
    end
    
    SRHTot(index,:)=SRHs;

    [maxi,posi]=max(SRHs);
    F0frame=posi;
        
    F0s(index)=F0frame;
    SRHVal(index)=SRHs(F0frame);
    
    start=start+shift;
    stop=stop+shift;
    index=index+1;
end

 NframesToAdd=round((100/1000*Fs)/(2*shift)); % commented by degottex
 
 F0s=[zeros(1,NframesToAdd) F0s zeros(1,NframesToAdd)]; % commented by degottex
 SRHVal=[zeros(1,NframesToAdd) SRHVal zeros(1,NframesToAdd)]; % commented by degottex

 for k=1:NframesToAdd     
     time=[time(1)-shift time time(end)+shift];
 end
 
return

