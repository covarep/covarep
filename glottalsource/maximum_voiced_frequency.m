% Maximum voiced frequency estimation
%
%
% Description
%  This technique estimates the maximum voiced frequency using the 5
%  proposed approaches: AS-IHPC, AS-IHPC-ICPC, AS, IHPC and ICPC. You can
%  see these 5 estimates as 5 alternatives. On overall, AS-IHPC and
%  AS-IHPC-ICPC provide the best estimates. These techniques are fully
%  described in [1].

%
% Inputs
%  wave         : [samples] [Nx1] input signal (speech signal)
%  Fs           : [Hz] Sampling frequency
%  f0           : [Hz] Vector containing the F0 estimates (0 values are
%                      provided in unvoiced parts).
%  t            : [s]  Analysis instants of the f0 values. These will
%                      also be the analysis instants for MVF estimates
%
%
% Outputs
%  [AS_IHPC,AS_IHPC_ICPC,AS,IHPC,ICPC]  : these are 5 vectors containing
%                           the MVF estimates. Values are provided
%                           synchronously with the f0 estimates.
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%  See also http://tcts.fpms.ac.be/~drugman/Toolbox/
%
% References
%  [1] T.Drugman, Y. Stylianou, "Maximum Voiced Frequency Estimation:
%  Exploiting Amplitude and Phase Spectra", IEEE Signal Processing Letters,
%  2014.
%  http://tcts.fpms.ac.be/~drugman/files/SPL-MVF.pdf
%
% Copyright (c) 2014 Toshiba Cambridge Research Laboratory
%
% License
%  This code will be part of the GLOAT toolbox (http://tcts.fpms.ac.be/~drugman/Toolbox/)
%  with the following licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is also be part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be

function [AS_IHPC,AS_IHPC_ICPC,AS,IHPC,ICPC] = maximum_voiced_frequency(wave,Fs,f0,t)

t = round(t*Fs)+1; % From [s] to [samples]

AS=[];
IHPC=[];
ICPC=[];
AS_IHPC=[];
AS_IHPC_ICPC=[];

for n=1:length(t)
    
    if f0(n)>0
        
        %% COMPUTE THE AMPLITUDE AND PHASE SPECTRA
        
        T0=round(Fs/f0(n));
        Start=t(n)-2*T0;
        Stop=t(n)+2*T0;
        
        if Start<1
            Start=1;
        end
        
        if Stop>length(wave)
            Stop=length(wave);
        end
           
        Seg=wave(Start:Stop);        
        Seg=Seg.*hanning(length(Seg));
        
        Spec1=fft(Seg,Fs);
        GD1=-diff(((unwrap(angle(Spec1(1:Fs/2))))));
        Phi1=(unwrap(angle(Spec1(1:Fs/2))));
        Spec1_a=20*log10(abs(Spec1(1:Fs/2)));
        
        Start=t(n)-1*T0;
        Stop=t(n)+3*T0;
        
        if Start<1
            Start=1;
        end
        
        if Stop>length(wave)
            Stop=length(wave);
        end
        
        
        Seg=wave(Start:Stop);
        Seg=Seg.*hanning(length(Seg));
        Spec=fft(Seg,Fs);
        Phi2=(unwrap(angle(Spec(1:Fs/2))));
        
        DPhi=(Phi1-Phi2);
        
        Freq=1:Fs/2;
        Del=2*pi*Freq*T0/Fs;
        DPhi=DPhi+Del';
        
        F0=f0(n);
        Harm=[];
        
        %% EXTRACT THE FEATURES
        
        Candidates=[];
        Candidates2=[];
        Candidates3=[];
        k=1;
        MidF0=round(k*F0);
        while MidF0<(Fs/2-2*F0)
            
            MidF0=round(k*F0);
            Vec=Spec1_a(MidF0-10:MidF0+10);
            [maxi,posi]=max(Vec);
            MidF0=MidF0-10+posi-1;
            Harm(k)=MidF0;
            F0=MidF0/k;
            
            Vec_N=[Spec1_a(MidF0-round(0.5*F0):MidF0-round(0.2*F0))' Spec1_a(MidF0+round(0.2*F0):MidF0+round(0.5*F0))'];            
            LevelN=mean(Vec_N);
            
            Vec_S=[Spec1_a(MidF0-round(0.2*F0):MidF0+round(0.2*F0))];            
            LevelS=mean(Vec_S);
            
            Candidates(k)=LevelS-LevelN;
            
            if k>1
                Candidates2(k)=(GD1(Harm(k))-GD1(Harm(k-1)))*(Harm(k)-Harm(k-1));
            else
                Candidates2(k)=0;
            end
            
            if k>1
                Candidates3(k)=(wrap(DPhi(Harm(k))-DPhi(Harm(k-1))))*(Harm(k)-Harm(k-1));
            else
                Candidates3(k)=0;
            end
            
            % The multiplications by (Harm(k)-Harm(k-1)) is because linear
            % phase has not been removed, and the linear phase depends upon
            % the window length, which is proportional to the pitch period.
            
            k=k+1;
            
        end
        
        %% CALCULATE THE LIKELIHOODS
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Using only the magnitude spectrum (AS)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Parameters of the Gaussian
        Std1=4;
        Mu1=9.61;
        Std2=4.35;
        Mu2=1.75;
        
        Probas1=[];
        for k=1:length(Candidates)
            p1=gaussmf_new(Candidates(k),[Std1 Mu1]);
            p2=gaussmf_new(Candidates(k),[Std2 Mu2]);
            p1=p1/(p1+p2);
            Probas1(k)=p1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Using only the inter-harmonics phase coherence (IHPC)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Std1=0.9;
        Mu1=0;
        Std2=3.66;
        Mu2=0;
        
        Alpha=gaussmf_new(0,[Std1 Mu1])/(gaussmf_new(0,[Std1 Mu1])+gaussmf_new(0,[Std2 Mu2]));
        % Alpha is here used to force Probas2 to be 1 in 0. It is optional
        % and has little impact on the results.
        
        Probas2=[];
        for k=1:length(Candidates2)
            p1=gaussmf_new(Candidates2(k),[Std1 Mu1]);
            p2=gaussmf_new(Candidates2(k),[Std2 Mu2]);
            
            p1=(1/Alpha)*p1/(p1+p2);
            Probas2(k)=p1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Using only the inter-cycle phase coherence (ICPC)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Std1=48.9;
        Mu1=0.72;
        Std2=179.75;
        Mu2=-1.54;
        
        Alpha=gaussmf_new(0,[Std1 Mu1])/(gaussmf_new(0,[Std1 Mu1])+gaussmf_new(0,[Std2 Mu2]));
        % Alpha is here used to force Probas2 to be 1 around 0. It is optional
        % and has little impact on the results.
        
        Probas3=[];
        for k=1:length(Candidates3)
            p1=gaussmf_new(Candidates3(k),[Std1 Mu1]);
            p2=gaussmf_new(Candidates3(k),[Std2 Mu2]);
            
            p1=(1/Alpha)*p1/(p1+p2);
            Probas3(k)=p1;
        end
        
        
        %% TAKING THE MVF DECISION ACCORDING TO THE ML CRITERION
        
        %% AS
        TotProba=[];
        for k=1:length(Candidates)
            TotProba(k)=1;
            for k2=1:k
                TotProba(k)=TotProba(k)*Probas1(k2);
            end
            for k2=k+1:length(Candidates)
                TotProba(k)=TotProba(k)*(1-Probas1(k2));
            end
        end
        
        [maxi,posi]=max(TotProba);
        if posi<=length(Harm)
            AS(n)=Harm(posi);
        else
            AS(n)=Harm(end);
        end
        
        %%%%%%%%%
        %% IHPC        
        TotProba=[];
        for k=1:length(Candidates)
            TotProba(k)=1;
            for k2=1:k
                TotProba(k)=TotProba(k)*Probas2(k2);
            end
            for k2=k+1:length(Candidates)
                TotProba(k)=TotProba(k)*(1-Probas2(k2));
            end
        end
        
        [maxi,posi]=max(TotProba);
        if posi<=length(Harm)
            IHPC(n)=Harm(posi);
        else
            IHPC(n)=Harm(end);
        end
        
        %%%%%%%%
        %% ICPC
        TotProba=[];
        for k=1:length(Candidates3)
            TotProba(k)=0;
            for k2=1:k
                TotProba(k)=TotProba(k)+log10(Probas3(k2));
            end
            for k2=k+1:length(Candidates3)
                TotProba(k)=TotProba(k)+log10((1-Probas3(k2)));
            end
        end
        
        [maxi,posi]=max(TotProba);
        if posi<=length(Harm)
            ICPC(n)=Harm(posi);
        else
            ICPC(n)=Harm(end);
        end
        
        %%%%%%%%
        %% AS-IHPC
        Probas12=(Probas1.*Probas2)./((Probas1.*Probas2)+((1-Probas1).*(1-Probas2)));
        
        TotProba=[];
        for k=1:length(Candidates)
            TotProba(k)=1;
            for k2=1:k
                TotProba(k)=TotProba(k)*Probas12(k2);
            end
            for k2=k+1:length(Candidates)
                TotProba(k)=TotProba(k)*(1-Probas12(k2));
            end
        end
        
        [maxi,posi]=max(TotProba);
        if posi<=length(Harm)
            AS_IHPC(n)=Harm(posi);
        else
            AS_IHPC(n)=Harm(end);
        end
                
        
        %%%%%%%
        %% AS-IHPC-ICPC
        Probas123=(Probas1.*Probas2.*Probas3)./((Probas1.*Probas2.*Probas3)+((1-Probas1).*(1-Probas2).*(1-Probas3)));
        
        TotProba=[];
        for k=1:length(Candidates)
            TotProba(k)=1;
            for k2=1:k
                TotProba(k)=TotProba(k)*Probas123(k2);
            end
            for k2=k+1:length(Candidates)
                TotProba(k)=TotProba(k)*(1-Probas123(k2));
            end
        end
        
        [maxi,posi]=max(TotProba);
        if posi<=length(Harm)
            AS_IHPC_ICPC(n)=Harm(posi);
        else
            AS_IHPC_ICPC(n)=Harm(end);
        end
        
    else
        AS(n)=0;
        IHPC(n)=0;
        ICPC(n)=0;
        AS_IHPC(n)=0;
        AS_IHPC_ICPC(n)=0;
    end
  
end


%% APPLY THE POST-PROCESS ON MVF TRAJECTORIES
[AS] = PostProcess(AS);
[IHPC] = PostProcess(IHPC);
[ICPC] = PostProcess(ICPC);
[AS_IHPC] = PostProcess(AS_IHPC);
[AS_IHPC_ICPC] = PostProcess(AS_IHPC_ICPC);



function [Starts,Stops] = DefineBoundsOfSequences(Sequence)

Starts=[];
Stops=[];

Ind=1;
if Sequence(1)>0;
    Starts(Ind)=1;
end


for k=2:length(Sequence)
    
    if (Sequence(k-1)==0)&&(Sequence(k)>0)
        Starts(Ind)=k;
    elseif (Sequence(k-1)>0)&&(Sequence(k)==0)
        Stops(Ind)=k-1;
        Ind=Ind+1;
    end
end

if length(Starts)>length(Stops)
    Stops=[Stops length(Sequence)];
end





function [Fm] = PostProcess(Fm)

pos0=find(Fm==0);

pos=find((Fm>0)&(Fm<1000));
Fm(pos)=1000;

[Starts,Stops] = DefineBoundsOfSequences(Fm);

for k=1:length(Starts)
    StartTmp=Stops(k)-6;
    StopTmp=Stops(k)+3;
    if StartTmp<1
        StartTmp=1;
    end
    if StopTmp>length(Fm)
        StopTmp=length(Fm);
    end    
    Fm(StartTmp:StopTmp)=mean(Fm(Starts(k):Stops(k)));
    
    StartTmp=Starts(k)-3;
    StopTmp=Starts(k)+6;
       if StartTmp<1
        StartTmp=1;
    end
    if StopTmp>length(Fm)
        StopTmp=length(Fm);
    end    
    Fm(StartTmp:StopTmp)=mean(Fm(Starts(k):Stops(k)));
end

Fm=smooth(Fm,15);

pos=find((Fm>0)&(Fm<1000));
Fm(pos)=1000;

Fm(pos0)=0;




function phase = wrap(phase)

    phase = phase - round(phase/2/pi)*2*pi;

    if phase>pi;        phase=phase-2*pi;
    elseif phase<-pi;   phase=phase+2*pi; end

return


function y = gaussmf_new(x, params)

if nargin ~= 2
    error('Two arguments are required by the Gaussian MF.');
elseif length(params) < 2
    error('The Gaussian MF needs at least two parameters.');
elseif params(1) == 0,
    error('The Gaussian MF needs a non-zero sigma.');
end

sigma = params(1); c = params(2);
y = (1/(sqrt(2*pi)*sigma))*exp(-(x - c).^2/(2*sigma^2));