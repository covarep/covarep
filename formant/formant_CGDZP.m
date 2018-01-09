% This script tracks the first 5 formants based on the chirp group delay.
%
% Octave compatible
%
% Description
% This formant tracking technique is fully described in [1]. The algorithm
% is based on processing the negative derivative of the argument of the
% chirp-z transform (termed as the differential phase spectrum) of a given
% speech signal. No modeling is included in the procedure but only peak
% picking on differential phase spectrum.
% CGDZP = Chirp Group Delay Zero Phase
%
%
% [formantPeaks,t_analysis]=formant_CGDZP(wave,fs,frameSize,frameShift)
%
% Inputs
%  wave             : [samples] [Nx1] input signal (speech signal)
%  fs               : [Hz]      [1x1] sampling frequency
%  frameSize        : [ms]      [1x1] size of the frames (default = 30)
%  frameShift       : [ms]      [1x1] duration between two successive
%                               frames. (default = 10)
%
% Outputs
%  formantPeaks    : vector containing the estimated frequencies (in Hz) of
%  the first 5 formants.
%  t_analysis      : [s] vector containing the analysis time instants.
%
% References
% [1] B. BOZKURT, B. DOVAL, C. D'ALESSANDRO, T. DUTOIT, 2004, "Improved
% differential phase spectrum processing for formant tracking",
% Proc. Icslp 2004, Jeju Island (Korea).
%
%  Publication available at the following link:
%  http://tcts.fpms.ac.be/publications/papers/2004/icslp2004_bbdcdtd3.pdf
%
% Copyright (c) 2011 University of Mons, FNRS
%
% License
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Baris Bozkurt / Thomas Drugman thomas.drugman@umons.ac.be
%

function [formantPeaks,t_analysis]=formant_CGDZP(wave,fs,frameSize,frameShift)

if nargin<4
    frameShift=10;
end

if nargin<3
    frameSize=30;
end
    
numFormants=5;
numFormatsFinal=numFormants;
frameSize=round(fs/1000*frameSize);
frameShift=round(fs/1000*frameShift);
fsLR=2048;%lower resolution
viewRange=round(fsLR/3.2);
R=1.12;Rfix=R;
maxFormantDelta=250;%discontinuity constraint(in Hz.) in formant tracks 

n=0:frameSize-2;%zeroPhase data has length of N-1

ScalingFreqAxis=0.5:1/(viewRange-1):1.5;
ScalingFreqAxis=ScalingFreqAxis-0.25;
CepsSmoothSelectPercent=3;
SelectTimeIndex=round(CepsSmoothSelectPercent*fsLR/2/100);

SIZEwave=length(wave);
numFrames=floor((SIZEwave-frameSize)/frameShift);

t_analysis=zeros(1,numFrames);
blackWin = blackman(frameSize);

formantPeaks=[];
for kk=0:numFrames-1
    
    speechData=wave(kk*frameShift+1:kk*frameShift+frameSize);
    windowedData=speechData.*blackWin;
    
    zeroPhaseData=real(ifft(abs(fft(diff(windowedData)))));%obtain zero-pha version
    
    numPeaks=0;R=Rfix;%the following loop searches for the R value where we have numFormants number of formats
    while(numPeaks~=numFormants && R>1.01 && R<1.25)
        %chirp z-transform calculation using fft,...multiplication with an exponential function is sufficient
        exponentialEnvelope=exp(log(1/R)*n)';%this is needed for computation of z-transform using fft
        fourierTrans=fft(zeroPhaseData.*exponentialEnvelope,fsLR);
        angFFT=angle(fourierTrans(1:viewRange));
        chirpGroupDelay=-diff(angFFT)';
        %chirpGroupDelay=chirpGroupDelay-sum(chirpGroupDelay)/viewRange;
        
        [peakIndex]=formantPeakPick(chirpGroupDelay,1);
        numPeaks=length(peakIndex);
        if(numPeaks>numFormants && R>=Rfix)
            R=R+0.01;peakIndex=peakIndex(1:numFormants);
        elseif(numPeaks<numFormants && R<=Rfix)
            R=R-0.01;peakIndex=[peakIndex zeros(1,numFormants-numPeaks)];
        else
            break;
        end
    end 
    
    %chirpGroupDelay_spectrogram(kk+1,:)=chirpGroupDelay;
    
    if(numPeaks>numFormants)
        peakIndex=peakIndex(1:numFormants);
    elseif(numPeaks<numFormants)
        peakIndex=[peakIndex zeros(1,numFormants-size(peakIndex,2))];
    end
    formantPeaks=[formantPeaks ; peakIndex];
    
    t_analysis(kk+1)=round((kk*frameShift+1+kk*frameShift+frameSize)/2)/fs;
    
end
formantPeaks=round(formantPeaks*fs/fsLR);

if(1)%filtering
    %form a matrix of distance cost
    formantPeaksCost=formantPeaks*0;
    for kk=3:numFrames-2
        prePrePeaks=formantPeaks(kk-2,:);postPostPeaks=formantPeaks(kk+2,:);
        prePeaks=formantPeaks(kk-1,:);postPeaks=formantPeaks(kk+1,:);
        currentPeaks=formantPeaks(kk,:);
        currentPeaksCost=currentPeaks*0;%cost will correspond to average distance from closest match to neighbours
        for mm=1:numFormants
            if(currentPeaks(mm)==0)
                currentPeaksCost(mm)=fs/2; 
            else
                %search for closest matches
                distanceArrayPre=sort(abs(prePeaks-currentPeaks(mm)));
                distanceArrayPost=sort(abs(postPeaks-currentPeaks(mm)));
                distanceArrayPrePre=sort(abs(prePrePeaks-currentPeaks(mm)));
                distanceArrayPostPost=sort(abs(postPostPeaks-currentPeaks(mm)));
                allDistances=[(distanceArrayPre(1)+distanceArrayPost(1))/2 distanceArrayPre(1) distanceArrayPost(1)];
                allDistances2=[(distanceArrayPrePre(1)+distanceArrayPostPost(1))/2 distanceArrayPrePre(1) distanceArrayPostPost(1)];
                currentPeaksCost(mm)=min([allDistances allDistances2]);
            end
        end
        formantPeaksCost(kk,:)=currentPeaksCost;
    end
    %make decision based on costs
    for kk=1:numFrames
        for mm=1:numFormants
            if(formantPeaksCost(kk,mm)>maxFormantDelta)
                formantPeaks(kk,mm)=0;
            end
        end
    end
end

if(1)%replace possible continuation values instead of zero values
    for kk=2:numFrames-1
        currentPeaks=formantPeaks(kk,:);
        [~,index]=find(currentPeaks);%finds non-zero elements
        nonZeroFormants=currentPeaks(index);numNonZeroFormants=length(index);
        numZeroFormants=numFormants-numNonZeroFormants;
        if(numNonZeroFormants<numFormants && ~isempty(index))
            %smoothArray=round(sort(formantPeaks(kk-1,:)+formantPeaks(kk+1,:))/2);
            possibleValues=sort([formantPeaks(kk-1,:) formantPeaks(kk+1,:)]);
            while((length(possibleValues)>0) && (possibleValues(1)==0)) %discard zero entries
                possibleValues=possibleValues(2:length(possibleValues));
            end
            possibleCandidates=[];%candidate for the zero valued formant
            for mm=1:length(possibleValues)
                distanceArray=sort(abs(nonZeroFormants-possibleValues(mm)));
                if(isempty(distanceArray))
                    possibleCandidates=[possibleCandidates possibleValues(mm)];
                elseif(distanceArray(1)>maxFormantDelta)%this possible value not found in the current vector
                    possibleCandidates=[possibleCandidates possibleValues(mm)];
                end
            end
            %choose among possible candidates
            lenPossibleCandidates=length(possibleCandidates);
            if(lenPossibleCandidates<=numZeroFormants)
                currentPeaks=sort([nonZeroFormants possibleCandidates zeros(1,numZeroFormants-lenPossibleCandidates)]);
            elseif(numZeroFormants==1)%the most common case
                [~,index]=sort(diff(possibleCandidates));
                currentPeaks=sort([nonZeroFormants possibleCandidates(index(1))]);
            elseif(numZeroFormants<lenPossibleCandidates)
                possibleCandidates=possibleCandidates(1:numZeroFormants);
                currentPeaks=sort([nonZeroFormants possibleCandidates]);
            end
            formantPeaks(kk,:)=currentPeaks;
        end
    end
end






%[peakIndex]=formantPeakPick(diffPhase,minPeakDist)
%
%a simple peak picking algorithm
%checks if the value of an entry is bigger than all neighbour values
%6 values on each side is controlled

function [peakIndex]=formantPeakPick(diffPhase,minPeakDist)

lendiffPhase=length(diffPhase);

peakIndex=[];
for kk=6:lendiffPhase-6
    if (diffPhase(kk)>=diffPhase(kk-1) && diffPhase(kk)>=diffPhase(kk+1)) ...%first neighbours
            && (diffPhase(kk)>=diffPhase(kk-2) && diffPhase(kk)>=diffPhase(kk+2)) ...%second neighbours
            && (diffPhase(kk)>diffPhase(kk-3) && diffPhase(kk)>diffPhase(kk+3)) ...%third neighbours
            && (diffPhase(kk)>diffPhase(kk-4) && diffPhase(kk)>diffPhase(kk+4)) ...%fourth neighbours
            && (diffPhase(kk)>diffPhase(kk-5) && diffPhase(kk)>diffPhase(kk+5))%fifth neighbours
        peakIndex=[peakIndex kk];
    end
end

lenPeakInd=length(peakIndex);
kk=2;
while(kk<=lenPeakInd)
    if(peakIndex(kk)-peakIndex(kk-1)<minPeakDist)
        
        peakIndex(kk)=round((peakIndex(kk)*diffPhase(kk)+peakIndex(kk-1)*diffPhase(kk-1))/(diffPhase(kk)+diffPhase(kk-1)));
        peakIndex=[peakIndex(1:kk-2) peakIndex(kk:lenPeakInd)];
        kk=kk-1;lenPeakInd=lenPeakInd-1;
    end
    kk=kk+1;
end
