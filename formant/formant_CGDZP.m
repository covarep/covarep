% This script tracks the first 5 formants based on the chirp group delay.
%
% Octave compatible
%
% Description
%
% This formant tracking technique is fully described in [1]. The algorithm
% is based on processing the negative derivative of the argument of the
% chirp-z transform (termed as the differential phase spectrum) of a given
% speech signal. No modeling is included in the procedure but only peak
% picking on differential phase spectrum.
%
%
%  %
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
%
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

function [formantPeaks]=formant_CGDZP(wave,fs,frameSize,frameShift)

if nargin<4
    frameShift=10;
end

if nargin<3
    frameSize=30;
end
    
numFormants=5;
numFormatsFinal=numFormants;
frameSize=fs/1000*frameSize;
frameShift=fs/1000*frameShift;
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

formantPeaks=[];
for kk=0:numFrames-1
    
    speechData=wave(kk*frameShift+1:kk*frameShift+frameSize);
    windowedData=speechData.*blackman(length(speechData));
    
    numPeaks=0;R=Rfix;%the following loop searches for the R value where we have numFormants number of formats
    while(numPeaks~=numFormants & R>1.01 & R<1.25)
        zeroPhaseData=real(ifft(abs(fft(diff(windowedData)))));%obtain zero-pha version
        %chirp z-transform calculation using fft,...multiplication with an exponential function is sufficient
        exponentialEnvelope=exp(log(1/R)*n)';%this is needed for computation of z-transform using fft
        fourierTrans=fft(zeroPhaseData.*exponentialEnvelope,fsLR);
        angFFT=angle(fourierTrans(1:viewRange));
        chirpGroupDelay=-diff(angFFT)';
        %chirpGroupDelay=chirpGroupDelay-sum(chirpGroupDelay)/viewRange;
        
        [peakIndex]=formantPeakPick(chirpGroupDelay,1);
        numPeaks=length(peakIndex);
        if(numPeaks>numFormants & R>=Rfix)
            R=R+0.01;peakIndex=peakIndex(1:numFormants);
        elseif(numPeaks<numFormants & R<=Rfix)
            R=R-0.01;peakIndex=[peakIndex zeros(1,numFormants-numPeaks)];
        else
            break;
        end
    end 
    
    chirpGroupDelay_spectrogram(kk+1,:)=chirpGroupDelay;
    
    if(numPeaks>numFormants)
        peakIndex=peakIndex(1:numFormants);
    elseif(numPeaks<numFormants)
        peakIndex=[peakIndex zeros(1,numFormants-numPeaks)];
    end
    formantPeaks=[formantPeaks ; peakIndex];
    
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
        [val,index]=find(currentPeaks);%finds non-zero elements
        nonZeroFormants=currentPeaks(index);numNonZeroFormants=length(index);
        numZeroFormants=numFormants-numNonZeroFormants;
        if(numNonZeroFormants<numFormants)
            %smoothArray=round(sort(formantPeaks(kk-1,:)+formantPeaks(kk+1,:))/2);
            possibleValues=sort([formantPeaks(kk-1,:) formantPeaks(kk+1,:)]);
            while(possibleValues(1)==0)%discard zero entries
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
                [val,index]=sort(diff(possibleCandidates));
                currentPeaks=sort([nonZeroFormants possibleCandidates(index(1))]);
            elseif(numZeroFormants<lenPossibleCandidates)
                possibleCandidates=possibleCandidates(1:numZeroFormants);
                currentPeaks=sort([nonZeroFormants possibleCandidates]);
            end
            formantPeaks(kk,:)=currentPeaks;
        end
    end
end
