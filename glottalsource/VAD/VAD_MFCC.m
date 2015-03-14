function [MFCC] = VAD_MFCC(wave,Fs)

%% Extracts MFCC features from the speech signal

frameL=30;
frameshift=10;

Nfft=1024;

Start=1;
Stop=round(frameL/1000*Fs);
Shift=round(frameshift/1000*Fs);
Ind=1;
Win=blackman(Stop);
Es=[];

wave=filter([1 -0.97],1,wave);

while Stop<length(wave)
    Seg=wave(Start:Stop);         
    Seg=Seg.*Win;
    
    Spec=fft(Seg,Nfft);
    
    mfcc = hspec2fwcep(abs(Spec(1:end/2+1)), Fs, 12);

    MFCC(Ind,:)=mfcc;
        
    Start=Start+Shift;
    Stop=Stop+Shift;
    Ind=Ind+1;    
end

for i = 1 : 13
    MFCC( :, i ) = MFCC( : , i ) - mean(MFCC( : , i));
end

MFCC=MFCC';