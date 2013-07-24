function [Ts,Es,ZCs_ms,Xpos,pos] = sil_unv_features(wave,Fs,frame_len_ms)


wave=wave.*(2^15);

if nargin < 3
    frame_len_ms=10;
end
frame_len=frame_len_ms/1000*Fs;
frame_shift_ms=5;
Shift=frame_shift_ms/1000*Fs;
 
Start=1; % intialise
%Stop=10/1000*Fs; 
Stop=Start+frame_len-1;
%Shift=5/1000*Fs;
pos=mean([Start Stop]);

[ZCs,Xpos] = get_zero_x_rate(wave,frame_len,Shift);
ZCs_ms=ZCs/1000*Fs;

Ind=1;
 
Es=[];
Ts=[];
 
while Stop<length(wave)
   
    Mid=round((Start+Stop)/2);
    Ts(Ind)=Mid;
   
    Sig=wave(Start:Stop); % Select frame segment
    Sig=Sig(:).*hanning(length(Sig)); % Window segment

   
    Es(Ind)=mean(Sig.^2); % Get energy value
   
    Start=Start+Shift;
    Stop=Stop+Shift;
    pos(Ind)=mean([Start Stop]);
    Ind=Ind+1;
end 
 
Es=log(Es);