function [peak_prom,peak_t,peak_idx] = get_res_peak_prom(rep,fs,winLen,x,Es)

% Function to find peak prominence values on a frame-by-frame basis.
% However, non-overlapping frames are used, and some shifting in position
% is done according to the maximum value
%
% REFERENCE:
%       Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%
%       Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', 27(4), pp. 1028-1047, Computer Speech and Language

%% Initial settings
Es=interp1(linspace(1,length(x),length(Es)),Es,1:length(x)); % Interpolated energy contour

ener_tol = 400/1000*fs;
ener_thresh=4;

peak_prom=zeros(1,round(length(rep)/winLen));
peak_t=zeros(1,round(length(rep)/winLen));
peak_idx=zeros(1,round(length(rep)/winLen));
winLenRem=winLen*.2; % just under half so periodic signals will have space for another peak

start=1;
stop=start+winLen-1;
n=1;
shift=20/1000*fs;

%% Do processing
while stop <= length(rep)
    
    frame_init = rep(start:stop);
    peak_t(n)=round(mean([start stop]));
    [~,maxIdx]=max(abs(frame_init));
    
    
    if maxIdx+start-round(winLen/2) < 1 || start+maxIdx+round(winLen/2) > length(rep)
        peak_prom(n)=0;
    else 
        start2=maxIdx+start-round(winLen/2);
        stop2=start+maxIdx+round(winLen/2);
        frame2=rep(start2:stop2);
        [~,maxIdx]=max(abs(frame2));
        if frame2(maxIdx)< 0 
            frame2=frame2*-1;
        end
        
        peak_idx(n) = maxIdx+start2-1;
        midpoint = round(winLen/2);
        frame2=frame2/max(frame2); % Normalise amplitude
        frame2(round(midpoint-winLenRem):round(midpoint+winLenRem))=0;
        frame2(frame2<0)=0;
        
        peak_prom(n)=1-max(frame2);
        
        if maxIdx+start-ener_tol < 1
            start_ener=1;
        else start_ener = maxIdx+start-ener_tol;
        end
        if maxIdx+start+ener_tol > length(x)
            stop_ener = length(x);
        else stop_ener = maxIdx+start+ener_tol; 
        end
        
        max_ener=max(Es(start_ener:stop_ener));
        
        if Es(maxIdx+start) < max_ener-ener_thresh
            peak_prom(n)=0;
        end      
    end
    
    start=start+shift;
    stop=start+winLen-1;
    n=n+1;
end

% Median filtering to remove outliers
peak_prom=medfilt1(peak_prom,5);

peak_prom(peak_t==0)=[];
peak_idx(peak_t==0)=[];
peak_t(peak_t==0)=[];
    