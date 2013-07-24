function [PwP,IFP,IPS,dec_cur,dec_orig,IPS_cur,time] = get_ishi_params_inter(x,fs)

% Function to derive acoustics features used in Ishi et al. (2008). The
% glottal synchronous parameters (i.e. PwP and IPS) are resampled to a
% fixed update rate.
%
% REFERENCE:
%       Ishi, C., Sakakibara, K-I, Ishiguro, H., (2008) `A method for 
%       automatic detection of vocal fry', IEEE TASLP, 16(1), 47-56.

%% Initial settings
% Thresholds from Ishi et al (2008)
PwP_thresh=7;
IFP_thresh=0.5;
IPS_thresh=0.5;
maxLen=35/1000*fs;

% Allocate space
IPS=zeros(1,length(x));
PwP.rise=zeros(1,length(x));
PwP.fall=zeros(1,length(x));

%% Extract parameters
[dec_orig,t_IFP,IFP_cur,PwP_cur,IPS_cur,t_pow] = ishi_creak_detection(x,fs,0);
time=[];
%% Resample
if isempty(IFP_cur) || length(IFP_cur) < 3
    IFP=zeros(1,length(x));
else
    IFP=interp1(t_IFP,IFP_cur,1:length(x));
    IFP(isnan(IFP))=0;
end

if isempty(PwP_cur.rise)==0 && length(PwP_cur.rise) > 2
    
    time=round(t_pow(PwP_cur.idx));
    
    for n=1:length(x)
        
        % Find nearest value
        time_cur=time;
        time_cur(time_cur<n)=0;
        [minDist,idx]=min(abs(time_cur-n));
        
        if minDist < maxLen
            PwP.rise(n)=PwP_cur.rise(idx);
            PwP.fall(n)=PwP_cur.fall(idx);
            IPS(n)=IPS_cur(idx);
        end
    end
           
end

%% Generate binary decision
dec_cur=zeros(1,length(x));
dec_cur(PwP.rise>PwP_thresh&IFP<=IFP_thresh&IPS>=IPS_thresh)=1;
        
        