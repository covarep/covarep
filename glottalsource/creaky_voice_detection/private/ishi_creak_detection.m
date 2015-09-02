% Function to do creaky voice detection using the algorithm described in
% Ishi et al 2008
%
% Description
% Function to do creaky voice detection using the algorithm described in
% Ishi et al 2008
%
%
% Inputs
%  x        : [samples] [Nx1]  Speech signal
%  fs       : [Hz]      [1x1]  Sampling frequency
%  plots    : [integer] [1x1]  Option for plots (1==do plot)
%
% Outputs
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] Ishi, C., Sakakibara, K-I, Ishiguro, H., (2008) `A method for 
%       automatic detection of vocal fry', IEEE TASLP, 16(1), 47-56.
%  [2] Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%  [3] Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', Computer Speech and Language 27(4), pp.
%       1028-1047.
%  [4] Drugman, T., Kane, J., Gobl, C., (2012) `Resonator-based creaky 
%       voice detection', Interspeech 2012, Portland, Oregon, USA.
%
% Copyright (c) 2013 University of Mons, FNRS & 2013 Trinity College Dublin
%
% License
%  This code is a part of the GLOAT toolbox with the following
%  licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Authors
%  Thomas Drugman <thomas.drugman@umons.ac.be> & John Kane <kanejo@tcd.ie>

function [creak,t_IFP,IFP,PwP,IPS,t_pow,veryShort_powCont_dB,peak_idx,peak_pow,peak_idx_samp] = ...
    ishi_creak_detection(x,fs,plots)


%% Thresholds from Ishi et al (2008)
localpeak_thresh = 2; % 2 dB threshold, 
PwPthresh = 7; % Power peaks above 7 dB for creak
IFPthresh = 0.5; % Infraframe perdiocity below 0.5 for creak
IPSthresh = 0.5; % Interpulse similarity above 0.5 for creak
Tmax = 100*(fs/1000); % Maximum interpulse interval 100 ms (i.e. F0 = 10 Hz)

%% Bandpass filter signal 100 ~ 1500 Hz
Rp = 3; Rs = 40;
Wp = [100 1500]/(fs/2);
Ws = [50 3000]/(fs/2);
[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[b,a]=butter(n,Wn);
x_filt=filtfilt(b,a,x); % Zero-phase filtering

%% (1) Get PwP values from "very short-term" power contour
veryShort_len = 4*(fs/1000); % 4ms frame length for "very short-term" analysis
veryShort_shift = 2*(fs/1000); % 2ms shift for "very short-term" analysis
veryShort_powCont = zeros(1,ceil((length(x)-veryShort_len)/veryShort_shift));
t_pow = zeros(1,ceil((length(x)-veryShort_len)/veryShort_shift));
start=1;
finish=start+veryShort_len-1;

% Get very short term power contour
n=1;
x_filt2 = x_filt.^2;
while finish <= length(x_filt)
    veryShort_powCont(n) = mean(x_filt2(start:finish));
    t_pow(n)=mean([start finish]);
    start = start + veryShort_shift;
    finish=start+veryShort_len-1;
    n=n+1;
end
clear x_filt2;

veryShort_powCont_dB = 20*log10(veryShort_powCont);

% Find peaks in power contour
[peak_idx,peak_pow]=v_findpeaks(veryShort_powCont_dB);

for n=1:length(peak_pow)
    
    if peak_idx(n)-3 > 0
        start=peak_idx(n)-3;
    else start=1;
    end
    if peak_idx(n)+3<=length(veryShort_powCont_dB)
        finish=peak_idx(n)+3;
    else finish=length(veryShort_powCont_dB);
    end
    
    if peak_pow(n)-localpeak_thresh < veryShort_powCont_dB(start) || ...
            peak_pow(n)-localpeak_thresh < veryShort_powCont_dB(finish)
        peak_pow(n)=NaN;
        peak_idx(n)=NaN;
    end
end

peak_pow(isnan(peak_pow))=[];
peak_idx(isnan(peak_idx))=[];

% Suggestion from Ishi
maxPow=max(peak_pow);
max2remove_dB=30;
peak_idx(peak_pow<maxPow-max2remove_dB)=[];
peak_pow(peak_pow<maxPow-max2remove_dB)=[];

% Get power rising and power falling values
PwP.rise =zeros(1,length(peak_idx));
PwP.fall =zeros(1,length(peak_idx));
PwP.idx=peak_idx;

for n=1:length(peak_idx)
    if peak_idx(n)-5 > 0
        start = peak_idx(n)-5;
    else start = 1;
    end
    if peak_idx(n)+5 <= length(veryShort_powCont_dB)
        finish = peak_idx(n)+5;
    else finish = length(veryShort_powCont_dB);
    end
    
    % power rising
    PwP.rise(n) = max(peak_pow(n)- ...
        veryShort_powCont_dB(start:peak_idx(n)-1));
    % power falling
    PwP.fall(n) = max(peak_pow(n)- ...
        veryShort_powCont_dB(peak_idx(n)+1:finish));
end
    
%% (2) Calculate Intraframe Periodicity Measure (IFP)
peak_idx_samp=round(t_pow(peak_idx));
[IFP,t_IFP,IFP_peaks] = get_ifp(x_filt,fs,IFPthresh,t_pow,peak_idx);

%% (3) Calculate interpulse similarity measure (IPS)
IPS = zeros(1,length(peak_idx));
N = 15*(fs/1000); % 15 ms frame length for similarity analysis
maxLags = 15*(fs/1000); % 15 ms frame length for similarity analysis

for n=2:length(IPS);
    if peak_idx_samp(n)-peak_idx_samp(n-1) > Tmax
        IPS(n) = 0;
    else

        if peak_idx_samp(n)+round(N/2)>length(x)
            break
        end
        if peak_idx_samp(n-1)-round(N/2) < 1 || peak_idx_samp(n)+round(N/2) > length(x_filt)
            IPS(n) = 0;
        else
            peak_seg_cur = x_filt(peak_idx_samp(n)-round(N/2):peak_idx_samp(n)+round(N/2));
            peak_seg_prev = x_filt(peak_idx_samp(n-1)-round(N/2):peak_idx_samp(n-1)+round(N/2));
            Rxx = xcorr(peak_seg_cur,peak_seg_prev,maxLags,'coeff'); % Normalised cross correlation
            IPS(n)=max(Rxx);
        end
    end
end

%% Decide on candidate peaks based on set threholds
creak_peak_idx=zeros(1,length(peak_idx_samp));
creak_peak_idx(PwP.rise>=PwPthresh & PwP.fall>=PwPthresh & ...
    IFP_peaks<=IFPthresh & IPS>=IPSthresh)=1; % binary vector
creak=zeros(1,length(x));
creak_dur_ms=[];
n=1;
cnt=1;
while n < length(creak_peak_idx)
    if creak_peak_idx(n)==1 && creak_peak_idx(n+1)==1
        start = peak_idx_samp(n);
        n=n+1;
        while creak_peak_idx(n)==1 && peak_idx_samp(n) ...
                - peak_idx_samp(n-1) < Tmax && n < length(creak_peak_idx)
            n=n+1;
        end
        finish= peak_idx_samp(n-1);
        
        creak(start:finish)=1;
        cnt=cnt+1;

    end
    n=n+1;
end

%creak=creak.*VUV; % Remove creak detected in unvoiced areas
IFP(IFP<0)=0;
IFP(IFP>1)=1;
%% Plots
if nargin > 2
    if plots==1
        t=(1:length(x))/fs;
       % figure,
       x=x/max(x);
        subplot(411), plot(t,x/max(x)), hold on, plot(t,(creak)*.8,'r'), hold off
         axis tight, ylabel('Amplitude','FontSize',14)
         title('Speech waveform with creak regions','FontSize',14,'FontWeight','bold'), ylim([min(x)*1.1 max(x)*1.1])
        subplot(412), plot(t(round(t_pow)),veryShort_powCont_dB), hold on, 
         plot(t(peak_idx_samp),peak_pow,'rx'), 
         plot(t(peak_idx_samp(PwP.rise>PwPthresh)),peak_pow(PwP.rise>PwPthresh),'gx'),
         hold off, axis tight, 
         ylabel('Power (dB)','FontSize',14)
         title('Very short power contour (dB) with Local peaks','FontSize',14,'FontWeight','bold')
        subplot(413), plot(t(t_IFP),IFP),hold on,% plot([0 length(x)], ...
           % [IFPthresh IFPthresh],'k'), plot(t_IFP,IFP_orig,'r')
            hold off, axis tight, title('Intraframe periodicity','FontSize',14,'FontWeight','bold')
            ylabel('IFP','FontSize',14)
        subplot(414), stem(t(peak_idx_samp),IPS), hold on,  plot(t([1 length(x)]), ...
            [IPSthresh IPSthresh],'k'), hold off, axis tight, 
        ylabel('IPS','FontSize',14),xlabel('Time (seconds)','FontSize',14)
        title('Interpulse similarity','FontSize',14,'FontWeight','bold')
    end
end