% Function to calculate the Intra-Frame Periodicity (IFP) contour) used as
% part of the creaky voice detection algorithm in Ishi et al (2008)
%
% Description
%  Function to calculate the Intra-Frame Periodicity (IFP) contour) used as
% part of the creaky voice detection algorithm in Ishi et al (2008)
%
%
% Inputs
%  x_filt   : [samples] [Nx1]  Bandlimited speech signal
%  fs       : [Hz]      [1x1]  Sampling frequency
%  IFPthresh: [samples] [1x1]  IFP threshold (default to 0.5 as per Ishi et al 2008)
%
% Outputs
%  IFP      : [samples] [Mx1] IFP contour
%  t_IFP    : [samples] [Mx1] Time locations of IFP contour
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

function [IFP,t_IFP,IFP_peaks] = get_ifp(x_filt,fs,IFPthresh,t_pow,peak_idx)

if nargin < 3
    IFPthresh=0.5;
end

%% Settings
N = round(32*(fs/1000)); % 32ms frame length for periodicity analysis
shift = round(10*(fs/1000)); % 10ms frame length for periodicity analysis
maxLag = round(15*(fs/1000)); % 15ms is maximum lag duration
peak_tol=1/1000*fs;
t_IFP = zeros(1,round((length(x_filt)-N)/shift));
IFP = zeros(1,round((length(x_filt)-N)/shift));

start=1;
finish=start+N-1;
n=1;

%% Do processing
while finish <= length(x_filt)

    frame = x_filt(start:finish); % select frame
    t_IFP(n)=round(mean([start finish]));
    [ACF,lags]=xcorr(frame,maxLag,'coeff'); % get normalised autocorrelation function
    ACF(lags<0)=[];
    lags(lags<0)=[];
    ACF(1:peak_tol)=0;
  %  plot(ACF)
   %
   ACF = (N./(N-(1:length(ACF))))'.*ACF(:);
   % hold on, plot(ACF,'r'), hold off, pause(.5)
    [ACF_max,max_idx]=max(ACF);
    
    m=2;
    while max_idx(1)*m < length(ACF)
        
        if (max_idx(1)*m)-peak_tol < 1
            begin=1;
        else begin=(max_idx(1)*m)-peak_tol;
        end
        if (max_idx(1)*m)+peak_tol > length(ACF)
            stop=length(ACF);
        else stop = (max_idx(1)*m)+peak_tol;
        end
        
        [ACF_max(m),max_idx(m)]=max(ACF(begin:stop));
        max_idx(m)=max_idx(m)+begin-1;
        m=m+1;
    end
    
    norm_fact = N./(N-max_idx); % noramlisation factor Eq q, Ishi et al (2008)
    
    if isempty(max_idx)
        IFP(n) = 0;
    else %IFP(n) = min(norm_fact.*ACF_max); % 
        IFP(n) = min(ACF_max); % 
    end
    
    start=start+shift;
    finish=start+N-1;
    n=n+1;
end

succFrames=3;

for n=1:length(IFP)-succFrames
    
    if IFP(n+1) < IFPthresh || IFP(n+2) < IFPthresh %|| IFP(n+3) < IFPthresh
        IFP(n)=0;
    end
end

if nargin > 3
    IFP_peaks=zeros(1,length(peak_idx));
    peak_idx_samp=round(t_pow(peak_idx));

    for n=1:length(IFP_peaks)

        [~,idx]=min(abs(t_IFP-peak_idx_samp(n)));
        IFP_peaks(n) = IFP(idx);
    end
else IFP_peaks=[];
end
    
