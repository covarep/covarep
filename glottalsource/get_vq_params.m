% Function to estimate the glottal parameters: NAQ, QOQ, H1-H2, HRF and PSP
%
% Octave compatible
%
% Description
%  This function can be used to estimate a range of conventional glottal
%  source parameters often used in the literature. This includes: the
%  normalised amplitude quotient (NAQ), the quasi-open quotient (QOQ), the
%  difference in amplitude of the first two harmonics of the differentiated
%  glottal source spectrum (H1-H2), the harmonic richness factor (HRF) and
%  the parabolic spectral parameter (PSP)
%
% Inputs
%  gf       : [samples] [Nx1] Glottal flow estimation
%  gfd      : [samples] [Nx1] Glottal flow derivative estimation
%  fs       : [Hz]      [1x1] sampling frequency
%  GCI      : [s]       [Mx1] Glottal closure instants 
%
% Outputs
%  NAQ      : [s,samples] [Mx2] Normalised amplitude quotient
%  QOQ      : [s,samples] [Mx2] Quasi-open quotient
%  H1H2     : [s,dB]      [Mx2] Difference in glottal harmonic amplitude
%  HRF      : [s,samples] [Mx2] Harmonic richness factor
%  PSP      : [s,samples] [Mx2] Parabolic spectral parameter
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] Alku, P., B ackstrom, T., and Vilkman, E. Normalized amplitude quotient 
%     for parameterization of the glottal flow. Journal of the Acoustical 
%     Society of America, 112(2):701–710, 2002.
%  [2] Hacki, T. Klassifizierung von glottisdysfunktionen mit hilfe der 
%     elektroglottographie. Folia Phoniatrica, pages 43–48, 1989.
%  [3] Alku, P., Strik, H., and Vilkman, E. Parabolic spectral parameter - 
%     A new method for quantification of the glottal flow. Speech 
%     Communication, 22(1):67–79, 1997.
%  [4] Hanson, H. M. Glottal characteristics of female speakers: Acoustic 
%     correlates. Journal of the Acoustical Society of America, 
%     10(1):466–481, 1997.
%  [5] Childers, D. G. and Lee, C. K. Voice quality factors: Analysis, 
%     synthesis and perception. Journal of the Acoustical Society of 
%     America, 90(5):2394–2410, 1991.
%
% Copyright (c) 2013 Trinity College Dublin
%
% License
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  John Kane kanejo@tcd.ie

function [NAQ,QOQ,H1H2,HRF,PSP] = get_vq_params(gf,gfd,fs,GCI)

GCI = round(GCI*fs)+1;

%% Initial settings
F0min=20;
F0max=500;
NAQ=zeros(1,length(GCI));
QOQ=zeros(1,length(GCI));
H1H2=zeros(1,length(GCI));
HRF=zeros(1,length(GCI));
PSP=zeros(1,length(GCI));

glot_shift=round(0.5/1000*fs);
qoq_level=0.5; % Threshold for QOQ estimation
T0_num=3; % Number of local glottal pulses to be used for harmonic spectrum
min_harm_num=5;
HRF_freq_max=5000; % Maximum frequency used for harmonic measurement
PSP_fft_size=2048; % As per Alku et al (1997)

%% Do processing
for n=1:length(GCI)
    
    % Get glottal pulse compensated for zero-line drift
    if n==1
        start=1;
        stop=GCI(n);
        T0=GCI(n+1)-GCI(n);
    else start=GCI(n-1);
        stop=GCI(n);
        T0=GCI(n)-GCI(n-1);
    end
    F0=fs/T0; 
    
    if isinf(F0)==0 && T0~=0 && F0 > F0min && F0<F0max 
       
        % Compensate for zero-line drift in the glottal flow pulse
        gf_comb=[gf(start) gf(stop)];
        if start~=stop
            if length(gf_comb)>1
                line=interp1(linspace(1,stop-start+1,2),gf_comb, ...
                    1:stop-start+1);
            else line=gf_comb;
            end
        else line=0;
        end
        gf_seg=gf(start:stop);
        gf_seg_comp=gf_seg(:)-line(:);
        
        if stop+glot_shift <= length(gfd)
            stop2=stop+glot_shift;
        else
            stop2=stop;
        end
        gfd_seg=gfd(start:stop2);
        
        % Get NAQ and QOQ 
        d_peak = max(abs(gfd_seg));
        [f_ac,max_idx] = max(gf_seg_comp);
        Amid=f_ac*qoq_level;
        [T1,T2] = findAmid_t(gf_seg_comp,Amid,max_idx);
        
        NAQ(n) = (f_ac/d_peak)/T0;
        QOQ(n)=(T2-T1)/(fs/F0);
        
         % Get frame positions for H1-H2 parameter
        if GCI(n)-round((T0*T0_num)/2) > 0
            f_start = GCI(n)-round((T0*T0_num)/2);
        else f_start=1;
        end
        if GCI(n)+round((T0*T0_num)/2) <= length(gfd)
            f_stop = GCI(n)+round((T0*T0_num)/2);
        else f_stop=length(gfd);
        end
        f_frame=gfd(f_start:f_stop);
        f_frame=f_frame.*(2^15);
        f_win=f_frame(:).*hamming(length(f_frame));
        f_spec = mag2db(abs(fft(f_win,fs)));
        
         % Get H1-H2/HRF measurements
        [h_idx,h_amp]=v_findpeaks(f_spec,[],F0/2);
        HRF_harm_num=floor(HRF_freq_max/F0);
        if length(h_idx) >= min_harm_num
            [~,f0_idx] = min(abs(bsxfun(@minus,h_idx,(1:HRF_harm_num)*F0)),[],1);
            H1H2(n)=h_amp(f0_idx(1))-h_amp(f0_idx(2));
            HRF(n) = sum(h_amp(f0_idx(2:end)))/h_amp(f0_idx(1));
        end
        
        % Get PSP
        if n>1
           start=GCI(n-1);
           stop=GCI(n);
           gf_seg=gf(start:stop);
           gf_seg=gf_seg-min(gf_seg); % Adjust minimum to zero
           gf_seg=gf_seg/max(gf_seg); % Scale to unity
           gf_seg(PSP_fft_size)=0; % Zero-pad
           X = mag2db(abs(fft(gf_seg))); % Spectrum
           X = X(linspace(1,fs,length(X))<=fs/2); % Use up to the nyquist
           a = psp_get_a(X);
           
           % Estimate measure of maximal spectral decay
           a_max_sig=1/round(fs/F0)*ones(1,round(fs/F0));
           a_max_sig(PSP_fft_size)=0; % Zero-pad
           X_a_max = mag2db(abs(fft(a_max_sig))); % Spectrum
           X_a_max = X_a_max(linspace(1,fs,length(X_a_max))<=fs/2); % Use up to the nyquist
           a_max = psp_get_a(X_a_max);
           
           PSP(n)=a/a_max;
        end
        
    end
   
    
end

%% Add in time to parameters
NAQ=[GCI(:)/fs NAQ(:)];
QOQ=[GCI(:)/fs QOQ(:)];
H1H2=[GCI(:)/fs H1H2(:)];
HRF=[GCI(:)/fs HRF(:)];
PSP=[GCI(:)/fs PSP(:)];


function [T1,T2] = findAmid_t(glot_adj,Amid,Tz)

% Function to find the start and stop positions of the quasi-open phase.

T1=0;
T2=0;
if Tz~=0
    n=Tz;

    while glot_adj(n) > Amid && n > 3
        n=n-1;
    end
    T1=n;
    n=Tz;

    while glot_adj(n) > Amid && n < length(glot_adj)-2
        n=n+1;
    end
    T2=n;
end

function a = psp_get_a(X)

% Function to calculate a coefficient, as per section 3.2.1 in Alku et al
% (1997)

%% Initial settings
NE_min = 0.01;
flag = 0;
N = 3; % An set to an initial value of 3, as per Alku et al (1997)
X=X(:);

% pre-allocate&compute
X2=X.^2;
k2 = ((0:100).^2)';

% init sum
k2_sum=sum(k2(1:N-1));
k4_sum=sum(k2(1:N-1).^2);
X_k1_sum = sum(X(1:N-1));
X_k1_k2_sum = sum(X((1:N-1)).*k2(1:N-1));
X2_k1_sum = sum(X2((1:N-1)));

%% Do processing
while flag == 0

    % re-allocate if necessary
    if N>numel(k2)
        k2 = ((0:numel(k2)*2).^2)';
    end
  
    % iteratively built sum
    k2_sum= k2_sum + k2(N);
    k4_sum= k4_sum + k2(N)*k2(N);
    X_k1_sum = X_k1_sum + X(N);
    X_k1_k2_sum = X_k1_k2_sum + X(N)*k2(N);
    X2_k1_sum = X2_k1_sum + X2(N);
    
    % Equation 5 
    a = (N*X_k1_k2_sum-X_k1_sum.*k2_sum)/(N*k4_sum-k2_sum.^2); 
    
    X_k1_a_k2 = X(1:N)-a*k2(1:N);
    
    % Equation 4
    b = 1/N*sum(X_k1_a_k2);
    NE = (sum((X_k1_a_k2-b).^2))/X2_k1_sum;

    % Increment or stop
    if NE < NE_min
        N = N+1;
    else
        flag = 1;
    end    
end
