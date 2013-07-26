% Short test script for code adding by John Kane on 28th June 2013

clear all;

% Load soundfile
[x,fs] = wavread('arctic_a0007.wav');

% Check the speech signal polarity
polarity = polarity_reskew(x,fs);
x=polarity*x;

% Extract the pitch and voicing information
[srh_f0,srh_vuv,srh_vuvc,srh_time] = pitch_srh(x,fs,80,500,10);

% Creaky probability estimation
[creak_pp,creak_bin] = detect_creaky_voice(x,fs); % Detect creaky voice
creak=interp1(creak_bin(:,2),creak_bin(:,1),1:length(x));
creak(creak<0.5)=0; creak(creak>=0.5)=1;

% GCI estimation
sd_gci = gci_sedreams(x,fs,median(srh_f0),1);        % SEDREAMS
se_gci = se_vq(x,fs,median(srh_f0),creak);           % SE-VQ

res = lpcresidual(x,25/1000*fs,5/1000*fs,fs/1000+2); % LP residual

m = mdq(res,fs,se_gci); % Maxima dispersion quotient measurement

ps = peakslope(x,fs);   % peakSlope extraction

gf_iaif = iaif_ola(x,fs);    % Glottal flow by the IAIF method

dgf_cc = complex_cepstrum(x,fs,sd_gci,srh_f0,srh_vuv); % Glottal flow derivative by the complex cesptrum-based method


% Plots
t=(0:length(x)-1)/fs;

fig(1) = subplot(511);
    plot(t,x, 'b');
    hold on
    plot(srh_time, srh_vuv, 'g');
    stem(sd_gci,ones(1,length(sd_gci))*-.1,'m');
    stem(se_gci,ones(1,length(se_gci))*-.1,'r');
    legend('Speech signal','Voicing (SRH)', 'GCI (SEDREAMS)','GCI (SE-VQ)');
    xlabel('Time [s]');
    ylabel('Amplitude');

fig(2) = subplot(512);
    plot(srh_time, srh_f0, 'r');
    legend('f0 (SRH)');
    xlabel('Time [s]');
    ylabel('Hz');

fig(3) = subplot(513);
    plot(t,x, 'b');
    hold on
    plot(ps(:,1), ps(:,2),'--b');
    plot(m(:,1),m(:,2),'m');
    legend('Speech signal','PeakSlope','MDQ');
    xlabel('Time [s]');

fig(4) = subplot(514);
    plot(t,x, 'b');
    hold on
    plot(t,norm(x)*gf_iaif./norm(gf_iaif), 'g');
    plot(creak_pp(:,2)/fs,creak_pp(:,1),'r')
    xlabel('Time [s]');
    ylabel('Amplitude');
    legend('Speech signal','Glottal flow (IAIF)','Creaky voice probability','Location','NorthWest');
    
    
fig(5) = subplot(515);
    plot(t,x, 'b');
    hold on
    plot(t,norm(x)*dgf_cc./norm(dgf_cc), 'r');    
    xlabel('Time [s]');
    ylabel('Amplitude');
    legend('Speech signal','Glottal flow derivative (CC)');
    

linkaxes(fig, 'x');
xlim([0 srh_time(end)]);
