function [] = HOWTO_formant()

clear all;

% Load soundfile
[x,fs] = wavread('arctic_a0007.wav');

[formantPeaks,t_analysis]=formant_CGDZP(x,fs,30,10);

t=(0:length(x)-1)/fs;
fig(1) = subplot(211);
    plot(t,x, 'b');
    legend('Speech signal');
    xlabel('Time [s]');
    ylabel('Amplitude');

    fig(2) = subplot(212);
    plot(t_analysis,formantPeaks(:,1), 'b');
    hold on
    plot(t_analysis,formantPeaks(:,2), 'r');
    plot(t_analysis,formantPeaks(:,3), 'g');
    plot(t_analysis,formantPeaks(:,4), 'm');
    plot(t_analysis,formantPeaks(:,5), 'c');
    legend('First formant','Second formant','Third formant','Fourth formant','Fifth formant');
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');


linkaxes(fig, 'x');
xlim([0 t(end)]);