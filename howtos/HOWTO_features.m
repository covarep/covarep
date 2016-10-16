% Extract a few features from a single sample
%
% alpha and gamma, which are set to 0.4 and 0.9, need to be tuned for better perfomance on your own data
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%

% Using modified_group_delay_feature.m -----------------------------------------

file_name = '0011.arctic_bdl1.wav';

[speech, fs]  = wavread(file_name);

[grp_phase, cep, ts] = modified_group_delay_feature(file_name, 0.4, 0.9, 12);

figure;
fig(1) = subplot(211);
    times = (0:length(speech)-1)/fs;
    plot(times, speech, 'k');
    ylim(0.6*[-1, 1]);
    xlabel('Time [s]');
    ylabel('Amplitude');

fig(2) = subplot(212);
    F = fs*(0:256)/512;
    imagesc(ts, F, grp_phase, [0.0, 0.05]);
    axis xy;
    title('Modified Group Delay');
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');

linkaxes(fig, 'x');
xlim([times(1), times(end)]);
