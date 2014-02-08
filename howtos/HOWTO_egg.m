% Demonstration script for processing ElectroGlottoGraphic (EGG) signals
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

clear all;

% Load sound and EGG files
[egg, fs] = wavread(['howtos' filesep 'DB-cresc-a-D4-m1.egg']);
[wav, fs] = wavread(['howtos' filesep 'DB-cresc-a-D4-m1.wav']);

degg=diff(egg);
degg(length(egg))=degg(length(egg)-1);

degg = degg./max(abs(degg));

[decom_Oqeggmess, decom_f0mess, decom_npicferms, decom_npicouvers, decom_atimes] = decom_sig(degg, fs);
idx = find(decom_f0mess==0);
decom_f0mess(idx) = NaN;
decom_Oqeggmess(idx) = NaN;


[howard_Oqeggmess, howard_f0mess, howard_atimes] = oq_egg_sig(egg, fs, [], 'howard');
idx = find(howard_f0mess==0);
howard_f0mess(idx) = NaN;
howard_Oqeggmess(idx) = NaN;

% Plots
t=(0:length(degg)-1)/fs;

fig(1) = subplot(411);
    plot(t, wav, 'k');
    hold on
    plot(t, 1+egg, 'b');
    plot(t, 3+degg, 'r');
    grid on;
    xlabel('Time [s]');
    ylabel('Amplitude');
    legend({'Waveform', 'EGG', 'DEGG'});

fig(2) = subplot(412);
    plot(howard_atimes, howard_f0mess, 'b');
    hold on
    plot(decom_atimes, decom_f0mess, 'r');
    grid on;
    ylim([200 400]);
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    legend({'f_0 (from DEGG)', 'f_0 (from EGG)'});

fig(3) = subplot(413);
    plot(decom_atimes, decom_Oqeggmess, 'b');
    hold on
    plot(howard_atimes, howard_Oqeggmess, 'r');
    grid on;
    ylim([0 0.8]);
    xlabel('Time [s]');
    ylabel('Oq');
    legend({'Estimated Oq using DECOM method', 'Estimated Oq using Howard method'});

fig(4) = subplot(414);
    plot(decom_atimes, decom_npicferms, 'b');
    hold on
    plot(decom_atimes, decom_npicouvers, 'r');
    grid on;
    ylim([0 5]);
    xlabel('Time [s]');
    ylabel('Number of peaks');
    legend({'Number of peaks at closure (DECOM method)', 'Number of peaks at opening (DECOM method)'});

linkaxes(fig, 'x');
