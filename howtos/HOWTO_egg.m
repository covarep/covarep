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
% Authors
%  Alexis Michaud <alexis.michaud@vjf.cnrs.fr> <michaud.cnrs@gmail.com>
%                  CNRS (Centre National de la Recherche Scientifique, France)
%  Gilles Degottex <gilles.degottex@gmail.com>
%                  Ircam, France
%

clear all;


% Example for the peakdet.m method ---------------------------------------------

% <peakdet.m> is intended for the semi-automatic analysis of the EGG signal for
% a set of continuously voiced items: for example, vowels, syllable rhymes, or
% sustained voiced sounds containing up to <MaxPerN> glottal cycles. 
%
% This demo shows the result for an example of type of signals for which peakdet
% was devised: a Vietnamese syllable with a glottalized ending.
%
% As of 2014, a guide to the use of <peakdet> is available from:
% http://voiceresearch.free.fr/egg/softwares.htm#peakdet
%

% reading the example file (electroglottographic signal)
[egg, fs] = audioread('glottalized-m1.egg');
[wav, fs] = audioread('glottalized-m1.wav');

% running the <peakdet.m> function
[results_matrix, SdSIG, SIG] = peakdet(egg, fs, 200, 3);

% Plotting the results
figure
% retrieving the first and last glottis-closure-instants detected by <peakdet>
firstclo = results_matrix(1,1);
lastclo = results_matrix(length(nonzeros(results_matrix(:,2))),2);

% EGG and dEGG in one figure, with indications on first and last detected
% Glottis-Closure-Instants (GCI)
fig(1) = subplot(311);
    times = (0:length(wav)-1)/fs;
    plot(times, wav, 'k');
    hold on
    plot(times, 1+SIG, 'b');
    times = (0:length(SdSIG)-1)/fs;
    plot(times, 2+SdSIG*30, 'r');
    plot([firstclo firstclo],[-1 3],'-r')
    plot([lastclo lastclo],[-1 3],'-r')
    xlabel('Time [s]');
    grid on
    legend({'EGG', 'DEGG', 'First and last GCI'});

% Second figure: estimated fundamental frequency (F0) and open quotient values
fig(2) = subplot(312);
    times = 0.5*(results_matrix(:,1)+results_matrix(:,2));
    plot(times, results_matrix(:,3), '-pb');
    grid on;
    ylim([0 250]);
    xlabel('Time [s]');
    ylabel('F0 [Hz]');
    legend({'Estimated F0'});

fig(3) = subplot(313);
    plot(times, results_matrix(:,5), '*g');
    hold on
    plot(times, results_matrix(:,7), '-pb');
    plot(times, results_matrix(:,8), 'or');
    plot(times, results_matrix(:,9), 'sk');
    grid on;
    ylim([0 80]);
    xlabel('Time [s]');
    ylabel('Oq');
    legend({'Oq from raw maximum','Oq from maximum after smoothing','Oq from peak detection','Oq from peak detection with smoothing'});

linkaxes(fig, 'x');
xlim([0, length(SIG)/fs]);



% Example for the DECOM method -------------------------------------------------

% Load sound and EGG files
[egg, fs] = audioread(['howtos' filesep 'DB-cresc-a-D4-m1.egg']);
[wav, fs] = audioread(['howtos' filesep 'DB-cresc-a-D4-m1.wav']);

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

figure
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
