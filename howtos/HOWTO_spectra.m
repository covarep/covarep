% Some example for running spectral analysis
%
% Please find the path management in the startup.m script in the root directory
% of this repository. Note that by starting matlab in the root directory, this
% script should automatically run. If it is not the case, you can also move to the
% root directory and run this script manually. 
%
% Copyright (c) 2013 University of Crete - Computer Science Department
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
%  Gilles Degottex <degottex@csd.uoc.gr>
%

clear all;

% Load the waveform
fname = '0011.arctic_bdl1';
[wav, fs] = audioread([fname '.wav']);
times = (0:length(wav)-1)'/fs;

% Load an f0 curve
f0s = load([fname '.f0.txt']);

dftlen = 4096;
winlen = round(fs*0.05/2)*2+1;
win = blackman(winlen); win=win./sum(win);
%  F = fs*(0:dftlen/2)/dftlen;
hopsize = round(fs*5/1000);

disp('Short Time Fourier Transform (DFT based)')
[S, freqs, ts] = spectrogram(wav, win, winlen-hopsize,dftlen,fs);
ts = ts - 0.5/fs;

disp('Short Time Fan-Chirp Transform (FChT based)')
[FC, freqs, ts, as] = fchtgram(wav, win, winlen-hopsize, dftlen, fs, f0s);

% Plot the waveforms and the envelopes
figure
fig(1) = subplot(311);
    plot(times, wav, 'k');
    ylim(0.6*[-1 1]);
    xlabel('Time [s]');
    ylabel('Amplitude');
    title('Waveform');
fig(2) = subplot(312);
    imagesc(ts, freqs, mag2db(abs(S)), [-100 -20]);
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Short Time Fourier Transform (DFT based)');
fig(3) = subplot(313);
    imagesc(ts, freqs, mag2db(abs(FC)), [-100 -20]);
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Short Time Fan-Chirp Transform (FChT based)');
linkaxes(fig,'x');
xlim([0 times(end)]);


