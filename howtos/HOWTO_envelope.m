% Some example for running the envelope estimators
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
[wav, fs] = wavread([fname '.wav']);
times = (0:length(wav)-1)'/fs;

% Load an f0 curve
f0s = load([fname '.f0.txt']);

disp('Estimate spectral peaks and spectra at regular intervals');
opt = sin_analysis();
opt.fharmonic  = false;
opt.use_ls     = false;
opt.dftlen     = 4096;  % Force the DFT length
opt.frames_keepspec = true;
frames = sin_analysis(wav, fs, f0s, opt);

disp('Compute Discrete All-Pole (DAP) envelope');
order = round(fs/1000+2);
Edap = zeros(numel(frames),opt.dftlen/2+1);
pb = progressbar(numel(frames));
for n=1:numel(frames)

    [g a Edap(n,:)] = env_dap(frames(n).sins, fs, order, opt.dftlen);

    pb = progressbar(pb, n);
end

disp('Compute the so-called "True-Envelope"');
Ete = zeros(numel(frames),opt.dftlen/2+1);
pb = progressbar(numel(frames));
for n=1:numel(frames)

    order = round(0.5*fs/frames(n).f0); % optimal cepstral order

    Ete(n,:) = env_te(hspec2spec(frames(n).S), order);

    pb = progressbar(pb, n);
end

% Plot the waveforms and the envelopes
figure
F = fs*(0:opt.dftlen/2)/opt.dftlen;
fig(1) = subplot(311);
    plot(times, wav, 'k');
    ylim(0.6*[-1 1]);
    xlabel('Time [s]');
    ylabel('Amplitude');
    title('Waveform');
fig(2) = subplot(312);
    imagesc([frames.t], F, mag2db(abs(Edap)).', [-100 -20]);
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Discrete All-Pole (DAP) envelope');
fig(3) = subplot(313);
    imagesc([frames.t], F, mag2db(abs(Ete)).', [-100 -20]);
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('"True-Envelope" (TE)');
linkaxes(fig, 'x');
xlim([0 times(end)]);
