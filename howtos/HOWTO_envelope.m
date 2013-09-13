% Some example for running the envelope estimators
% 
% Please find the path management in the startup.m script in the root directory
% of this repository. Note that by starting matlab in the root directory, this
% script should automatically run. If it is not the case, you can also move to the
% root directory and run this script manually. 
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
opt.fharmonic  = false; % Use sinusoidal model (no constrains on the frequencies)
opt.use_ls     = false; % Use Peak Picking
opt.dftlen     = 4096;  % Force the DFT length
opt.frames_keepspec = true; % Keep the computed spectra in the frames structure
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
Etec = zeros(numel(frames),opt.dftlen/2+1);
pb = progressbar(numel(frames));
for n=1:numel(frames)

    order = round(0.5*fs/frames(n).f0); % optimal cepstral order

    Ete(n,:) = env_te(hspec2spec(frames(n).S), order);

    mfcc = spec2mfcc(hspec2spec(Ete(n,:)), fs, 24);
    Etec(n,:) = mfcc2hspec(mfcc, fs, opt.dftlen);

    pb = progressbar(pb, n);
end

opt = phase_rpspd();
opt.harm2freq = true;
opt.pd_method = 1;
PD = phase_rpspd(frames, fs, opt);
opt.pd_method = 2;
opt.pd_vtf_rm = false;
opt.polarity_inv = true;
RPS = phase_rpspd(frames, fs, opt);

% Plot the waveforms and the envelopes
figure

F = fs*(0:opt.dftlen/2)/opt.dftlen;
fig(1) = subplot(421);
    plot(times, wav, 'k');
    grid on;
    ylim(0.6*[-1 1]);
    ylabel('Amplitude');
    title('Waveform');
fig(2) = subplot(423);
    imagesc([frames.t], F, mag2db(abs(Edap)).', [-100 -20]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('Discrete All-Pole (DAP) envelope');
fig(3) = subplot(425);
    imagesc([frames.t], F, mag2db(abs(Ete)).', [-100 -20]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('"True-Envelope" (TE)');
fig(4) = subplot(427);
    imagesc([frames.t], F, mag2db(abs(Etec)).', [-100 -20]);
    colormap(jet); freezeColors;
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Compressed/Decompressed TE envelope through MFCC');

fig(5) = subplot(422);
    plot(times, wav, 'k');
    grid on;
    ylim(0.6*[-1 1]);
    ylabel('Amplitude');
    title('Waveform');
fig(6) = subplot(424);
    F = fs*(0:opt.dftlen/2)/opt.dftlen;
    imagesc([frames.t], F, RPS', [-pi pi]);
    colormap(circmap); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('Relative Phase Shift');
fig(7) = subplot(426);
    F = fs*(0:opt.dftlen/2)/opt.dftlen;
    imagesc([frames.t], F, PD', [-pi pi]);
    colormap(circmap); freezeColors;
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Phase Distortion');

linkaxes(fig, 'x');
xlim([0 times(end)]);
