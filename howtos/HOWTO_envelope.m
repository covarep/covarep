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

disp('Compute WLP, SWLP, and XLP (LP variants)');
fftlen = 4096;
order = 2*round(ceil(fs/1000+4)/2);
N = round(0.025*fs);
Nshift = round(0.005*fs);
Nframes = ceil(length(wav)/Nshift - N/Nshift - 1);
E_wlp = zeros(Nframes,fftlen);
E_swlp = zeros(Nframes,fftlen);
E_xlp = zeros(Nframes,fftlen);
time_lp = zeros(Nframes,1);
i = 1;
Flp = fs*(0:fftlen/2)/fftlen;
pb = progressbar(Nframes);
while i*Nshift+N < length(wav)
    
    % Get frame
    frame = wav((i-1)*Nshift+1:i*Nshift+N);
    
    % Save time
    time_lp(i) = (i*Nshift+N/2)/fs;
    
    % Estimate AR coefficients
    a1 = env_wlp_ste(frame,order);
    a2 = env_swlp_ste(frame,order);
    a3 = env_xlp_avs(frame,order);
    
    % Estimate envelopes
    E_wlp(i,:) = db(abs(freqz(1,a1,fftlen,fs)));
    E_swlp(i,:) = db(abs(freqz(1,a2,fftlen,fs)));
    E_xlp(i,:) = db(abs(freqz(1,a3,fftlen,fs)));

    % Progress
    i = i + 1;
    pb = progressbar(pb, i);
end

% Phase
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
fig(1) = subplot(4,3,1);
    plot(times, wav, 'k');
    grid on;
    ylim(0.6*[-1 1]);
    ylabel('Amplitude');
    title('Waveform');
fig(2) = subplot(4,3,4);
    imagesc([frames.t], F, mag2db(abs(Edap)).', [-100 -20]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('Discrete All-Pole (DAP) envelope');
fig(3) = subplot(4,3,7);
    imagesc([frames.t], F, mag2db(abs(Ete)).', [-100 -20]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('"True-Envelope" (TE)');
fig(4) = subplot(4,3,10);
    imagesc([frames.t], F, mag2db(abs(Etec)).', [-100 -20]);
    colormap(jet); freezeColors;
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Compressed/Decompressed TE envelope through MFCC');

    
fig(5) = subplot(4,3,2);
    plot(times, wav, 'k');
    grid on;
    ylim(0.6*[-1 1]);
    ylabel('Amplitude');
    title('Waveform');
fig(6) = subplot(4,3,5);
    F = fs*(0:opt.dftlen/2)/opt.dftlen;
    imagesc([frames.t], F, RPS', [-pi pi]);
    colormap(circmap); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('Relative Phase Shift');
fig(7) = subplot(4,3,8);
    F = fs*(0:opt.dftlen/2)/opt.dftlen;
    imagesc([frames.t], F, PD', [-pi pi]);
    colormap(circmap); freezeColors;
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Phase Distortion');
    
    
fig(8) = subplot(4,3,3);
    plot(times, wav, 'k');
    grid on;
    ylim(0.6*[-1 1]);
    ylabel('Amplitude');
    title('Waveform');
    
fig(9) = subplot(4,3,6);
    imagesc(time_lp, Flp, E_wlp.', [-25 65]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('WLP envelope');
    
fig(10) = subplot(4,3,9);
    imagesc(time_lp, Flp, E_swlp.', [-25 65]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('SWLP envelope');
    
fig(11) = subplot(4,3,12);
    imagesc(time_lp, Flp, E_xlp.', [-25 65]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('XLP envelope');  
    
    
    
linkaxes(fig, 'x');
xlim([0 times(end)]);
