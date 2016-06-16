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
[wav, fs] = audioread([fname '.wav']);
times = (0:length(wav)-1)'/fs;

% Load an f0 curve
f0s = load([fname '.f0.txt']);
idx = f0s(:,2)>0.0;
f0s(:,2) = interp1(f0s(idx,1), f0s(idx,2), f0s(:,1), 'linear', 'extrap');

% Alternatively extract f0 with the SRH algorithm
[srh_f0,srh_vuv,srh_vuvc,srh_time] = pitch_srh(wav,fs,70,500,10);

dftlen = 4096;
hopdur = median(diff(f0s(:,1))); % Should be 5ms
ar_order = round(fs/1000+2);

disp('Estimate spectral peaks and spectra at regular intervals');
opt = sin_analysis();
opt.fharmonic  = false; % Use sinusoidal model (no constrains on the frequencies)
opt.use_ls     = false; % Use Peak Picking
opt.dftlen     = dftlen;  % Force the DFT length
opt.frames_keepspec = true; % Keep the computed spectra in the frames structure
frames = sin_analysis(wav, fs, f0s, opt);

srh_f0 = interp1(srh_time, srh_f0, [frames.t]);
srh_f0 = interp1_fixnan(srh_f0);
srh_vuv = interp1(srh_time, srh_vuv, [frames.t], 'nearest');
srh_vuv = interp1_fixnan(srh_vuv);
srh_f0s = [[frames.t]', srh_f0(:)];

disp('Compute Discrete All-Pole (DAP) envelope');
Edap = zeros(numel(frames),opt.dftlen/2+1);
pb = progressbar(numel(frames));
for n=1:numel(frames)

    [g a Edap(n,:)] = env_dap(frames(n).sins, fs, ar_order, opt.dftlen); % TODO DC has been dropped ???

    pb = progressbar(pb, n);
end

disp('Compute the so-called "True-Envelope"');
Ete = zeros(numel(frames),opt.dftlen/2+1);
Etec = zeros(numel(frames),opt.dftlen/2+1);
pb = progressbar(numel(frames));
for n=1:numel(frames)

    order = round(0.5*fs/frames(n).f0); % optimal cepstral order

    Ete(n,:) = env_te(hspec2spec(frames(n).S), order);

    fwcep = hspec2fwcep(Ete(n,:), fs, 24);
    Etec(n,:) = fwcep2hspec(fwcep, fs, opt.dftlen);

    pb = progressbar(pb, n);
end

disp('Compute WLP, SWLP, and XLP (LP variants), FIHR and FIHRZP');
order = 2*round(ceil(fs/1000+4)/2);
N = round(0.025*fs);
Nshift = round(hopdur*fs);
Nframes = ceil(length(wav)/Nshift - N/Nshift - 1);
win = hanning(N);
win = win./sum(win);
E_wlp = zeros(Nframes,dftlen/2+1);
E_swlp = zeros(Nframes,dftlen/2+1);
E_xlp = zeros(Nframes,dftlen/2+1);
E_fihr = zeros(Nframes,dftlen/2+1);
E_fihrzp = zeros(Nframes,dftlen/2+1);
time_lp = zeros(Nframes,1);
i = 1;
Flp = fs*(0:dftlen/2)/dftlen;
pb = progressbar(Nframes);
while i*Nshift+N < length(wav)
    
    % Get frame
    start = (i-1)*Nshift+1;
    stop = start+N-1;
    seg = wav(start:stop);

    f0 = interp1td(srh_f0s, (0.5*(start+stop)-1)/fs);

    % Save time
    time_lp(i) = (i*Nshift+N/2)/fs;

    % Estimate AR coefficients
    a = env_wlp_ste(seg,order);
    E_wlp(i,:) = mag2db(abs(gba2hspec(1,1,a,dftlen)));
    a = env_swlp_ste(seg,order);
    E_swlp(i,:) = mag2db(abs(gba2hspec(1,1,a,dftlen)));
    a = env_xlp_avs(seg,order);
    E_xlp(i,:) = mag2db(abs(gba2hspec(1,1,a,dftlen)));

    [a, e] = env_fihr(win.*seg, fs, f0, order);
    E_fihr(i,:) = mag2db(abs(gba2hspec(e,1,a,dftlen)));
    [a, e] = env_fihrzp(win.*seg, fs, f0, order);
    E_fihrzp(i,:) = mag2db(abs(gba2hspec(e,1,a,dftlen)));

    % Progress
    i = i + 1;
    pb = progressbar(pb, i);
end

disp('Discrete-Cepstral-Envelope Multi-Frame-Analysis (DCE-MFA)');
E_dcemfa = cell(1,numel(frames));
mfa_dur = 0.050; % [s]
mfa_size = round(0.5*mfa_dur/median(diff([frames.t])))*2+1;
ms = -(mfa_size-1)/2:(mfa_size-1)/2;
ns = 1-ms(1):numel(frames)-ms(end);
pb = progressbar(length(ns));
for n=ns
    order = round(0.5*fs/mean([frames(n+ms).f0]));

    [dum dum dum E] = env_dce_mfa(frames(n+ms), fs, order, true, [], 3000*(fs/16000), 0, dftlen);

    E_dcemfa{n} = E.';
    pb = progressbar(pb, n);
end
E_dcemfa = mag2db(abs(cell2mat(E_dcemfa')));


% Phase
disp('Compute the Relative Phase Shift and Phase Distortion');
opt = sin_analysis();
opt.fharmonic  = true; % Use harmonic model
opt.use_ls     = false; % Use Peak Picking
frames = sin_analysis(wav, fs, f0s, opt);
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
fig = [];
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
    title('Compressed/Decompressed TE envelope through mel freq. scale');

fig(5) = subplot(4,3,5);
    imagesc(time_lp, Flp, E_wlp.', [-10 50]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('WLP envelope');

fig(6) = subplot(4,3,8);
    imagesc(time_lp, Flp, E_swlp.', [-10 50]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('SWLP envelope');

fig(7) = subplot(4,3,11);
    imagesc(time_lp, Flp, E_xlp.', [-10 50]);
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('XLP envelope'); 
    xlabel('Time [s]');
    
fig(8) = subplot(4,3,6);
    imagesc([frames.t], Flp, E_fihr', [-120 -40])
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('FIHR envelope');

fig(9) = subplot(4,3,9);
    imagesc([frames.t], Flp, E_fihrzp', [-120 -40])
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('FIHRZP envelope');

fig(10) = subplot(4,3,12);
    imagesc([frames.t], Flp, E_dcemfa', [-100 -20])
    colormap(jet); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('DCE-MFA envelope');

linkaxes(fig, 'x');
xlim([0 times(end)]);


figure
fig = [];
ts = [frames.t];
F = fs*(0:opt.dftlen/2)/opt.dftlen;
fig(1) = subplot(2,1,1);
    imagesc(ts, F, RPS', [-pi pi]);
    colormap(circmap); freezeColors;
    axis xy;
    ylabel('Frequency [Hz]');
    title('Relative Phase Shift');
fig(2) = subplot(2,1,2);
    imagesc(ts, F, PD', [-pi pi]);
    colormap(circmap); freezeColors;
    axis xy;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('Phase Distortion');
