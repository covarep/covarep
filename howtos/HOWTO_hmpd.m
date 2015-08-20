% Harmonic Model + Phase Distortion (HMPD) vocoder
%
% Copyright (c) 2014 University of Crete - Computer Science Department(UOC-CSD)/ 
%                    Foundation for Research and Technology-Hellas - Institute
%                    of Computer Science (FORTH-ICS)
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

fname = '0011.arctic_bdl1';
[wav, fs] = audioread([fname '.wav']);


% Analysis ---------------------------------------------------------------------

hmpdopt = hmpd_analysis();
hmpdopt.f0min   = 75; % To adapt to the analyzed voice
hmpdopt.f0max   = 220;% To adapt to the analyzed voice

% Speed up options
%  hmpdopt.sin.use_ls=false;   % For analysis step, use simple Peak Picking 
%  hmpdopt.sin.fadapted=false; % For analysis step, use the stationary DFT 
%  hmpdopt.usemex = true; % Use mex linear interp (affects mainly the synth.)

% Compression options
%  hmpdopt.amp_enc_method=2; hmpdopt.amp_log=true; hmpdopt.amp_order=39;
%  hmpdopt.pdd_log=true; hmpdopt.pdd_order=12;% MFCC-like phase variance
%  hmpdopt.pdm_log=true; hmpdopt.pdm_order=24;% Number of log-Harmonic coefs

[f0s, AE, PDM, PDD] = hmpd_analysis(wav, fs, [], hmpdopt);


% Synthesis --------------------------------------------------------------------

synopt = hmpd_synthesis();
synopt.enc = hmpdopt; 
synopt.usemex = hmpdopt.usemex; % Speed up with mex function
syn = hmpd_synthesis(f0s, AE, [], PDD, fs, length(wav), synopt);
audiowrite([fname '.hmpd-pdd.wav'], syn, fs);

if ~hmpdopt.pdm_log
    syn = hmpd_synthesis(f0s, AE, PDM, PDD, fs, length(wav), synopt);
    audiowrite([fname '.hmpd-pdmpdd.wav'], syn, fs);
end


% Plot -------------------------------------------------------------------------
figure
fig(1) = subplot(411);
    times = (0:length(wav)-1)'/fs;
    plot(times, wav, 'k');
    hold on;
    plot((0:length(syn)-1)/fs, syn, 'b');
    plot(f0s(:,1), log2(f0s(:,2))-log2(440), 'r');
    title('Waveform');
    xlabel('Time [s]');
fig(2) = subplot(412);
    F = fs*(0:hmpdopt.dftlen/2)/hmpdopt.dftlen;
    if ~hmpdopt.amp_log; imagesc(f0s(:,1), F, mag2db(AE)', [-120 -20]);
    else;                imagesc(f0s(:,1), (0:hmpdopt.amp_order), AE'); end
    colormap(jet);
    freezeColors;
    axis xy;
    title('Amplitude Envelope');
    xlabel('Time [s]');

fig(3) = subplot(413);
    F = fs*(0:hmpdopt.dftlen/2)/hmpdopt.dftlen;
    if ~hmpdopt.pdd_log; imagesc(f0s(:,1), F, PDM', [-pi pi]);
    else;                imagesc(f0s(:,1), (0:hmpdopt.pdm_order), PDM', [-pi pi]); end
    colormap(circmap);
    freezeColors;
    axis xy;
    title('Phase Distortion Mean (PDM)');
    xlabel('Time [s]');

fig(4) = subplot(414);
    if ~hmpdopt.pdd_log; imagesc(f0s(:,1), F, PDD', [0 2]);
    else;                imagesc(f0s(:,1), (0:hmpdopt.pdd_order), PDD'); end
    colormap(jet);
    freezeColors;
    axis xy;
    title('Phase Distortion Deviation (PDD)');
    xlabel('Time [s]');

linkaxes(fig, 'x');
xlim([times(1) times(end)]);
