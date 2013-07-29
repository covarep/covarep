% Some example for running the sinusoidal and harmonic estimators
%
% Copyright (c) 2012 University of Crete - Computer Science Department
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

% Load a rough f0 curve
% Can be replaced by any f0 estimate which is robust to octave jumps
f0s = load([fname '.f0.txt']);

disp(['Compute: Sinusoidal Model (SM) using Peak Picking']);
opt = sin_analysis();
opt.fharmonic  = false;
opt.use_ls     = false;
opt.resyn      = true; % Use the internal OLA method for the resynthesis
[frames syn_sm] = sin_analysis(wav, fs, f0s, opt);
wavwrite(syn_sm, fs, [fname '.sm.wav']);

disp(['Compute: Harmonic Model (HM)']);
opt = sin_analysis();
opt.fharmonic  = true;
opt.fadapted   = false;
opt.use_ls     = true;
frames = sin_analysis(wav, fs, f0s, opt);
syn_hm = hm_synthesis4(frames, length(wav), fs); % Use the harmonic resynthesis
wavwrite(syn_hm, fs, [fname '.hm.wav']);

disp(['Compute: Adaptive Harmonic Model (aHM) using the Adaptive Iterative Refinment (AIR)']);
optair = ahm_air_analysis2();
optair.do_final_ahm_step = false;
[f0sair, frames] = ahm_air_analysis2(wav, fs, f0s, optair);

disp('Compute the last step with uniform analysis instants');
f0sair = interp1td(f0sair, f0s(:,1));
opt = sin_analysis();
opt.fharmonic  = true;
opt.fadapted   = true;
opt.use_ls     = true;
frames = sin_analysis(wav, fs, f0sair, opt);
syn_ahm = hm_synthesis4(frames, length(wav), fs); % Use the harmonic resynthesis
wavwrite(syn_ahm, fs, [fname '.ahm-air.wav']);
%  save([fname '.ahm_frames.mat'], 'frames');

figure
plot(times, wav, 'k');
hold on;
plot(f0s(:,1), log2(f0s(:,2)), '--k');
plot(times, syn_sm, 'g');
plot(times, syn_hm, 'b');
plot(f0sair(:,1), log2 (f0sair(:,2)), 'r');
plot(times, syn_ahm, 'r');
xlabel('Time [s]');

legend({'Waveform', 'Input f_0', 'Sinusoidal Model (SM)', 'Harmonic Model', 'AIR refined f_0', 'Adaptive Harmonic Model'});
