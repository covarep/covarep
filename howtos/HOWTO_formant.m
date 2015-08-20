% Some example for running formant tracking algorithms
%
% Please find the path management in the startup.m script in the root directory
% of this repository. Note that by starting matlab in the root directory, this
% script should automatically run. If it is not the case, you can also move to the
% root directory and run this script manually. 
%
% Copyright (c) 2013 University of Mons, FNRS
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
%  Thomas Drugman thomas.drugman@umons.ac.be
%

clear all;

% Load soundfile
[x,fs] = audioread('arctic_a0007.wav');

[formantPeaks,t_analysis]=formant_CGDZP(x,fs,30,10);

t=(0:length(x)-1)/fs;
fig(1) = subplot(211);
    plot(t,x, 'b');
    legend('Speech signal');
    xlabel('Time [s]');
    ylabel('Amplitude');

    fig(2) = subplot(212);
    plot(t_analysis,formantPeaks(:,1), 'b');
    hold on
    plot(t_analysis,formantPeaks(:,2), 'r');
    plot(t_analysis,formantPeaks(:,3), 'g');
    plot(t_analysis,formantPeaks(:,4), 'm');
    plot(t_analysis,formantPeaks(:,5), 'c');
    legend('First formant','Second formant','Third formant','Fourth formant','Fifth formant');
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');


linkaxes(fig, 'x');
xlim([0 t(end)]);
