% Function to automatically detect creaky voice (otherwise known as
% vocal/glottal fry)
%
% Description
%  Function to detect creaky voice using novel acoustic features developed
%  by Kane & Drugman (KD features) as well as previous acoustic features
%  developed by Carlos Ishi and colleagues. Detection is done using artificial
%  neural networks
%
%  This version is a slightly later one than that described in the
%  above published 2013 CSL paper [2]. The algorithm here rather than using binary
%  decision trees using artificial neural networks and combines the
%  features used in the CSL paper with those proposed in Ishi et al.
%  (2008). This updated version has been submitted to CSL for a special
%  issue on glottal source processing on April 14th 2013. It will have
%  reference [1].
%
% Inputs
%  x        : [samples] [Nx1] Speech signal
%  fs       : [Hz]      [1x1] sampling frequency
%
% Outputs
%  creak_pp : [samples] [Nx2] Estimated posterior probability of creaky
%                             voice, with time (samples)
%  creak_bin: [samples] [Nx2] Binary creaky voice decision, with time
%                             (samples)
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%  [2] Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', Computer Speech and Language 27(4), pp.
%       1028-1047.
%  [3] Drugman, T., Kane, J., Gobl, C., (2012) `Resonator-based creaky 
%       voice detection', Interspeech 2012, Portland, Oregon, USA.
%
% Octave Compatibility
%  As this algorithm relies on the Matlab neural networks toolbox it is not
%  compatible with Octave. Furthermore, this algorithm requires Matlab
%  R2011 or later
%
% Copyright (c) 2013 University of Mons, FNRS & 2013 Trinity College Dublin
%
% License
%  This code is a part of the GLOAT toolbox with the following
%  licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Authors
%  Thomas Drugman <thomas.drugman@umons.ac.be> & John Kane <kanejo@tcd.ie>
%

function [creak_pp,creak_bin] = detect_creaky_voice(x,fs)


%% Load files
ANN.net=load('system_net_creak');
ANN.Maxis=load('maxis_creak.mat');
ANN.Minis=load('minis_creak.mat');

%% Do feature extraction
FeatMat = get_all_creak_features(x,fs);

%% Do classification
[creak_pp,creak_bin] = creaky_voice_do_detection(FeatMat,ANN);
t=(1:length(creak_pp))*10/1000;
creak_pp=[creak_pp(:) round(t(:)*fs)];
creak_bin=[creak_bin(:) round(t(:)*fs)];

