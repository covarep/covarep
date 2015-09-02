% Function to automatically detect creaky voice 
%
% Description
%  Function to automatically detect creaky voice 
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
%  FeatTot  : [samples] [Nx36] Feature matrix
%  ANN      : [struct]  [1x1]  Trained ANN
%
% Outputs
%  creak_pp : [samples] [Nx1] Estimated posterior probability of creaky
%                             voice
%  creak_bin: [samples] [Nx1] Binary creaky voice decision
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

function [creak_pp,creak_bin] = creaky_voice_do_detection(FeatTot,ANN)

%% Initial settings
Maxis=ANN.Maxis.Maxis;
Minis=ANN.Minis.Minis;
net=ANN.net.net;

ANN_decision_threshold=0.3; % Optimised in publication
med_len=3; % Length of median filtering

%% Do Z-score normalisation
m=size(FeatTot);
X=zeros(m(2),m(1));
for k=1:m(2)
    X(k,:)=FeatTot(:,k)';    
end

m=size(X);

for k=1:m(1)
    vec=X(k,:);
    
    mini=Minis(k);
    maxi=Maxis(k);
    
    pos=isnan(vec);
    vec(pos)=mini;
    
    X(k,:)=-1+((vec-mini)/(maxi-mini))*2;
end

%% Simulate ANN network
creak_pp=sim(net,X);

creak_pp2=medfilt1(creak_pp,med_len); % Median filtering to remove transients

creak_bin=zeros(1,length(creak_pp));
creak_bin(creak_pp2>ANN_decision_threshold)=1;
