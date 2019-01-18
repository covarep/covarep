% Voice activity detection.
%
% Description
%  This is the voice activity detection (VAD) algorithm which is described
%  in [1]. The VAD exploits 3 sets of features: MFCCs as filter-based
%  features, and 2 sets of source-related features: 4 features presented in
%  Sadjadi's paper [2] and 3 proposed features. These features are also
%  used for robust ASR in [3].
%
% Inputs
%  wave            : [samples] [Nx1] input signal (speech signal)
%  Fs              : [Hz]      [1x1] sampling frequency
%  doPlot          : flag to show the results (=0 no plot, else plot the
%                    results). Default value=0.
%   
% Outputs
%  [Outs_Final,Outs_MFCC,Outs_Sadjadi,Outs_New] : these are 4 vectors containing
%                   the VAD posteriors using respectively: i) the combined
%                   system which makes use of a decision fusion strategy
%                   and is based on the 3 feature sets, ii) the system
%                   using only MFCCs, iii) the system using only Sadjadi's
%                   features, iv) the system using only the proposed
%                   features.
%  t            : [seconds] Instants of the VAD posteriors.
%
% Example
%  Please see the HOWTO_VAD.m example file.
%
% References
%  [1] T.Drugman, Y. Stylianou, Y. Kida, M. Akamine: "Voice Activity Detection:
%  Merging Source and Filter-based Information", IEEE Signal Processing Letters,
%  2015.
%  http://tcts.fpms.ac.be/~drugman/files/SPL_VAD.pdf
%  [2] S.O. Sadjadi, J. Hansen: "Unsupervised Speech Activity Detection Using
%  Voicing Measures and Perceptual Spectral Flux", IEEE Sig. Pro. Letters,
%  vol. 20, pp. 197-200, 2013.
%  [3] T. Drugman, Y. Stylianou, L. Chen, X. Chen, M. Gales, Robust
%  Excitation-based Features for Automatic Speech Recognition, IEEE 
%  International Conference on Acoustics, Speech and Signal Processing (ICASSP),
%   2015: http://tcts.fpms.ac.be/~drugman/files/ICASSP15_ASR.pdf
%
% Copyright (c) 2014 Toshiba Cambridge Research Laboratory
%
% License
%  This code will be part of the GLOAT toolbox
%   (http://tcts.fpms.ac.be/~drugman/Toolbox/)
%  with the following licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is also be part of the Covarep project:
%   http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be

function [Outs_Final,Outs_MFCC,Outs_Sadjadi,Outs_New,t] = VAD_Drugman(wave,Fs,doPlot)

if nargin<3
    doPlot=0;
end

if max(abs(wave))>1
    wave=wave/max(abs(wave));
end

% Feature Extraction
[Feat,t] = VAD_Drugman_FeatureExtraction(wave,Fs);

% The feature trajectories are smoothed using a median-filter and the first
% and second derivatives are appended
Feat_MFCC=Feat(1:13,:);
[Feat_MFCC] = Feature_FilteringAndDerivatives(Feat_MFCC);
Feat_Sadjadi=Feat(14:17,:);
[Feat_Sadjadi] = Feature_FilteringAndDerivatives(Feat_Sadjadi);
Feat_New=Feat(18:20,:);
[Feat_New] = Feature_FilteringAndDerivatives(Feat_New);

% Classification using only the MFCCs
load(['Minis_MFCC.mat'])
load(['Maxis_MFCC.mat'])
load(['ANNSystem_MFCC.mat'])
X=Feat_MFCC';
for k=1:size(X,1)
    vec=X(k,:);
    mini=Minis(k);
    maxi=Maxis(k);
    X(k,:)=-1+((X(k,:)-mini)/(maxi-mini))*2;
end
Outs=sim(net,X);
Outs_MFCC=medfilt1(Outs,11);

% Classification using only Sadjadi's features
load(['Minis_Sadjadi.mat'])
load(['Maxis_Sadjadi.mat'])
load(['ANNSystem_Sadjadi.mat'])
X=Feat_Sadjadi';
for k=1:size(X,1)
    vec=X(k,:);
    mini=Minis(k);
    maxi=Maxis(k);
    X(k,:)=-1+((X(k,:)-mini)/(maxi-mini))*2;
end
Outs=sim(net,X);
Outs_Sadjadi=medfilt1(Outs,11);

% Classification using only the new features
load(['Minis_New.mat'])
load(['Maxis_New.mat'])
load(['ANNSystem_New.mat'])
X=Feat_New';
for k=1:size(X,1)
    vec=X(k,:);
    mini=Minis(k);
    maxi=Maxis(k);
    X(k,:)=-1+((X(k,:)-mini)/(maxi-mini))*2;
end
Outs=sim(net,X);
Outs_New=medfilt1(Outs,11);

% Final classification using a decision fusion (simple geometrical mean)
Outs_Final=(Outs_MFCC.*Outs_Sadjadi.*Outs_New).^(1/3);

if doPlot
    t2=(1:length(wave))/Fs;
    plot(t2,wave)
    hold on
    plot(t,Outs_MFCC,'k')
    plot(t,Outs_Sadjadi,'g')
    plot(t,Outs_New,'c')
    plot(t,Outs_Final,'r')
    xlabel('Time (s)')
    legend('Speech waveform','Posterior from the VAD using MFCCs','Posterior from the VAD using Sadjadi''s features','Posterior from the VAD using proposed fetaures','Posterior from the proposed VAD using the 3 features sets')
end


function [MatFeat] = Feature_FilteringAndDerivatives(Feat)

Nfeats=size(Feat,1);
Nframes=size(Feat,2);

%Median filtering
for k=1:Nfeats
    Feat(k,:)=medfilt1(Feat(k,:),11);
end

MatFeat=zeros(Nframes,3*Nfeats);
MatFeat(:,1:Nfeats)=Feat';

% Computing the derivatives
Naround=10;
for k=1:Nfeats
    for l=Naround+1:Nframes
        Vec=Feat(k,l-Naround:l);
        Val=0;
        for p=1:length(Vec)-1
            Val=Val+(Vec(end)-Vec(p));
        end
        MatFeat(l,k+Nfeats)=Val;
    end
end

for k=1:Nfeats
    for l=Naround+1:Nframes
        Vec=MatFeat(l-Naround:l,k+Nfeats);
        Val=0;
        for p=1:length(Vec)-1
            Val=Val+(Vec(end)-Vec(p));
        end
        MatFeat(l,k+2*Nfeats)=Val;
    end
end