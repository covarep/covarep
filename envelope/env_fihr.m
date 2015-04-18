% Fast Inter-Harmonic Reconstruction as a pre-process for LPC-based
% spectral envelope estimation.
%
%
% Description
%  This technique reconstructs inter-harmonics of the voice signal,
% as a pre-process for an efficient spectral envelope estimation in
% high-pitched voices. The technique is fully described in [1].

%
% Inputs
%  wave             : [samples] [Nx1] input signal (speech signal)
%  Fs               : [Hz]      [1x1] sampling frequency
%  f0               : F0 estimate (a 0 value is for an unvoiced frame)
%  order            : Order used for the LP analysis
%
%
% Outputs
%  LP               : LP coefficients
%  Energy           : energy of the frame
% 
%
% Example
%  Please see the HOWTO_envelope.m example file.
%  Please see http://tcts.fpms.ac.be/~drugman/Toolbox/ for more details.
%
% References
%  [1] T.Drugman, Y. Stylianou, "Fast Inter-Harmonic Reconstruction for
%      Spectral Envelope Estimation in High-Pitched Voices", vol. 21, pp.
%      1418-1422, IEEE Signal Processing Letters, 2014.
%      http://tcts.fpms.ac.be/~drugman/files/SPL-FIHR.pdf
%
% Copyright (c) 2014 Toshiba Cambridge Research Laboratory
%
% License
%  This code will be part of the GLOAT toolbox with the following
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
% This function will also be part of the Covarep project: http://covarep.github.io/covarep
%
% Author
%  Thomas Drugman thomas.drugman@umons.ac.be

function [LP,Energy] = env_fihr(Seg,Fs,F0_local,order)


F0_target=100;
% F0* in the paper. Here set to 100 Hz.

Energy=sum(Seg.^2);

if F0_local>0
    
    TmpSignal=ones(1,length(Seg));
    t_tmp=1:length(TmpSignal);
    
    I=ceil(F0_local/F0_target);
    for k=1:I-1
        W=(I-k)/I;
        % We use here a linear weighting function. Other functions are
        % possible, as long as they meet the properties mentioned in
        % the paper.
        TmpSignal=TmpSignal+2*W*cos(2*pi*(k/I)*F0_local/Fs*t_tmp);
        
        % Equation (2) in the paper
    end
    
    Seg2=Seg.*TmpSignal';
else
    Seg2=Seg;
end

LP=lpc(Seg2,order);
