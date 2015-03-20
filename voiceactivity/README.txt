This a Matlab code for Voice Activity Detection (VAD).

The "HOWTO_VAD.m" script gives an example on how to use the package.
It builds upon "VAD_Drugman.m" which provides the posterior probability of voice 
activity, using any of 3 feature sets independently, or combining them in a 
decision fusion strategy.
These 3 feature sets are: MFCCs, 4 voicing measurements proposed in Sadjadi's 
paper in SPL 2013 (see reference below), and 3 new proposed source-related 
features.

According to our experiments, the combined decisions (called "Outs_Final" in the 
code), provides the best results.

All details are provided in the paper:
T.Drugman, Y. Stylianou, "Voiced Activity Detection: Merging Source and Filter-
based Information", IEEE Signal Processing Letters, 2015.

These features are also the basis of:
T. Drugman, Y. Stylianou, L. Chen, X. Chen, M. Gales, Robust Excitation-based 
Features for Automatic Speech Recognition, IEEE International Conference on 
Acoustics, Speech and Signal Processing (ICASSP), 2015 

Please refer to these papers in your publication.

Sadjadi's features are described in:
S.O. Sadjadi, J. Hansen: "Unsupervised Speech Activity Detection Using Voicing 
Measures and Perceptual Spectral Flux", IEEE Sig. Pro. Letters, vol. 20, pp. 
197-200, 2013.


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
% This function will also be part of the Covarep project:
%  http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be