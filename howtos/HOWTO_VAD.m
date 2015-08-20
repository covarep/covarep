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
%  Thomas Drugman <thomas.drugman@umons.ac.be>

% Read any speech audio file
[x,fs] = audioread(['howtos' filesep 'arctic_a0007.wav']);

% VAD algorithm
% The third input is just a flag to plot the results or not
% The outputs are the following:
%  [Outs_Final,Outs_MFCC,Outs_Sadjadi,Outs_New]  : these are 4 vectors containing
%                   the VAD posteriors using respectively: i) the combined
%                   system which makes use of a decision fusion strategy
%                   and is based on the 3 feature sets, ii) the system
%                   using only MFCCs, iii) the system using only Sadjadi's
%                   features, iv) the system using only the proposed
%                   features.
%  t            : [seconds] Instants of the VAD posteriors.

try
    [Outs_Final,Outs_MFCC,Outs_Sadjadi,Outs_New,t] = VAD_Drugman(x,fs,1);
catch
    disp('Version or toolboxes do not support neural network object used in VAD. VAD skipped.')
end
