% Recover the half spectrum corresponding to a Frequency Warped CEPstrum (FWCEP)
%
% Description
%  This function allow to recover the half spectrum which has been compressed
%  using Frequency Warped CEPstrum (FWCEP)
%
%  By default, the amplitude spectrum is unwarped using the frq2mel.m function
%  of the VOICEBOX toolbox. However, any other warping function can be used and
%  specified in the warpfn argument.
%  
%  See documentation of hspec2fwcep for more information about the compression.
%  
% Input
%  fwcep     : The frequency warped cepstral coefs to recover the spec from
%  fs        : [Hz] Sampling frequency.
%  dftlen    : Length of the DFT of the recovered spectrum.
% [warpfn]   : The frequency warping function.
%              Any functions of the "Frequency Scale Conversion" section of the
%              VOICEBOX can be used.
%              (def. frq2mel)
% [varargin] : Any additionnal argument for the warpfn function.
%
% Output
%  C         : The decompressed spectrum.
%
% See also
%  hspec2fwcep
%
% Reference
%  About the warping functions, please have a look at the references given in the
%  documentation of the corresponding functions in VOICEBOX (e.g. frq2mel.m)
%
% Copyright (c) 2011 University of Crete - Computer Science Department
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

function C = fwcep2hspec(fwcep, fs, dftlen, warpfn, varargin)

    if nargin<4 || isempty(warpfn); warpfn=@frq2mel; end

    % Compute the warping function
    freqlin = (0:dftlen/2)'*fs/dftlen;
    freqmel = 0.5*fs*frq2mel(freqlin)/frq2mel(fs/2);

    % Decode the spectrum from the cepstral coefficients
    Cwrapl = exp(fft(fwcep, dftlen));

    % Unwarp the spectrum
    C = interp1q(freqlin, abs(Cwrapl(1:end/2+1,:)), freqmel);

    if 0
        V3spec(Cwrapl, fs, 'g');
        V3spec(C, fs, 'b');
        keyboard
    end    

return