% Compute Mel-Frequency Cepstral Coefficients (MFCC) from a spectrum
%
% Description
%  This function allow to compute the well-known MFCC of a given spectrum (e.g.
%  an envelope).
%  By default, the amplitude spectrum is warped using the frq2mel.m function of
%  the VOICEBOX toolbox. However, any other warping function can be used and
%  specified in the warpfn argument.
%  Note that the energy wich is lost by droping the negative cepstral coefficients
%  is compensated by multiplying the positive cepstral coefficients by 2, thus
%  preserving the energy between the spectrum and the MFCCs.
%
% Input
%  C         : The full spectrum to compress (of the size of the DFT)
%  fs        : [Hz] Sampling frequency
%  order     : Order of the cepstrum (without counting the energy coefficient)
%              The size of the mfcc vector is always 1+order !
% [warpfn]   : The frequency warping function.
%              Any functions of the "Frequency Scale Conversion" section of the
%              VOICEBOX can be used.
%              (def. frq2mel, in order to compute the MFCC)
% [varargin] : Any additionnal argument for the warpfn function.
%
% Output
%  mfcc      : Frequency warped coefficients
%
% See also
%  mfcc2hspec
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

function mfcc = spec2mfcc(C, fs, order, warpfn, varargin)

    if nargin<4 || isempty(warpfn); warpfn=@frq2mel; end

    C = C(:);

    % The input C is assumed to be a full spectrum
    dftlen = length(C);

    % Compute the warping function
    freqlin = (0:dftlen/2)'*fs/dftlen;
    freqmel = 0.5*fs*warpfn(freqlin, varargin{:})/warpfn(0.5*fs, varargin{:});

    % Warp the spectrum
    env = interp1q(freqmel, abs(C(1:end/2+1)), freqlin);
    if isnan(env(end)); env(end)=env(end-1); end
    
    % Symmetrize the warped spectrum prior to cepstral computation
    Cwrap = [abs(env(1)); env(2:end-1); abs(env(end)); conj(env(end-1:-1:2))];

    % Compute the cepstrum
    mfcc = ifft(log(abs(Cwrap)));
    
    % Drop the negative quefrencies and compensate the loss of energy
    mfcc = [mfcc(1); 2*mfcc(2:1+order)];
    
    if 0
        Cmel = mfcc2spec(mfcc, fs, dftlen, warpfn, varargin{:});

        V3spec(C, fs, 'k');
        V3spec(Cwrap, fs, 'r');
        V3spec(Cmel, fs, 'b');
        keyboard
    end
return
