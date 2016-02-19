% Compute Frequency Warped CEPstrum (FWCEP) from a half spectrum
%
% Description
%  This function allows to compute frequency warped cepstral coefficients
%  of a given spectrum (e.g. an envelope).
%  
%  By default, the amplitude spectrum is warped using the frq2mel.m function of
%  the VOICEBOX toolbox. However, any other warping function can be used and
%  specified in the warpfn argument.
%
%  Note that the energy wich is lost by droping the negative cepstral coefficients
%  is compensated by multiplying the positive cepstral coefficients by 2, thus
%  preserving the energy between the spectrum and the cepstral coefficients.
%
%  By using the mel frequency scale, this implementation is very similar to the
%  initial step of the Mel CEPstral computation (MCEP) [2]. The main difference
%  is that, here, the frequency scale is compressed in frequency domain using
%  simple linear interpolation, whereas in [2], it is compressed in cepstral
%  domain.
%
% Input
%  C         : The half spectrum to compress. If C is a vector, C should have
%              the size [dftlen/2+1 1]. C can be a matrix of size [dftlen/2+1 n]
%  fs        : [Hz] Sampling frequency
% [order]    : Order of the cepstrum
%              (Without counting the 0-coefficient. The size of the fwcep vector
%              is thus always 1+order)
% [warpfn]   : The frequency warping function.
%              Any functions of the "Frequency Scale Conversion" section of the
%              VOICEBOX can be used.
%              (def. frq2mel, in order to compute the MCEP-like coefficients)
% [varargin] : Any additionnal argument for the warpfn function.
%
% Output
%  fwcep      : Frequency warped coefficients
%
% See also
%  fwcep2hspec
%
% Reference
%  [1] K. Tokuda, T. Kobayashi, T Masuko and S. Imai, "Mel-generalized cepstral
%      analysis - A unified approach to speech spectral estimation", Proceedings
%      of International Conference on Spoken Language Processing, vol.3,
%      pp.1043-1046, 1994.
%  [2] T. Fukada, K. Tokuda, T. Kobayashi, and S. Imai, "An adaptive algorithm
%      for mel-cepstral analysis of speech," in IEEE International Conference on
%      Acoustics, Speech and Signal Processing (ICASSP), vol. 1, p. 137-140
%      vol.1, 1992.
%  * About the warping functions, please have a look at the references given in
%    the documentation of the corresponding functions in VOICEBOX (e.g. frq2mel.m)
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

function fwcep = hspec2fwcep(C, fs, order, warpfn, varargin)

    if size(C,1)==1, C = C'; end

    % The input C is assumed to be a half spectrum
    dftlen = (size(C,1)-1)*2;

    if nargin<3 || isempty(order); order=dftlen/2; end
    if nargin<4 || isempty(warpfn); warpfn=@frq2mel; end

    % Compute the warping function
    freqlin = (0:dftlen/2)'*fs/dftlen;
    freqmel = 0.5*fs*warpfn(freqlin, varargin{:})/warpfn(0.5*fs, varargin{:});

    % Warp the spectrum
    env = interp1q(freqmel, abs(C), freqlin);
    idx = isnan(env(end,:));
    if any(idx)
        env(end,idx)=env(end-1,idx);
    end
    
    % Symmetrize the warped spectrum prior to cepstral computation
    Cwrap = [env; conj(env(end-1:-1:2,:))];

    % Compute the cepstrum
    fwcep = ifft(log(abs(Cwrap)),'symmetric');

    % Drop the negative quefrencies and compensate the loss of cepstral energy
    fwcep = [fwcep(1,:); 2*fwcep(2:1+order,:)];
end
