% Minimum-phase log spectrum of a given spectrum using the complex cepstrum
%
% Description
%  Compute the minimum-phase spectrum of a given spectrum X by dropping the
%  negative quefrencies of X.
%
% Input
%  X   : The spectrum to convert (full DFT length)
%
% Output
%  lM  : The minimum-phase log spectrum (half DFT length)
%
% References
%  [1] Alan V. Oppenheim and Ronald W. Schafer, "Digital Signal Processing",
%      Prentice-Hall, 2nd edition, 1978.
%      Note that this 2nd edition contains a chapter about complex cepstrum which
%      has been removed in the 3rd edition (and possibly came back in the 4th ed.)
%
% Copyright (c) 2012 University of Crete - Computer Science Department (UOC-CSD)
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

function lM = hspec2minphaseloghspec(X)

    X = X(:).';

    rcc = ifft(hspec2spec(log(abs(X)))); % Compute the real cepstrum

    if mod(length(rcc),2)==0
        % For even DFT length
        lM = fft([rcc(1),2*rcc(2:end/2),rcc(end/2+1)], length(rcc));% [1]
    else
        % For odd DFT length
        lM = fft([rcc(1),2*rcc(2:(end-1)/2+1)], length(rcc));       % [1]
    end

    lM = lM(1:end/2+1);

return
