% Minimum-phase spectrum of a given spectrum using the complex cepstrum
%
% Description
%  Compute the minimum-phase spectrum of a given spectrum X by dropping the
%  negative quefrencies of X.
%
% Input
%  X   : The half spectrum to convert (half DFT length)
%
% Output
%  M   : The corresponding minimum-phase spectrum (half DFT length)
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

function M = hspec2minphasehspec(X)

    M = exp(hspec2minphaseloghspec(X));

return
