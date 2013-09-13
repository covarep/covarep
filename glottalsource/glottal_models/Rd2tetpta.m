% Gives estimate of T0-normalized LF parameters given an Rd value
%
% Description
%  The Rd shape parameter [1] of the Liljencrants-Fant (LF) glottal model [2]
%  expresses a regression of the original shape parameters {te,tp,ta}. This
%  function gives the shape parameters {te,tp,ta} which correspond to a given Rd
%  value.
%
% Input
%  Rd  : The Rd shape parameter to convert.
%
% Output
%  te  : Glottal shape parameter, see [1]p.6 (assuming T0=1 ! (T0-normalized))
%  tp  : Glottal shape parameter, see [1]p.6 (assuming T0=1 ! (T0-normalized))
%  ta  : Glottal shape parameter, see [1]p.6 (assuming T0=1 ! (T0-normalized))
%
% See Also
%  gfm_spec_lf.m
%  
% References
%  [1] G. Fant, "The LF-model revisited. Transformations and frequency domain
%      analysis", STL-QPSR 36(2-3):119-156, 1995.
%  [2] G. Fant, J. Liljencrants and Q. Lin, "A four-parameter model of glottal
%      flow", STL-QPSR, vol. 4, pp. 1-13, 1985.
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

function [te tp ta] = Rd2tetpta(Rd)

	Rap = (-1+4.8.*Rd)/100;                             % [1](2)
	Rkp = (22.4+11.8.*Rd)/100;                          % [1](3)
	Rgp = 1./(4*((0.11.*Rd./(1/2+1.2.*Rkp))-Rap)./Rkp); % [1] indirectly (4)

	tp = 1./(2.*Rgp); % [1]p.121
	te = tp.*(Rkp+1); % [1]p.121
	ta = Rap;         % [1]p.121

return
