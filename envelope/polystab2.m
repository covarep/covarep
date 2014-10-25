% Move the roots of a polynom inside the unit circle
% 
% Description
%  Stabilize a polynom.
%  Compared to polystab, polystab2 adds 2 elements
%   * Use the LEJA's order to recover the stable polynom from its roots in order
%     to minimize the numerical errors
%   * Even though the result should be stable in a single iteration, numerical
%     errors can still prevent the stability. Therefore, the stabilization step
%     is runned until the final result is finally stable.
%
% Inputs
%  a        : The AR coefficients to stabilize
%  [modmax] : Maximum modulus of the roots.
%             Any root with a bigger modulus than this value will be moved
%             towards the origin.
%             (def. 1)
%
% Outputs
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

function b = polystab2(a, modmax)

    if nargin<2; modmax=1; end

    b = a;

    v = roots(a);

    idx = find(abs(v)>modmax);
    while ~isempty(idx)

        v(idx) = modmax*exp(1j*angle(v(idx)));

        ind = find(a~=0);
        b = a(ind(1))*poly(leja(v)); % Use LEJA's order to minimize num errors

        % Return only real coefficients if input was real:
        if ~any(imag(a))
            b = real(b);
        end

        a = b;
        v = roots(a);

        idx = find(abs(v)>modmax);
    end

    b = b(:);
return
