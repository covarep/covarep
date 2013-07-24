% Move the roots of a polynom inside the unit circle
%  
% Stabilize a polynom.
% Compared to polystab, polystab2 adds 2 elements
%  * Use the LEJA's order to recover the stable polynom from its roots in order
%    to minimize the numerical errors
%  * Even though the result should be stable in a single iteration, numerical
%    errors can still prevent the stability. Therefore, the stabilization step
%    is runned until the final result is finally stable.
%  
% $Id: polystab2.m 201 2012-12-14 09:40:31Z degottex $

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
