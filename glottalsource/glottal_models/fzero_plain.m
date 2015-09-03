%% This function contains copies of algorithms from MATLAB's fzero function.
% MATLAB's copyright may apply to this files: Copyright 1984-2013 The MathWorks, Inc.
%
% fzero_plain is a simple version of MATLAB's fzero function. Since it does
% not need to process several options and set default values, it is faster
% than MATLAB's fzero (for simple problems). fzero_plain only contains the
% plain algorithm's of fzero. 
function b = fzero_plain(fun,x)

tol = eps;

if numel(x)==1
    %% find interval
    fx = fun(x);
    
    if fx == 0
        b = x;
        return;
    end
    
    if x ~= 0,
        dx = x/50;
    else
        dx = 1/50;
    end
    
    % Find change of sign.
    twosqrt = sqrt(2);
    a = x; fa = fx; b = x; fb = fx;
    
    while (fa > 0) == (fb > 0)
        dx = twosqrt*dx;
        a = x - dx;  fa = fun(a);
        
        if (fa > 0) ~= (fb > 0) % check for different sign
            break
        end
        
        b = x + dx;  fb = fun(b);
    end
    
else
    %% intervals are provided
    assert(numel(x) == 2)
    a = x(1);
    b = x(2);
    
    fa = fun(a);
    fb = fun(b);
end


%% intervals have been found
fc = fb;
assert(sign(fa) ~= sign(fb))

if fa == 0
    b = a;
    return
end

if fb == 0
    return
end

%% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0)
        break
    end
    
    
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end;
        if p > 0, q = -q; else p = -p; end;
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
        else
            d = m;  e = m;
        end;
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler, b = b + d;
    elseif b > c, b = b - toler;
    else b = b + toler;
    end
    fb = fun(b);
end

end