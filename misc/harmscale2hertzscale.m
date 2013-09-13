% Resample feature across frequency, from harmonic scale to linear Hz scale
function [Xr F] = harmscale2hertzscale(X, f0s, fs, dftlen, fn, ifn, varargin)

    if nargin<6 || isempty(fn)
        fn = @(x)x;
        ifn = @(x)x;
    end

    % ... then across frequency
    F = fs*(0:dftlen/2)/dftlen;
    Xr = NaN*ones(size(f0s,1), length(F));
    for n=1:size(f0s,1)
        idx = find(~isnan(X(n,:)));
        if length(idx)>1
            Xr(n,:) = ifn(interp1(f0s(n,2)*(idx-1), fn(X(n,idx)), F, varargin{:}));
        end
%          plot(Xr(n,:));
%          pause
    end

return
