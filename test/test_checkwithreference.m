function test_checkwithreference(X, R, hard)
    if nargin<3; hard = true; end

    ok = true;

    if isstruct(X)
        [ok, why] = structeq(X, R);
        if ~ok
            disp(why);
        end
    else
        d = mean(mean(abs(X-R)./(max(abs(R),eps))));
        ok = d==0;
        if ~hard
            disp(['Mean relative error ' num2str(d)]);
        end
    end

    if hard & ~ok
        warning('Values differ from the reference. The method does not reproduce the same results as in previous COVAREP versions.');
    end

return
