
function x = interp1_extrapbounds(x)

    % Check if bounds are defined, and replace nan values
    if isnan(x(1))
        firstnonnan = min(find(~isnan(x)));
        x(1:firstnonnan-1) = x(firstnonnan);
    end
    if isnan(x(end))
        lastnonnan = max(find(~isnan(x)));
        x(lastnonnan+1:end) = x(lastnonnan);
    end

return

