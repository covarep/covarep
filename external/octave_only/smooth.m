% Function to do moving average filtering on a one-dimensional signal x. N should be an odd number, and otherwise it is reduced by one.

% Copyright (c) 2013 Trinity College Dublin
%
% Author 
%  John Kane kanejo@tcd.ie

function y = smooth(x,N)

    if rem(N,2)==0
    disp('smooth - reducing N by 1')
    N=N-1;
    end

    y=x;

    halfLen=(N-1)/2;
    start=halfLen+1;
    stop=length(x)-halfLen;

    for m=start:stop
    y(m) = mean(x(m-halfLen:m+halfLen));
    end

return
