% 
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

% Estimate the paramaters of a wrapped normal distribution
% http://en.wikipedia.org/wiki/Wrapped_normal_distribution
function [m s] = wrappednormestim(angles, dim)

    if nargin<2; dim=1; end

    expjangle = exp(1j*angles);

    z = mean(expjangle, dim);

    m = angle(z);

    s = sqrt(-2*log(abs(z)));

    if nargout==0
        hx = -pi:0.1:pi;
        hy = hist(angles, hx);
        hy = hy./sum(hy);

        hold off;
        plot(hx, (hy), 'k');
        hold on;
        p = wrappednormpdf(hx, m, s);
        p = p./sum(p);
        plot(hx, (p), 'b');
        ylim([0 max(p)+0.1]);
        title(['m=' num2str(m) ' s=' num2str(s)]);
    end
    
return
