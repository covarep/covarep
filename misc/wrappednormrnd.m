% This is an approximation of random numbers wrt to wrapped Normal Distrib
%
% Any neater, faster, more precise and accurate generator is very welcome.
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

function wgr = wrappednormrnd(mu, sigma, n, pregen)

    if nargin<4; pregen=false; end

    % The CDF is symetric around the mean. Thus, just need zero-mean random values
    % of unitary std and scale them wrt sigma.
    wgr = wrappednormrndunitystd(n, pregen);

    wgr = sigma.*wgr;

    wgr = wrap(wgr + mu);

return


function wgr = wrappednormrndunitystd(n, pregen)

    if nargin<2; pregen=0; end

    if pregen==1
        global wrappednormrnd;

        % If run for the first time, initialize the cumulative distrib fns
        if isempty(wrappednormrnd) || ~isfield(wrappednormrnd, 'icdf')
            fprintf('wrappednormrnd.m: Pre-compute the cumulative distribution function once for all ... ');

            xs = -pi:(2*pi)/1e6:pi;
            xs = xs(1:end-1);
            cdf = wrappednormcdf(xs, 0, 1);
            % Resample the inverse CDF with 1 million uniform steps
            % This also assumes that we don't need a thinner resolution of the distribution
            wrappednormrnd.icdf_ys = (0:1e-4:1);
            wrappednormrnd.icdf = interp1(cdf, xs, wrappednormrnd.icdf_ys, 'linear', 'extrap').';

            fprintf('done.\n');
            if nargin==0; return; end
        end

        % wgr = interp1(wrappednormrnd.icdf_ys, sigma*wrappednormrnd.icdf, rand(n,1), 'linear', 'extrap');
        % The following is a speed-up version of the above line.
        ri = (length(wrappednormrnd.icdf)-1)*rand(n,1)+1;
        rfi = floor(ri);
        rci = ceil(ri);
        wgr = wrappednormrnd.icdf(rfi);
        wgr = wgr + (ri-rfi).*(wrappednormrnd.icdf(rci)-wgr);

        % % Or a faster sampling (implies only 1e4 different random values can be generated).
        % idx = round(length(wrappednormrnd.icdf)*rand(n,1))+1;
        % idx = max(1,min(length(wrappednormrnd.icdf),idx));
        % wgr = wrappednormrnd.icdf(idx);
    else
        % The non-optimized random generator
        xs = -pi:(2*pi)/1e6:pi;
        xs = xs(1:end-1);
        cdf = wrappednormcdf(xs, 0, 1);

        wgr = interp1(cdf, xs, rand(n,1), 'linear', 'extrap');
    end

return
