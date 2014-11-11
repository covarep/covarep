% Resample a feature across time
%
% Input
%
% Output
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

function Xr = irregsampling2uniformsampling(T, X, nT, fn, ifn, method, def, usemex)

    if isempty(fn) || isempty(ifn); fn=@(x)x; ifn=@(x)x; end

    if nargin<8 || isempty(usemex) || ~usemex
        interp1fn = @(x, y, xi, method, yid) interp1(x, y, xi, method, yid).';
    else
        interp1fn = @(x, y, xi, method, yid) interp1ordered(x.', y.', xi.', yid).';
    end

    % Interpolate X according to nT
    Xr = nan(length(nT), size(X,2));
    for h=1:size(X,2) % Interpolate each bin independently
        Xr(:,h) = ifn(interp1fn(T, fn(X(:,h)), nT, method, fn(def)));
    end

return
