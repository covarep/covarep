% Interpolate a feature wich is on a harmonic scale into a lin freq scale in Hertz
%
% Description
%  See end of the code of phase_rpspd.m for an example.
%
% Input
%  X          : A matrix containing the harmonic data to interpolate.
%  f0s        : [time,f0] The f0 corresponding to the harmonic data
%               at each analysis time.
%  fs         : [Hz] The sampling frequency of the frequency scale.
%  dftlen     : DFT length of the envelope (the output will be half spectrums !).
%  [fn]       : A function converting the data into a better suited scale.
%  [ifn]      : The inverse function of fn.
%  [varargin] : Extra arguments to provide to the interpolation function.
%
% Output
%  Xr         : Interpolated values on a linear frequency scale in Hertz
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

function [Xr F] = harmscale2hertzscale(X, f0s, fs, dftlen, dcval, fn, ifn, varargin)

    if strcmp(varargin{1},'linear') && any(strcmp(varargin, 'usemex'))
        varargin = varargin(~strcmp(varargin, 'usemex'));
        if numel(varargin)<2; error('For the mex version, varargin needs at least two elements, the second one being the default value'); end
        interp1fn = @(x, y, xi, varargin) interp1ordered(x.', y.', xi.', varargin{2});
    else
        interp1fn = @(x, y, xi, varargin) interp1(x, y, xi, varargin{:});
    end

    if nargin<5
        dcval = [];
    end
    if nargin<6 || isempty(fn)
        fn = @(x)x;
        ifn = @(x)x;
    end

    F = fs*(0:dftlen/2)/dftlen;
    Xr = nan(size(f0s,1), length(F));
    for n=1:size(f0s,1)
        idx = find(~isnan(X(n,:)));
        if length(idx)>1
            if isempty(dcval)
                Xr(n,:) = ifn(interp1fn(f0s(n,2)*(idx-1), fn(X(n,idx)), F, varargin{:}));
            else
                Xr(n,:) = ifn(interp1fn([0 f0s(n,2)*0.95 f0s(n,2)*idx], fn([dcval dcval X(n,idx)]), F, varargin{:}));
            end
        end
    end

return
