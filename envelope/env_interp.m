% Compute a simple spectral envelope based on interpolation
%
% Input
%  sins        : [Hz;amp] [2xN] Spectral peaks with frequency and linear amplitudes.
%                (as provided by sin_analysis.m)
%  fs          : [Hz] Sampling frequency
%  dftlen      : DFT length of the envelope (the output will be a half spectrum !)
%  extrap_dcny : If true, extrapolate the peaks towards high frequencies in order
%                to avoid any hole just before Nyquist.
%
% Output
%  E  : Amplitude envelope (half spectrum)
%
% Copyright (c) 2011 University of Crete - Computer Science Department
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

function E = env_interp(sins, fs, dftlen, extrap_dcny, varargin)

    if nargin<5; extrap_dcny=false; end

    fks = sins(1,:);
    aks = sins(2,:);

    if nargin<4; method='linear'; end

    if extrap_dcny
        % Add freq up to Nyquist if necessary
        if fks(end)~=fs/2
            mfd = median(diff(fks));
            if mfd<1; warning('The median frequency difference between sinusoidal components is smaller than 1Hz. Something is surely wrong in the sinusoidal parameters.'); end
            while fks(end)<fs/2-mfd
                fks = [fks, fks(end)+mfd];
                aks = [aks, aks(end)];
            end
        end
    end

    % Manage behavior around DC
    % Use a few points in order to ensure symmetry at DC when using splines
    if fks(1)>0
        fks = [-fks(4:-1:1) fks];
        aks = [aks(4:-1:1) aks];
    elseif fks(1)==0
        fks = [-fks(5:-1:2) fks];
        aks = [aks(5:-1:2) aks];
    end

    % Manage Nyquist the same way as the DC
    if fks(end)<fs/2
        dftony = fs/2-fks(end:-1:end-4);
        fks = [fks fs/2+dftony]; % Add the symetrical to Nyquist
        aks = [aks aks(end:-1:end-4)];
    end

    % Do the interpolation
    bins = fs*(0:dftlen/2)'/dftlen;
    [dum idx] = unique(fks);
    fks = fks(idx);
    aks = aks(idx);
    E = interp1(fks, log(abs(aks)), bins, varargin{:});

    if 0
        figure
        plot(fks, ld(aks), 'xr');
        hold on;
        stem([0 fs/2], [-180 -180], 'xk');
        binsext = fs*(-dftlen/8:dftlen)/dftlen;
        Eext = interp1(fks, log(abs(aks)), binsext, method);
        plot(binsext, ld(exp(Eext)), 'r');
        plot(bins, ld(exp(E)), 'b');
        keyboard
    end

    E = exp(E);

return
