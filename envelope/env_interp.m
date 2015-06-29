% Compute a spectral envelope based on a simple interpolation
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

    if nargin<4 || isempty(extrap_dcny); extrap_dcny=0; end

    if strcmp(varargin{1}, 'interp1fn')
        interp1fn = varargin{2};
        varargin = varargin(3:end);
    else
        interp1fn = @interp1;
    end

    fks = sins(1,:);
    aks = sins(2,:);

    if sum(extrap_dcny)
        if length(extrap_dcny)==1 || (length(extrap_dcny)==2 && extrap_dcny(1))
            if fks(1)~=0 % If there is no DC, add one ...
                a0 = aks(1);
                % a0 = af(mi).a(1)+(af(mi).a(1)-af(mi).a(2));
                fks = [0, fks];
                aks = [a0, aks];
            else                  % If there is a DC, put it to |H1|
                a0 = aks(2);
                % a0 = af(mi).a(1)+(af(mi).a(1)-af(mi).a(2));
                aks(1) = a0;
            end
        end

        if length(extrap_dcny)==1 || (length(extrap_dcny)==2 && extrap_dcny(2))
            % Add freq up to Nyquist if necessary
            if fks(end)~=fs/2
                mfd = median(diff(fks));
                if mfd<0.01; warning('The median frequency difference between sinusoidal components is smaller than 1Hz. Something is surely wrong in the sinusoidal parameters.'); end
                while fks(end)<fs/2-mfd
                    fks = [fks, fks(end)+mfd];
                    aks = [aks, aks(end)];
                end
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
    elseif fks(end)==fs/2
        dftony = fs/2-fks(end-1:-1:end-5);
        fks = [fks fs/2+dftony]; % Add the symetrical to Nyquist
        aks = [aks aks(end-1:-1:end-5)];
    end

    % Do the interpolation
    bins = fs*(0:dftlen/2)'/dftlen;
    [dum idx] = unique(fks);
    fks = fks(idx);
    aks = aks(idx);
    E = interp1fn(fks, log(abs(aks)), bins, varargin{:});

    E = exp(E);

    if 0
        hold off;
        plot(fks, ld(aks), 'xr');
        hold on;
        stem([0 fs/2], [-180 -180], 'xk');
        binsext = fs*(-dftlen/8:dftlen)/dftlen;
        Eext = interp1(fks, log(abs(aks)), binsext, varargin{:});
        plot(binsext, ld(exp(Eext)), 'r');
        plot(bins, ld(E), 'b');
        xlim([0 fs/2]);
        ylim([-250 -80]);
        pause
%          keyboard
    end

return
