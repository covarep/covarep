% HMPD: Decompress the ampitude envelope into harmonic amplitudes
%
% This function is always used, even if the scale of the amplitude envelope is
% linear.
%
% Inputs
%  AE    : [NxD] A matrix containing the amplitude envelope.
%          D is either opt.dftlen/2+1 (from hmpd_features_compute.m), or
%          opt.amp_order+1, depending if compression is disabled or enabled.
%  f0s   : [s, Hz] [Nx2] A time/data column vector, containing the f0 curve.
%  fs    : [Hz] The sampling rate of the analyzed waveform
%  opt   : Additional options (see hmpd_features_compute.m)
%
% Outputs
%  AH    : [NxH] A marix with a row for each time instant containing the
%          harmonic amplitudes.
%          H is given by the maximum harmonic number among all frames.
%  APH   : The minimum-phase response corresponding to the sampled amplitude
%          envelope.
%
% Copyright (c) 2012 University of Crete - Computer Science Department(UOC-CSD)
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

function [AH, APH] = hmpd_amplitude_decompress(AE, f0s, fs, opt)

    if opt.usemex
        interp1fn = @(x, y, xi, yid) interp1ordered(x.', y.', xi.', yid).';
    else
        interp1fn = @(x, y, xi, yid) interp1(x, y, xi, 'linear', yid).';
    end

    Hmax = ceil(0.5*fs/min(f0s(:,2)));
    AH = log(db2mag(opt.defdbamp))*ones(size(f0s,1),1+Hmax);
    APH = zeros(size(f0s,1),1+Hmax);
    F = fs*(0:opt.enc.dftlen/2)/opt.enc.dftlen;
    for n=1:size(f0s,1)

        fks = f0s(n,2)*(0:floor(0.5*fs/f0s(n,2)));

        % Decode the envelope
        if opt.enc.amp_enc_method==0 % No encoding, it's harmonic amplitudes
            AH(n,:) = AE(n,:);
            APH(n,:) = zeros(1,1+Hmax);
        else
            if opt.enc.amp_enc_method==1     % Non-encoded envelope
                E = AE(n,:)';
            elseif opt.enc.amp_enc_method==2 % Cepstral coefficients
                if opt.enc.amp_log           % Mel-Frequency Cepstral Coefs
                    E = fwcep2hspec(AE(n,:), fs, opt.enc.dftlen, opt.enc.amp_logfn);
                else                         % Linear cepstral coefficients
                    E = cc2hspec(AE(n,:), opt.enc.dftlen);
                end
            end
            lE = hspec2minphaseloghspec(E);

            if opt.enc.amp_f0norm; lE = lE+log(sqrt(f0s(n,2))); end

            % Sample amplitude envelope with harmonics frequencies
            AH(n,1:length(fks)) = interp1fn(F, real(lE), fks, log(db2mag(opt.defdbamp)));
            APH(n,1:length(fks)) = interp1fn(F, imag(lE), fks, 0);
        end
    end
    
return
