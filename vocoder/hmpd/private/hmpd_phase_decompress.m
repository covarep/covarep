% HMPD: Create a Relative Phase Shift corresponding to given PDM and PDD
%  
% Inputs
%  PDM   : [NxD] A matrix containing the Phase Distortion Mean.
%          D is either opt.dftlen/2+1 or opt.pdm_order+1, depending if
%          compression is disabled or enabled.
%  PDD   : [NxD] A matrix containing the Phase Distortion Deviation.
%          D is either opt.dftlen/2+1 or opt.pdd_order+1, depending if
%          compression is disabled or enabled.
%  f0s   : [s, Hz] [Nx2] A time/data column vector, containing the
%          analysis instants and the f0 curve. 
%  fs    : [Hz] The sampling rate of the analyzed waveform
%  opt   : Additional options (see hmpd_features_compute.m)
%  
% Outputs
%  RPS   : A Relative Phase Shift (RPS) which can be used for synthesis.
%
% Copyright (c) 2013 University of Crete - Computer Science Department(UOC-CSD)/ 
%                    Foundation for Research and Technology-Hellas - Institute
%                    of Computer Science (FORTH-ICS)
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

function [RPS] = hmpd_phase_decompress(PDM, PDD, f0s, fs, opt)

    Hmax = ceil(0.5*fs/min(f0s(:,2)));

    % Phase info will be decoded into a Relative Phase Shift (RPS),
    % (phase without the linear phase and without VTF phase)
    RPS = zeros(size(f0s,1),1+Hmax);
    for n=1:size(f0s,1)

        % PDM, the mean (pulse's shape)
        if isempty(PDM)
            a = zeros(1,opt.enc.dftlen);
        else
            % If PDM is in linear-freq scale, put it back to a harmonic scale
            if opt.enc.pdm_log
                % Compressed PDM has not been tested for synthesis
                error('Compressed PDM not supported for synthesis');
            else
                F = fs*(0:opt.enc.dftlen/2)/opt.enc.dftlen;
                fks = f0s(n,2)*(1:floor(0.5*fs/f0s(n,2)));
                if opt.usemex
                    a = wrap(interp1ordered(F.', unwrap(PDM(n,:)).', fks.', 0));
                else
                    a = wrap(interp1(F, unwrap(PDM(n,:)), fks, 'linear'));
                end
            end
        end

        % PDD, the variance (phase randomness)
        if isempty(PDD)
            % If PDD is not provided, do not add anything.
            % (randomness of unvoiced segments still depend on opt.rps_randunvoiced)
            hmin = size(RPS,2);
        else
            if opt.enc.pdd_log % Decompress the variance
                PDDv = fwcep2hspec(PDD(n,:), fs, opt.enc.dftlen);
                PDDv = abs(PDDv);
            else
                PDDv = PDD(n,:);
            end

            F = fs*(0:length(PDDv)-1)/(2*(length(PDDv)-1));
            fks = f0s(n,2)*(1:floor(0.5*fs/f0s(n,2)));
            PDDv = max(PDDv, 10^-100);
            if opt.usemex
                av = exp(interp1ordered(F.', log(PDDv).', fks.', log(2)));
            else
                av = exp(interp1(F, log(PDDv), fks, 'linear'));
            end

            av = hmpd_phase_pdd_correction(av, opt.pdv_corrthresh);

            hmin = min([length(a), length(av), size(RPS,2)]);

            a(1:hmin) = a(1:hmin) + wrappednormrnd(0, av(1:hmin).', hmin, true).';

            a(hmin+1:end) = a(hmin+1:end) + rand(1,length(a)-hmin);
        end

        a = wrap(a);

        % Retreive the corresponding RPS
        rp = wrap(cumsum(a)).';
        rp = [0; rp(1:end-1)]; % Add the DC phase (lost through diff)
        RPS(n,1:hmin) = rp(1:hmin);
        RPS(n,hmin+1:end) = 2*pi*rand(1,size(RPS,2)-hmin); % Fill with noise
        RPS(n,:) = wrap(RPS(n,:));
    end

return
