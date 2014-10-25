% Harmonic Model + Phase Distortion (HMPD)
%
% Copyright (c) 2012 University of Crete - Computer Science Department
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

function [AEC, PDMC, PDDC] = hmpd_features_compress(f0s, AE, PDM, PDD, fs, opt)

    % Amplitude
    if opt.amp_enc_method==1
        AEC = AE;
    else
        AEC = zeros(size(AE,1), 1+opt.amp_order);
        for n=1:size(AE,1)
            if opt.amp_enc_method==2 % Cepstral coefficients
                if opt.amp_log       % With mel frequency scale
                    AEC(n,:) = hspec2fwcep(AE(n,:), fs, opt.amp_order, opt.amp_logfn);
                else
                    cc = spec2rcc(hspec2spec(AE(n,:).'))';
                    AEC(n,:) = [cc(1:min(length(cc),1+opt.amp_order)); zeros(1+opt.amp_order-length(cc),1)];
                end
            end
        end
    end


    % Phase Distortion Mean (PDM)
    if ~opt.pdm_log
        PDMC = PDM;
    else
        PDMC = zeros(size(PDM,1), 1+opt.pdm_order);
        % Now, the phase is on a harmonic scale, with varying size cross time.
        % If asked, compress it using an f0 relative log scale with fixed size.
        for n=1:size(PDM,1)
            PDMC(n,:) = philin2philog(PDM(n,:), opt.pdm_log_hb, opt.dftlen/2, 1+opt.pdm_order);

            if 0
                hold off;
                plot(1:length(PDM(n,:)), PDM(n,:), 'k');
                hold on;
                % plot(1:length(PEMC(n,:)), PEMC(n,:), 'b');
                phidec = philog2philin(PDMC(n,:), opt.pdm_log_hb, opt.dftlen/2, 1+opt.pdm_order);
                plot(1:length(phidec), phidec, 'r');
                xlim([0, 50]);
                ylim([-pi, pi]);
                pause
                % keyboard
            end
        end
    end


    % Phase Distortion Deviation (PDD)
    if ~opt.pdd_log
        PDDC = PDD;
    else
        PDDC = zeros(size(PDD,1), 1+opt.pdd_order);
        for n=1:size(PDD,1)
            X = PDD(n,:);
            X(1) = X(2);
            idx = find(X==0);
            X(idx) = 0.001;
            PDDC(n,1:1+opt.pdd_order) = hspec2fwcep(X, fs, opt.pdd_order);
            if 0
                B = abs(mfcc2spec(mfcc, fs, opt.dftlen));
                hold off;
                F = fs*(0:opt.dftlen/2)/opt.dftlen;
                plot(F, X(1:end/2+1), 'k');
                hold on;
                plot(F, B(1:end/2+1), 'b');
                ylim([0 4]);
                pause
            % PEV(n,:) = B(1:end/2+1);
            end
        end
    end

return
