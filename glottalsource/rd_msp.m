% Rd parameter estimation of the LF glottal model using Mean Squared Phase (MSP)
%
% Description
%  Given sinusoidal parameters, this function estimates the Rd shape parameter [1]
%  of the Liljencrants-Fant (LF) glottal model [2] using the Mean Squared Phase
%  (MSP) method based on MSPD2 [3-4].
% 
% Input
%  frames : A vector of N frames structures as given by sin_analysis.m
%           (see the HOWTO_glottal_source for an example).
%  fs     : [Hz] The sampling frequency
%
% Output
%  rds    : Nx[time, Rd, conf] a 3 column vector with the analysis instant,
%           the estimated Rd value and a confidence value between 0 (lowest
%           confidence) to 1 (best confidence). This last value describes how well
%           the glottal model fits the signal.
%
% References
%  [1] G. Fant, "The LF-model revisited. Transformations and frequency domain
%      analysis", STL-QPSR 36(2-3):119-156, 1995.
%  [2] G. Fant, J. Liljencrants and Q. Lin, "A four-parameter model of glottal
%      flow", STL-QPSR, vol. 4, pp. 1-13, 1985.
%  [3] G. Degottex, A. Roebel and X. Rodet, "Phase Minimization for Glottal Model
%      Estimation", IEEE Transactions on Audio, Speech, and Language Processing
%      19(5):1080-1090, 2011.
%  [4] G. Degottex, A. Roebel and X. Rodet, "Function of phase-distortion for
%      glottal model estimation", Proc. IEEE Int. Conf. on Acoustics, Speech, and
%      Signal Processing (ICASSP), 4608-4611, 2011.
%  [5] G. Degottex, "Glottal source and vocal-tract separation", Ph.D. thesis,
%      University Pierre and Marie Curie - Institut de Recherche et Coordination
%      Acoustique/Musique (Ircam) - CNRS-UMR9912-STMS, Paris, France, 2010.
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
% TODO
%  MSPD2IX: S. Huber, A Roebel, Degottex, G.: Glottal source shape parameter
%           estimation using phase minimization variants, Proc. Interspeech, 2012.
%

function [rds] = rd_msp(frames, fs, opt)

    if nargin<3
        % Options
        opt.gmodel = 1;        % 1:LF; 2:ALM
        opt.search_method = 2; % 1:Grid search; 2:Brent
        opt.debug = 0;
    end
    if nargin==0; rds=opt; return; end

    Rdrg = 0.3:0.1:2.5;

    % Allocate the output: [Time, Rd, confidence factor]
    rds = nan(length(frames), 3);
    rds(:,1) = [frames.t];

    pb = progressbar(numel(frames));
    for n=1:numel(frames)

        shmopt = sin2shm();
        shmopt.nbh = size(frames(n).sins,2);

        f0 = frames(n).f0;

        M = sin2shm(frames(n).sins, shmopt);
        M = M(1:end/2+1).';

        if length(M)<9;
            warning('The order of the harmonic model is too low for a reliable Rd estimate (Is the f0 estimate too high for the used sampling frequency?).');
        end

        if opt.search_method==1
            % Grid search
            freqs = f0*(0:length(M)-1);
            E = ones(length(Rdrg),1);
            for Rdi=1:length(Rdrg)
                E(Rdi) = optimfnRd(Rdrg(Rdi), fs, f0, M, opt);
            end
            [err, mini] = min(E);
            rds(n,2) = Rdrg(mini);
            rds(n,3) = 1 - sqrt(err)/pi; % [5](5.9)

            if opt.debug
                subplot(2,1,1);
                    [s, t] = sin2sig(frames(n).sins, fs, 2*round(3*fs/f0/2)+1);
                    plot(t, s, 'k');
                    title(['t=' num2str(frames(n).t)]);
                subplot(2,1,2);
                    plot(Rdrg, E);
                pause
            end

        elseif opt.search_method==2
            % Brent's search
            Options.Display = 'off';
            Options.MaxIter = 32;
            Options.TolX = 0.01;
            [Rd, err] = fminbnd(@(x)optimfnRd(x, fs, f0, M, opt), min(Rdrg), max(Rdrg), Options);
            rds(n,2) = Rd;
            rds(n,3) = 1 - sqrt(err)/pi; % [5](5.9)
        end

        pb = progressbar(pb, n);
    end

return

function err = optimfnRd(Rd, fs, f0, M, opt)
    
    [te, tp, ta] = Rd2tetpta(Rd);

    % Compute the glottal model frequency response
    G = gfm_spec_lf(f0*(1:length(M)-1), fs, 1/f0, 1, te, tp, ta);
    G = [eps G]; % Add a fake DC (replaced below anyway)

    % Compute the residual
    N = M./G.';
    N(end) = abs(N(end-1));
    N(1) = abs(N(2));
    N = hspec2spec(N);
    R = N./spec2minphasespec(N(:));

    % Extract the phase
    o = min([length(R)-1, 8]);
    P = angle(R(1:1+o));

    % Compute the MSPD2 error according to [3](last eq. p.5)  or [4](5)
    PD2I = wrap(diff(P) - P(2));  % Equal to [3](last eq. p.5) (and a lot simpler)
    err = mean(PD2I(2:o).^2);

return
