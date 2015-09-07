% Spectrum of the Liljencrants-Fant (LF) glottal flow derivative model
%
% Description
%  This function generates the spectrum of the LF glottal pulse model of the
%  glottal flow derivate. See [1] for the origianl publication of the LF model
%  and [2] for the analytical expressions in the spectral domain.
%   
% Input
%  f   : The frequency values where the spectrum has to be estimated.
%  fs  : The sampling frequency (used for proper normalization of the
%        pulse amplitude).
%  T0  : [s] The fundamental period (the duration of the glottal cycle), 1/f0.
%  Ee  : The amplitude of the pulse (e.g. 1).
%  te  : Glottal shape parameter, see [1]p.6 (assuming T0=1 ! (T0-normalized))
%  tp  : Glottal shape parameter, see [1]p.6 (assuming T0=1 ! (T0-normalized))
%  ta  : Glottal shape parameter, see [1]p.6 (assuming T0=1 ! (T0-normalized))
%
% Output
%  G   : LF spectrum model of the  Glottal Flow Derivative.
%        Same length as the f parameter
%
% References
%  [1] G. Fant, J. Liljencrants and Q. Lin, "A four-parameter model of glottal
%      flow", STL-QPSR, vol. 4, pp. 1-13, 1985.
%  [2] B. Doval and C. d'Alessandro, "Spectral correlates of glottal waveform
%      models: an analytic study", ICASSP, 1997.
%      Erratum: There is a missing parenthesis in (1), please see the code below
%               for the proper equation.
%  [3] B. Doval, C. d'Alessandro and N. Henrich, "The spectrum of glottal flow
%      models", Acta acustica united with acustica, 92(6), 1026-1046, 2006.
%
% Copyright (c) 2008 University of Crete - Computer Science Department, SigProcLab
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
%  Georgos P. Kafentzis <kafentz@csd.uoc.gr>
%  Gilles Degottex <degottex@csd.uoc.gr>
%
% TODO
%  Stabilize the numerical computation, it doesn't work for Rd<0.3 or Rd>2.5
%  Make the computation independent from T0 (for better numerical stability)
%

function G = gfm_spec_lf(f, fs, T0, Ee, te, tp, ta)

    Te = te*T0;
    Tp = tp*T0;
    Ta = ta*T0;

    wg = pi/Tp;                                     % [1](2)

    % e is expressed by an implicit equation
    fb = @(e) 1 - exp(-e.*(T0-Te)) - e.*Ta;         % [1](12) (or [3](p.18) )
    e = fzero(fb,1/(Ta+eps));

    % a is expressed by another implicit equation   % based on [3]p.18
    % integral{0, T0} ULF(t) dt, where ULF(t) is the LF model equation
    A = (1-exp(-e*(T0-Te)))/(e^2*Ta) - (T0-Te)*exp(-e*(T0-Te))/(e*Ta);
    fa = @(a) (a.^2+wg^2)*sin(wg*Te)*A + wg*exp(-a.*Te) + a.*sin(wg*Te) - wg*cos(wg*Te);
    % find smaller interval than [0, 1e9]
    x = [0 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9];
    idx = find(diff(sign(fa(x))),1);
    a = fzero(fa, x([idx idx+1]));

    % E0 parameter
    E0 = -Ee/(exp(a*Te)*sin(wg*Te));                % [1](5)

    % LF spectrum formula                           % [2](1)
    P1 = E0*(1./((a - 1i*2*pi.*f).^2 + wg^2));
    P2 = (wg + exp((a - 1i*2*pi.*f)*Te).*((a - 1i*2*pi.*f)*sin(wg*Te) - wg*cos(wg*Te)));
    P3 = (Ee*(exp( - 1i*2*pi.*f*Te)./((e*Ta*1i*2*pi*f).*(e + 1i*2*pi.*f))));
    P4 = (e*(1 - e*Ta).*(1 - exp(- 1i*2*pi.*f*(T0 - Te))) - e*Ta*1i*2*pi.*f);
    G = P1.*P2 + P3.*P4;

    % Fix the amplitude so as Ee is respected
    G = fs*G; % The max frequency is assumed to be the band limit

    if 0
        figure
        subplot(211);
        G(1) = 0;
        G = hspec2spec(G);
        g = real(ifft(G));

        ts = (0:length(g)-1)/fs;
        hold off;
        plot(ts, g, 'k');
        hold on;
        stem(Te, min(g), 'xb');
        stem(Tp, 0, 'xb');
        grid on;

        subplot(212);
        GI = G./spec_derivative(length(G),1);
        GI(1) = abs(GI(2));
        GI(end/2+1) = abs(GI(end/2));
        g = real(ifft(GI));
        ts = (0:length(g)-1)/fs;
        hold off;
        plot(ts, g, 'k');

        keyboard
    %      pause
    end
end
