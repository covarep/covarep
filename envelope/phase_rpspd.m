% Compute the Phase Distortion or the Relative Phase Shift given a sinsoidal model
%
% Description
%  Based on ampltiude and phase measurements of a sinusoidal model, this function
%  compute the Phase Distortion (PD) [1-3] or the Relative Phase Shift (RPS) [4-5]
%  
%  Theoretically speaking, the PD has been shown to be related to the glottal
%  pulse shape only [1]. In practice, this property holds while both the sinsoidal
%  model and the minimum-phase estimate of the speech spectrum is accurate enough.
%
%  Phase distortion and relative phase shift are closely related. Indeed, at a
%  given time, for each harmonic k, the RPS is [4-5]:
%  
%     RPS(k) = phi(k) - k*phi(1)
%  
%  with phi being the instantaneous phase. Even though the equation in [1] is
%  rather complicate, the PD is nothing but the harmonic-by-harmonic phase
%  difference of the RPS where the minimum-phase component of the sigmal has
%  been first removed from the instantaneous phase values:
%  
%     PD(k) = diff(phitilde(k) - k*phitilde(1))
%           = diff(phitilde) - phitilde(1)
%  
%  where
%  
%     phitilde(k) = phi(k) - angle(V(k*f0))
%  
%  where V(omega) is the minimum-phase frequency response which is computed 
%  from the estimated amplitude envelope of the spectrum.
%
%  Erratum: In [3], the PD is called, by mistake, Function of Phase Distortion.
%  This latter should be actually used to describe the PD when it dependent on a
%  glottal model as a function of one of its shape parameter, as illusrated in
%  [2](p.91,Fig.6.5).
% 
%  For vizualization purpose, it is better to interpolate the PD on a frequency
%  scale in Hertz, as demonstrated in the HOWTO_envelope.m, and, thus, turning
%  the pd_harm2freq otpion to true. However, this option is false by default
%  in order to correspond to the work published in [3].
%  
% Input
%  frames : A vector of N frames structures as given by sin_analysis.m
%           (see the HOWTO_glottal_source for an example).
%  fs     : [Hz] The sampling frequency
%
% Output
%  PE  : [rad] The Phase "Envelope" representation. It is the Phase Distortion
%        (PD) [1-3] or the Relative Phase Shift (RPS) [4-5].
%        (see options pd_method in the code below to chose between the two)
%        Depending on the option pd_harm2freq, there is one value per harmonic or
%        the harmonic phase values are interpolated on a scale in Hertz.
%  AE  : The Amplitude Envelope used to compute the minimum-phase component
%        that can be removed from the phase representation.
%
% References
%  [1] G. Degottex, A. Roebel and X. Rodet, "Function of phase-distortion for
%      glottal model estimation", Proc. IEEE Int. Conf. on Acoustics, Speech, and
%      Signal Processing (ICASSP), 4608-4611, 2011.
%  [2] G. Degottex, "Glottal source and vocal-tract separation", Ph.D. thesis,
%      University Pierre and Marie Curie - Institut de Recherche et Coordination
%      Acoustique/Musique (Ircam) - CNRS-UMR9912-STMS, Paris, France, 2010.
%  [3] M. Tahon, G. Degottex and L. Devillers, "Usual voice quality features and
%      glottal features for emotional valence detection", Proc. International
%      Conference on Speech Prosody, 693-696, 2012.
%  [4] I. Saratxaga, I. Hernaez, D. Erro, E. Navas, J. Sanchez, "Simple
%      representation of signal phase for harmonic speech models", Electronics
%      Letters 45(7):381-383, 2009.
%  [5] I. Saratxaga, I. Hernaez, M. Pucher, I. Sainz, "Perceptual Importance of
%      the Phase Related Information in Speech", Proc. Interspeech, 2012.
%  [6] M. Koutsogiannaki, O. Simantiraki, G. Degottex and Y. Stylianou, "The
%      Importance of Phase on Voice Quality Assessment", In Proc. Interspeech,
%      Singapore. International Speech Communication Association (ISCA), September
%      2014.
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

function [PE, AE, opt] = phase_rpspd(frames, fs, opt)

    % Options
    if nargin<3
        % Phase
        opt.pd_vtf_rm    = true; % Remove the VTF phase from the inst. phase

        opt.dc_phase     = 0;    % if empty, does not change it; otherwise, set the given value
        opt.polarity_inv = false;% For vizualisation purpose, the phase of the
                                 % inverted signal might be more convenient.
                                 % (applied after DC's phase is set using dc_phase)
        opt.pd_method    = 1;    % 1:Phase Distortion (PD) [1-3]
                                 % 2:Relative Phase Shift (RPS) [4-5]
        opt.harm2freq    = false;% Convert the harmonic values on a hertz scale

        opt.dftlen    = 4096;    % The DFT length used for envelope estimation

        opt.usemex    = false; % Use mex fn, faster but use linear interpolation
        opt.debug     = 0;
    end
    if nargin==0; PE=opt; return; end

    AE = [];

    % A 2 column-vector for storing time and f0
    f0s = [[frames.t]', [frames.f0]'];

    % The maximum harmonic number among all frames
    Hmax = 0;
    for n=1:size(f0s,1); Hmax = max(Hmax,size(frames(n).sins,2)-1); end
    
    if opt.pd_vtf_rm

        % Estimate the amplitude envelope
        F = fs*(0:opt.dftlen/2)/opt.dftlen;      % The freq scale in Hertz
        AE = zeros(size(f0s,1),1+opt.dftlen/2); % AE will store the amplitude env
        for n=1:size(f0s,1)

            % Compute a very simple amplitude envelope
            % using linear interpolation
            E = env_interp(frames(n).sins(:,2:end), fs, opt.dftlen, true, 'linear');

            AE(n,:) = E(1:opt.dftlen/2+1);

            % If asked, remove the min-phase response of the envelope
            if opt.pd_vtf_rm
                E = spec2minphasespec(hspec2spec(E));
                E = E(1:end/2+1);
                fks = f0s(n,2)*(0:size(frames(n).sins,2)-1);
                vtfp = interp1(F, unwrap(angle(E)), fks, 'linear', 0);
                frames(n).sins(3,:) = wrap(frames(n).sins(3,:) - vtfp);
            end

            if 0
                disp(['t=' num2str(f0s(n,1)) ' i=' num2str(n)]);
                hold off;
                stem(frames(n).sins(1,2:end), mag2db(abs(frames(n).sins(2,2:end))), 'xk');
                hold on;
                F = fs*(0:opt.dftlen/2)/opt.dftlen;
                plot(F, mag2db(abs(E)), 'b');
                pause
            end
        end
    end

    % PE will store the phase information
    % by default, the missing values are random phase values in [-pi pi]
    PE = wrap((2*pi)*rand(length(f0s(:,1)), 1+Hmax));
    for n=1:size(f0s,1)

        Ks = (0:size(frames(n).sins,2)-1); % including DC

        pk = frames(n).sins(3,:); % The instantaneous phase
                        % (might be processed by the amplitude env estimation)

        % Force the DC's phase value, if given
        % By default, assume zero, since it is meaningless for acoustic sigals
        if ~isempty(opt.dc_phase);
            pk(1) = opt.dc_phase;
        end

        % If asked, reverse the signal, thus, rotate the phase half a turn
        if opt.polarity_inv;
            pk = pk + pi;
        end

        if opt.pd_method==0; % Instantaneous phase
            ppk = pk;
        elseif opt.pd_method==1;  % PD
            ppk = pk - Ks.*pk(2); % including DC
            ppk = diff(unwrap(ppk));
        elseif opt.pd_method==2;  % RPS
            ppk = pk - Ks.*pk(2); % including DC
        end

        mino = min(length(ppk),1+Hmax);
        PE(n,1:mino) = ppk(1:mino);
    end

    PE = wrap(PE);

    if opt.harm2freq
        % If asked, use a linear freq scale in Hz and not the harmonic scale
        args = {'linear', NaN};
        if opt.usemex; args{end+1} = 'usemex'; end
        PE = harmscale2hertzscale(PE, f0s, fs, opt.dftlen, opt.dc_phase, @unwrap, @wrap, args{:});
        idx = find(isnan(PE));
        PE(idx) = wrap(2*pi*rand(length(idx),1));
    end

return
