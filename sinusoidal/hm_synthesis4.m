% Synthesize a harmonic signal
%
% Octave compatible
%
% Description
%  Synthesize a harmonic signal given sinusoidal parameters [1].
%  It assumes that the sinusoids are ordered and can be linked from one time
%  frame to the next.
%
% Inputs
%  frames  : N structures containing sinusoidal parameters.
%            For each frame, the sinusoid parameters are in a matrix with format:
%               [4xK] for each column: the frequency [Hz], the linear amplitude,
%               the instantaneous phase [rad] and the harmonic number of each
%               sinusoidal component.
%               The DC is ALWAYS included at the beginning of the matrix.
%
%  [wavlen]: [samples] Length of the waveform
%  [fs]    : [Hz] Sampling frequency of the synthesized waveform
%  [opt]   : Additional options (see code below)
%
% Outputs
%  syn     : The synthesized waveform
%  fs      : The sampling frequency of the synthesized waveform
%
% Example
%  Please se the HOWTO_sinusoidal example
%
% References
%  [1] G. Degottex and Y. Stylianou, "Analysis and Synthesis of Speech using an
%      Adaptive Full-band Harmonic Model," IEEE Transactions on Acoustics,
%      Speech and Language Processing, Accepted May 2013.
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

function [syn, fs, opt] = hm_synthesis4(frames, wavlen, fs, opt)

    if nargin<4
        % Options
        opt.syn_dc    = true;   % Synthesize the DC component


        opt.debug     = 1;

        opt.usemex    = 0; %[]; % Use mex, faster but use linear interpolations
                                % instead of splines.
        opt.fmapin{1} = {{'frames', 'wavlen', 'fs'}, 'sin'};
        opt.fmapout{1}  = {{'syn', 'fs'}, 'snd'};
    end
    if nargin==0; syn=opt; return; end

    if isempty(opt.usemex); opt.usemex=exist('interp1ordered')==3; end

    % ==========================================================================

    if nargin<2; wavlen=round(fs*frames(end).t)+1; end
    if nargin<3; fs=44100; end

    f0sin = [[frames.t]', [frames.f0]'];
    [f0sin, aks] = sins2ak({frames.sins}, f0sin, fs);

    defamp = log(db2mag(-300));

    T = f0sin(:,1); % Analysis times given by the provided f0 curve
    nt = length(T); % Number of anchors (and equal to number of frames)
    times = (0:wavlen-1)'/fs;

    % Estimate instantaneous fundamental frequency along the whole recording
    if opt.usemex
        f1 = interp1ordered(f0sin(:,1), f0sin(:,2), times, NaN)';
    else
        f1 = interp1(f0sin(:,1), f0sin(:,2), times, 'spline');
    end
    f1 = interp1_extrapbounds(f1); % Check if bounds are defined, and replace nan values

    % Fundamental phase from fundamental frequency
    p1 = filter(1, [1 -1], 2*pi*f1/fs);
    p1at = interp1(times, p1, T);
    p1at = interp1_extrapbounds(p1at); % T can be outside of times (const extrap is not good, but ...)

    % Count the maximum number of harmonic in the whole recording
    Hmax = 0;
    for n=1:numel(aks);
        Hmax=max(Hmax, length(aks{n})-1);
    end

    if opt.debug>0; disp(['Harmonic synthesis (Hmax=' num2str(Hmax) ', usemex=' num2str(opt.usemex) ')']); end
    syn = zeros(wavlen,1);

    if opt.syn_dc
        % Add the DC first
        % Get instantaneous amplitudes in a vector
        aa = zeros(nt,1);
        for n=1:nt
            aa(n) = real(aks{n}(1));
        end
        % Interpolate on a linear scale
        if opt.usemex
            am = interp1ordered(T, aa, times, 0)';
        else
            % Add extra values at boundaries
            Text = [times(1)-(T(2)-T(1)); T; times(end)+(T(end)-T(end-1))];
            aa = [0; aa; 0];
            am = interp1q(Text, aa, times);
            am(find(isnan(am))) = 0;
        end
        syn = syn + am;
    end
    
    if opt.debug>1
        plot(times, syn);
        keyboard
    end

    if opt.debug>0; pb = progressbar(Hmax); end
    % Then for each harmonic
    for h=1:Hmax

        % Phase
        if 1
            % Get instantaneous phase values in a vector
            % Put zero phase where there is no values. whatever the phase, the
            % amplitude is -300dB there, so it doesn't matter.
            phi = 0*ones(nt,1);
            for n=1:nt
                if h<=length(aks{n})-1
                    phi(n) = angle(aks{n}(1+h));
                end
            end
            % Get the Relative Phase using first harmonic phase
            rp = phi - h*p1at;
            % Interpolate the Relative Phases
            % This interpolation should be twice differentiable such as the first
            % derivative (the frequency) is continous.
            if opt.usemex
                % Because this interpolation is linear, step exist on the
                % frequency curve ! (see plot(diff(ph)))
                % But it doesn't seem to deteriorate the signal percetively.
                rpr = interp1ordered(T, cos(rp), times, 1)';
                rpi = interp1ordered(T, sin(rp), times, 0)';
                rp = atan2(rpi, rpr);
            else
                % Using linear, the derivative (the frequency) has steps, very bad
                % Using spline, the curve can degenerate, mainly at the bounds
                % Using pchip, the derivative (the frequency) can have strong peaks
                rp = angle(interp1(T, exp(1j*rp), times, 'spline'));
            end
            % Get the final phase from the RP and the harmonic phase
            ph = rp + h*p1;
            
%              plot(fs/(2*pi)*diff(unwrap(ph)), 'k');
%              keyboard
        else
            % Pantazis' "Adaptive AM-FM ..." Eq(30)
            ph = zeros(wavlen,1);
            for i=1:length(T)-1
                if h<=length(aks{i})-1 && h<=length(aks{i+1})-1
                    idxs = round(T(i)*fs)+1;
                    idxe = round(T(i+1)*fs)+1;
    %                  idxs-(T(i)*fs+1)
    %                  idxe-(T(i+1)*fs+1)

                    fm_hat = 2*pi/fs*h*f1(idxs:idxe);

                    % intergation of frequency
                    pm_inst = filter(1, [1 -1], fm_hat);

                    % correct phase value at nz(i)
                    pm_inst = pm_inst + repmat(angle(aks{i}(1+h))-pm_inst(1), length(pm_inst), 1);

                    % correction for the phase value at nz(i+1)
                    M = round((pm_inst(end) - angle(aks{i+1}(1+h)))/(2*pi));
                    er = pi*(pm_inst(end)-angle(aks{i+1}(1+h))-2*pi*M)/(2*(idxe-idxs));
                    t = 0:idxe-idxs;
                    ft = sin(pi*t/(idxe-idxs));
                    pm_inst = pm_inst - filter(1, [1 -1], ft'*er);

                    % save the interpolated phase
                    ph(idxs:idxe) = pm_inst;
                end
            end

%              plot(fs/(2*pi)*diff(unwrap(ph)));
%              keyboard
        end        

        % Amplitude
        % Get instantaneous amplitudes in a vector
        aa = defamp*ones(nt,1);
        for n=1:nt
            if h<=length(aks{n})-1
                aa(n) = log(abs(aks{n}(1+h)));
            end
        end
        % Interpolate on a log scale
        if opt.usemex
            am = interp1ordered(T, aa, times, defamp)';
        else
            % Add extra values at boundaries
            Text = [times(1)-(T(2)-T(1)); T; times(end)+(T(end)-T(end-1))];
            aa = [defamp; aa; defamp];
            am = interp1q(Text, aa, times);
            am(find(isnan(am))) = defamp;
        end
        am = exp(am); % Retrieve the amplitude on the linear scale

        % Keep indices only where the frequency doesn't go above the Nyquist
        idx = find(f1*h<fs/2);

        if opt.debug>1
            plot(times(idx), 2*am(idx).*cos(ph(idx)));
            keyboard
        end
        
        syn(idx) = syn(idx) + 2*am(idx).*cos(ph(idx));

        if opt.debug>0; pb = progressbar(pb,h); end
    end

    if opt.debug>1
        plot(times, syn);
        keyboard
    end

return

% Convert matrices [freq; amps; phase; #harm] to (f0,ak) pairs
function [f0s, aks] = sins2ak(parts, f0s, fs)

    aks = cell(size(f0s,1),1);

    for ind=1:size(f0s,1)
        if ~isempty(parts{ind})
            aks{ind} = parts{ind}(2,:).*exp(1j*parts{ind}(3,:));
        end
    end

return


