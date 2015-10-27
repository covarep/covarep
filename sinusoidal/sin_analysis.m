% Estimate sinusoid parameters using peak picking or LS solutions
%
% Octave compatible
% 
% Description
%  Model a given waveform using sinusoids whose parameters are estimated
%  using various methods. Depending on the combinations of the options:
%  fharmonic, fquasiharm, fadapted and use_ls, the methods below are available:
%                    | fharmonic | fquasiharm | fadapted  | use_ls
%  Peak Peaking [1]     false       false        false      false
%  Peak Peaking [1]     true        false        false      false
%  Peak Peaking [1]     true        true         false      false
%  HM [2] **            true        false        false      true    (default)
%  QHM (partly [3]*)    true        true         false      true
%  aHM [5]              true        false        true       true
%  (iQHM [3] to come)
%  (aQHM [4] to come)
%
%  ** This is a full-band HM, this is NOT the Harmonic+Noise Model (HNM)
%  *  Note that [3] is NOT fully implemented. The QHM estimation is available,
%     but NOT the iterative algorithm presented in [3].
%
%  The analysis instants need to be always provided through the f0s argument.
%
%  The window length is always odd so as the sample at the middle of the
%  window correspond to the sample of the analysis instant.
%
% Inputs
%  wav    : The waveform
%  fs     : [Hz] The sampling frequency
%  f0s    : [s, Hz] [Nx2] A temporal vector with time instants and fundamental
%           frequency f0 estimated at the given time instants.
%           The sinusoid parameters are always estimated at these time instants.
%           For the Peak Peaking method, the f0 (second column) can be omitted.
%  [opt]  : Additional options (see code below)
%
% Outputs
%  frames : N structures containing the estimated sinusoid parameters and extra
%           information (e.g. window length, the f0 used).
%           For each frame, the sinusoid parameters are in a matrix with format:
%              [5xK] for each column: the frequency [Hz], the linear amplitude,
%              the instantaneous phase [rad], the harmonic number of each
%              sinusoidal component and a boolean specifying if the sinusoidal
%              parameters are from a spectral peak or through sampling.
%              The DC is ALWAYS included at the beginning of the matrix.
%  syn    : if asked, the resynthesized waveform using an Overlap-Add method.
%  opt    : The options structure which might have been altered for consistency
%           purpose.
%
% Example
%  Please se the HOWTO_sinusoidal example
%
% References
%  [1] McAulay, R., Quatieri, T.: Speech analysis/Synthesis based on a sinusoidal
%      representation, IEEE Transactions on Acoustics, Speech and Signal
%      Processing 34(4):744-754, 1986.
%  [2] Stylianou, Y.: Harmonic plus Noise Models for Speech combined with
%      Statistical Methods, for Speech and Speaker Modification, TelecomParis,
%      PhD Thesis, 1996.
%  [3] Pantazis, Y., Rosec, O., Stylianou, Y.: Iterative Estimation of Sinusoidal
%      Signal Parameters, Signal Processing Letters, IEEE 17(5):461-464, 2010.
%  [4] Pantazis, Y., Rosec, O., Stylianou, Y.: Adaptive AM-FM Signal Decomposition
%      With Application to Speech Analysis, IEEE Transactions on Audio, Speech,
%      and Language Processing 19(2):290-300, 2010.
%  [5] Degottex, G., Stylianou, Y.: Analysis and Synthesis of Speech using an
%      Adaptive Full-band Harmonic Model, IEEE Transactions on Acoustics, Speech
%      and Language Processing, 21(10):2085-2095, 2013.
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

function [frames, syn, opt] = sin_analysis(wav, fs, f0s, opt)

    if nargin<4
        % Options

        % Window
        opt.win_durf0sync = true;
        opt.win_durnbper  = 3;       % Number of period per window (def. 3)
                                     % (used only if win_durf0sync==true)
        opt.win_dur       = 30/1000; % [s] Duration of the analysis window
                                     % (used only if win_durf0sync==false)
        opt.win_fn        = @blackman; % Window type
        opt.win_dropoutside=true;    % Drop windows which are partly outside of
                                     % the signal (and thus drop the 
                                     % corresponding analysis instants).

        % Partials estimation
        opt.fharmonic = true; % Use harmonic or quasi-harmonic frequencies
        opt.fquasiharm= false;% Use quasi-harmonic frequencies
                              % It is possible to estimate the parameters using
                              % the QHM model. HOWEVER, the iterative algorithm
                              % in [3] is not implemented yet !
        opt.fadapted = false; % Adapt the frequency basis to the f0 curve
                              % (currently, works only with LS solution)

        % For LS solution [2-5]
        opt.use_ls   = true;  % Use the Least Square solution (LS) [2,3]
                              % Otherwise, use Peak Picking (PP) [1]
        opt.win_ls_marg = 0.1;% For the LS solution, add 10% to the
                              % theoretical minimum window length to ensure the
                              % stability of the LS solution

        % For Peak Picking [1]
        opt.normstr       = 'sum(win)'; % Normalize the window using the given
                                % expression.
                                % Among the available terms, the following exist
                                % win, winlen and dftlen 
        opt.dftlen        = []; % Force the DFT to a given length
                                % if empty, the DFT length is adapted to the 
                                % window length + the oversampling factor below.
                                % (ignored when using the LS solution)
        opt.osf           = 2;  % Frequency OverSampling Factor, according to:
                                %    dftlen=2^(nextpow2(winlen)+opt.osf)
        opt.frames_keepspec = false; % Keep the halfspec
        opt.funique       = true;% Keep only unique frequencies
                                 % Considered only if harmonic=false

        % Resynthesis
        opt.resyn    = false; % Do OverLap-Add (OLA) resynthesis
        
        opt.debug    = 1;

        opt.fmapin{1}  = {{'wav', 'fs'}, 'snd'};
        opt.fmapin{2}  = {{'f0s'}, 'yin'};
        opt.fmapout{1} = {{'frames'}, 'mat'};
    end
    if nargin==0; frames=opt; return; end

    if opt.debug>0; disp('Sinusoidal Analysis'); end
    if opt.fquasiharm; opt.fharmonic=true; end
    if opt.use_ls; opt.fharmonic=true; end % Currently, the LS solution is only possible with a harmonic model
%      if ~opt.use_ls && opt.fharmonic==true; opt.fquasiharm=true; end % commented: otherwise we can't do harmonic PP
    if ~opt.use_ls && opt.fquasiharm==true; opt.fharmonic=true; end
    if size(f0s,2)>1 && any(f0s(:,2)<=0); error('If a fundamental frequency is specified, it cannot be zero.'); end
    if opt.fharmonic
        if size(f0s,2)<2; error('An input f0 curve is necessary for harmonic model'); end
        opt.win_durf0sync=true;
    end
    if opt.fadapted
        if size(f0s,2)<2; error('An input f0 curve is necessary for adaptivity'); end
        if ~opt.use_ls; error('Adaptivity with peak picking is not possible. Should use LS solution.'); end
    end
    if nargout<2; opt.resyn=false; end
    if opt.debug>0; disp(opt); end

    wav = wav(:);
    if opt.resyn
        syn  = zeros(size(wav));
        wins = zeros(size(wav));
    else
        syn = [];
    end  

    T = f0s(:,1); % Get analysis time instants

    if opt.fadapted
        times = (0:length(wav)-1)'/fs;
        % Sample the first harmonic all along the signal
        f1 = interp1(f0s(:,1), f0s(:,2), times, 'spline');
        f1 = interp1_extrapbounds(f1); % Check if bounds are defined, and replace nan values
        % Compute the fundamental phase
        p1 = filter(1, [1 -1], 2*pi*f1/fs);
        % p1 = 2*pi*cumtrapz(f1)/fs; % Create Bad conditionned matrices !
    end

    % Use a constant winlen if winlen is not pitch sync
    if ~opt.win_durf0sync
        if size(f0s,2)>1
            f0m = exp(median(log(f0s(:,2))));
            winlen = round(opt.win_durnbper*fs/f0m/2)*2+1;
        else
            winlen = round(opt.win_dur*fs/2)*2+1;
        end
        win = opt.win_fn(winlen);
        if ~isempty(opt.normstr)
            eval(['d = ' opt.normstr ';']);
            win = win./d;
        end
        if ~isempty(opt.dftlen); dftlen=opt.dftlen;
        else                     dftlen=2^(nextpow2(winlen)+opt.osf); end
        W = delay2spec((winlen-1)/2, dftlen);
    end

    if opt.debug>0; pb = progressbar(length(T)); end
    for ind=1:length(T)

        % Be sure the analysis instant is on the sample of the window center
        T(ind) = round(T(ind)*fs)/fs;
        fr.t = T(ind);

        if opt.win_durf0sync
            f0 = f0s(ind,2);
            fr.f0 = f0;

            if opt.fadapted; winlen = get_optimal_winlen(f0s, fs, ind, opt);
            else             winlen = round(opt.win_durnbper*fs/f0/2)*2+1; end
            win = opt.win_fn(winlen);

            if ~isempty(opt.dftlen); dftlen=opt.dftlen;
            else                     dftlen=2^(nextpow2(winlen)+opt.osf); end
            W = delay2spec((winlen-1)/2, dftlen);
        end

        fr.winlen = winlen;
        fr.dftlen = dftlen;

        winids = -(winlen-1)/2:(winlen-1)/2; % Indices relative to the center
        idscenter = round(T(ind)*fs)+1;
        ids = idscenter + winids; % Indices of the window in the signal

        if opt.win_dropoutside
            if ids(1)<1 || ids(end)>length(wav);
                T(ind) = NaN;
                continue;
            end
            iddx = (1:winlen);
            idsb = ids;
            wavsel = wav(ids);
        else
            iddx = find(ids>=1 & ids<=length(wav)); % Valid indices of the window
            idsb = ids(iddx); % Indices of the win in the sig bounded by the sig limits
            wavsel = zeros(winlen,1);
            wavsel(iddx) = wav(idsb);
        end

        if opt.use_ls
            % Use Least Squares (LS) solution [2]
            if opt.fadapted
                % Adapt the frequency basis to the f0 curve
                Ho = floor(((fs/2)-max(f1(idsb))/2)/max(f1(idsb)));

                % Use interpolated phase (through interpolated frequencies)
                if length(idsb)<winlen;
                    % Extrapolate p1 if the window is partly outside of the sig
                    p1sel = interp1(idsb, p1(idsb), ids, 'nearest', 'extrap')';
                else
                    p1sel = p1(ids);
                end
                pm = p1sel - p1(idscenter);
                pm = pm*(-Ho:Ho);
                Nbk = size(pm,2);

                if opt.fquasiharm % Adaptive Quasi-Harmonic Model (aQHM) [4]
                    error('Adaptivity + Quasi-harmonicity (aQHM) not implemented ! ... yet');
                    
                else % Adaptive Harmonic Model (aHM) [5]
                    % Build matrices to compute the LS solution of ak
                    E = cos(pm)+1j*sin(pm); % dimension of E: (2N+1)x(Ho*2+1)
                    Ew = repmat(win,1,Nbk).*E;
                    R = Ew'*Ew;
                    fr.RCN = rcond(R); % Estimate the matrix condition number

                    x = R\(Ew'*(wavsel.*win)); % The LS solution
                    ak = x(Ho+1:end);     % amplitudes (skip the negative freqs)
                end

                % Get the f0 and the center of the window
                cf0 = f1(idscenter);
                fr.sins = [cf0*(0:Ho); abs(ak)'; angle(ak)'; (0:Ho); zeros(1,Ho+1)];

            else
                % Use stationary components
                Ho = floor(((fs/2)-f0/2)/f0);
                fk = f0*(-Ho:Ho);
                Nbk = length(fk);

                if opt.fquasiharm % Quasi-harmonics
                    pm = winids'*2*pi*fk/fs;
                    E = cos(pm)+1j*sin(pm);
                    E = [E repmat(winids',1,Nbk).*E];
                    Ew = repmat(win,1,2*Nbk).*E;
                    R = Ew'*Ew;

                    fr.RCN = rcond(R);   % Estimate the matrix condition number

                    x = R\(Ew'*(wavsel.*win));   % The LS solution
                    ak = x(Ho+1:Nbk);       % amplitudes

                else % simple-harmonics [2]
                    pm = winids'*2*pi*fk/fs;
                    E = cos(pm)+1j*sin(pm);
                    Ew = repmat(win,1,Nbk).*E;
                    R = Ew'*Ew;

                    fr.RCN = rcond(R); % Estimate the matrix condition number

                    ak = R\(Ew'*(wavsel.*win)); % The LS solution

                    ak = ak(Ho+1:end);     % amplitudes (skip the negative freqs)

                    %  y = real(E*ak); % reconstructed signal
                end

                fr.sins = [fk(Ho+1:end); abs(ak)'; angle(ak)'; (0:Ho); zeros(1,Ho+1)];
            end
        else
            % Use Peak Picking (PP) from a spectrum [1]
            if opt.win_durf0sync && ~isempty(opt.normstr)
                eval(['d = ' opt.normstr ';']);
                win = win./d;
            end
            
            % Window the signal segment
            s = wavsel.*win;

            % Compute the spectrum and compansate the window delay
            S = fft(s, dftlen).';
            S = S.*W;

            % TODO add FChT

            if opt.fharmonic
                % Select only peaks around harmonic frequencies
                fr.sins = spec_getsins_f0(S, fs, f0);

                if ~opt.fquasiharm
                    % Force harmonic frequencies after partial estimation
                    mf0 = median(diff(fr.sins(1,:)));
                    fr.sins(1,:) = mf0*(0:size(fr.sins,2)-1);
                end
            else
                % If not quasi-harmonic => unconstrained sinusoids
                fr.sins = spec_getsins(S, fs);

                if opt.funique
                    [~, idx] = unique(fr.sins(1,:));
                    fr.sins = fr.sins(:,idx);
                end
            end

            if opt.frames_keepspec; fr.S=spec2hspec(S); end
        end

        if opt.resyn
            y = sin2sig(fr.sins, fs, winlen);
            y = 2*y;
            fr.SNR = mag2db(std(wavsel)/std(wavsel-y));
            syn(idsb) = syn(idsb) + y(iddx).*win(iddx);
            wins(idsb) = wins(idsb) + win(iddx);
        end

        if ind==1
            frames(length(T)) = fr; % pre-allocate with correct fieldnames
        end
        frames(ind)=fr;

        if 0 && T(ind)>0.3
            hold off;
            plot(wavsel, 'k');
            hold on;
            plot(y, 'b');
            keyboard
        end

        if 0
            % Compute the spectrum and compansate the window delay
            % TODO Check the behavior of the harmonic structure with the aDFT.
            win = win./sum(win);
            S = fft(wavsel.*win, dftlen);
            F = fs*(0:length(S)-1)/length(S);
            hold off;
            plot(F, ld(S), 'k');
            hold on;
            stem(fr.sins(1,:), ld(fr.sins(2,:)), 'xr');
            xlim([0 fs/2]);
%              ylim([-140 40]);
%              keyboard
            pause
        end

        if 0 && opt.debug>1 && T(ind)>0.3
            V3clear();
            V3spec(S, fs, 'k');
            V3part(frames(end).sins, fs);
            keyboard
        end

        if opt.debug>0; pb = progressbar(pb, ind); end
    end
    if opt.debug>0; pb = progressbar(pb, length(T)); end

    % Drop the necessary frames
    idx = find(~isnan(T));
    if length(idx)<length(T) && opt.debug>0
        disp(['    Some windows were outside of the signal. ' num2str(length(T)-length(idx)) ' frames dropped. (use opt.win_dropoutside=false if you want to keep all windows and zero-pad the necessary ones at signal boundaries).']);
    end
    T = T(idx);
    frames = frames(idx);

    if opt.resyn
        idx = find(wins>0);
        syn(idx) = syn(idx)./wins(idx);
        if opt.debug>0; disp(['Mean SNR=' num2str(mean([frames.SNR]))]); end
    end
    
    if opt.debug>1
        Hmax = 0;
        for ind=1:numel(frames)
            Hmax = max(Hmax,size(frames(ind).sins,2));
        end
        harmstruct = NaN*ones(numel(frames), Hmax);
        for ind=1:numel(frames)
            M = min(length(frames(ind).sins(1,:)),Hmax);
            harmstruct(ind,1:M) = frames(ind).sins(1,1:M);
        end
        hold off;
        [mf0, winlen, dftlen] = spec_info(fs, f0s);
        [X, Fs, Ts] = spectrogram(wav, blackman(winlen), round(0.9*winlen), dftlen, fs);
        imagesc(Ts, Fs, lin2db(X));
        axis xy;
        hold on;
        plot([0 (length(wav)-1)/fs], 0.5*fs*[1 1], '--k');
        plot(T, harmstruct, 'k');
        keyboard
    end

return

