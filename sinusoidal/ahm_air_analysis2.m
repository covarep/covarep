% aHM-AIR iterative estimation method
%  
% Octave compatible
% 
% Description
%  Refine a fundamental frequency curve using the Adaptive Iterative Refinement
%  method and the Adaptive Harmonic Model (aHM) [1][2].
%
% Inputs
%  wav   : [Nx1] The waveform
%  fs    : [Hz] Sampling frequency
%  f0sin : [s,Hz] [Nx2] A temporal-data vector containing a rough estimate
%          of the fundamental frequency and their corresponding analysis instants.
%  opt   : Additionnal options (see code below)
%
% Outputs
%  f0sout : [s,Hz] [Mx2] The temporal-data vector containing the refined
%           fundamental frequency estimate.
%  frames : N structures containing the estimated sinusoid parameters and extra
%           information (e.g. window length, the f0 used).
%           For each frame, the sinusoid parameters are in a matrix with format:
%              [4xK] for each column: the frequency [Hz], the linear amplitude,
%              the phase [rad] and the harmonic number of each sinusoidal
%              component.
%              The DC is ALWAYS included at the beginning of the matrix.
%
% Example
%  Please se the HOWTO_sinusoidal example
%
% References
%  [1] G. Degottex and Y. Stylianou, “Analysis and Synthesis of Speech using an
%      Adaptive Full-band Harmonic Model,” IEEE Transactions on Acoustics,
%      Speech and Language Processing, Accepted May 2013.
%  [2] G. Degottex and Y. Stylianou, "A Full-Band Adaptive Harmonic
%      Representation of Speech," in Proc. Interspeech, ISCA, 2012.
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

function [f0sout, frames, opt] = ahm_air_analysis2(wav, fs, f0sin, opt)

    if nargin<4
        % Options
        opt.win_durnbper = 3;% Number of periods in the window
        opt.win_ls_marg= 0.1;% For LS, add 10% to the theoretical minimum size
        opt.win_fn    = @blackman; % Window function
        opt.use_f0time= false; % If true, use the time instants given in
                               % the f0sin argument.
        opt.nbat      = 1;   % Number of analysis time per period (def. 1)
                             % Use 1, otherwise the params mean nothing
                             % (the params accuracy and precision degrades
                             %  substantially)
                             % 1 or more makes the quality of the resynthesis
                             % almost perfect.

        opt.maxit   = 24;    % Maximum number of iterations
        opt.Kimin   = 4;     % Minimal starting number of harmonic
        opt.Kimax   = 12;    % Maximal starting number of harmonic
        opt.f0min   = 50;    % [Hz] Lower bound of the f0
        opt.f0max   = 2000;  % [Hz] Higher bound of the f0
        opt.df0     = 20;    % [Hz] maximum error assumed on the f0 curve
        opt.f0_prec = 0.1;   % adapt while df0*Hmax<prec*f0 (max is 0.5)
        opt.snr_prec= 0.1;   % [dB] Ready to stop if the SNR is not improving
                             % better than snr_prec

        opt.do_final_ahm_step = true;

        opt.debug = 1;
        
        opt.use_ls = true;    % Do not change
        opt.fquasiharm = true;% Do not change

        opt.fmapin{1}  = {{'wav', 'fs'}, 'snd'};
        opt.fmapin{2}  = {{'f0sin'}, 'yin'};
        opt.fmapout{1} = {{'f0sout', 'frames'}, 'AIR-aHM'};
    end
    if nargin==0; f0sout=opt; return; end

    % ==========================================================================

    disp(['aHM-AIR analysis']);
    if opt.debug>0; disp(opt); end

    if any(f0sin(:,2)<=0); error('The fundamental frequency cannot be zero. A value has to be specified at all time instants. In noisy segments, a consant value can be used, or an interpolation bewteen the surrounding quasi-period segments.'); end
    f0sin = f0sin(:,1:2); % Drop possible extra info
    % Clip f0
    opt.f0max = min(opt.f0max, fs/2);
    f0sin(:,2) = max(opt.f0min*ones(size(f0sin,1),1), f0sin(:,2));
    f0sin(:,2) = min(opt.f0max*ones(size(f0sin,1),1), f0sin(:,2));

    disp('    Estimate main-lobe bandwidth w.r.t. winlen');
    maxwl = ceil(opt.win_durnbper*fs/min(f0sin(:,2)));
    wls = 7:(sqrt(maxwl)-3)/100:sqrt(maxwl);
    wls = round(wls.^2);
    Bws = zeros(size(wls));
    for i=1:length(wls)
        Bws(i) = mainlobebw(opt.win_fn(wls(i)).^2, fs);
    end

    % Generate analysis time instants
    tmargin = 0.51*opt.win_durnbper*1/opt.f0min;
    if opt.use_f0time; atmethod = 1;
    else               atmethod = 2; end
    T = gen_analysis_times(wav, fs, tmargin, true, atmethod, f0sin, opt.nbat);
    f0sin = interp1td(f0sin, T);
    nt = length(T);   % Number of anchors (and equal to number of frames)
    f0s = f0sin(:,2); % (shortcut)
    times = (0:length(wav)-1)/fs;

    disp('    Compute window length and initial frequency limit');
    df0s = Inf*ones(nt,1);
    f0s_new = f0s;
    SNRs = zeros(nt,opt.maxit);
    frames = [];
    aks = cell(nt,1);
    for ind=1:nt
        fr.t = T(ind);

        fr.winlen = get_optimal_winlen([T, f0s], fs, ind, opt);

        fr.Hmax = floor((0.5*fs-0.5*f0s(ind))/f0s(ind));

        % Estimate Main-lobe bandwidth
        Bw = interp1(wls, Bws, fr.winlen);
        Bw = min(Bw/3,0.5*f0s(ind));

        % Starting harmonic number
        % TODO adapt to the freq where the interferances can start
        H = floor(Bw/opt.df0);
        H = min(opt.Kimax, H); % first step shouldn't be higher
        H = max(opt.Kimin, H);
        fr.H = min(H, fr.Hmax);

        fr.RCN = 1;
        fr.f0 = []; % Filled at the end
        fr.sin = [];% Filled at the end

        if length(frames)==0; frames = fr;
        else                  frames(end+1)=fr; end
    end

    aks_new = aks;
    f0s_new = f0s;
    toups = 1:nt;
    isfinished=zeros(nt,1); hasconverged=zeros(nt,1); isimproving=ones(nt,1); isdiverging=zeros(nt,1);
    U = eye(nt);
    adapt = 1;
    while length(toups)>0 && adapt<=opt.maxit

        % Sample the first harmonic all along the recording
        f1 = interp1(T, f0s, times, 'linear')';
        f1 = interp1_extrapbounds(f1);
        p1 = filter(1, [1 -1], 2*pi*f1/fs);

        disp(['    Adaptation ' num2str(adapt) ': ' num2str(length(toups)) '/' num2str(nt) ' frames to update']);
        if opt.debug>0; pb = progressbar(nt); end
        for toupsi=1:length(toups)
            ind = toups(toupsi);

            % Update window size here because f0 may have changed
            frames(ind).winlen = get_optimal_winlen([T, f0s], fs, ind, opt);

            % Generate the window
            win = opt.win_fn(frames(ind).winlen);

            % Get indices of the window (like spectrogram.m)
            winids = -(frames(ind).winlen-1)/2:(frames(ind).winlen-1)/2;
            ids = round(T(ind)*fs)+1 + winids;

            if ids(1)<1 || ids(end)>length(wav); error('Window out of signal'); end

            s = wav(ids);

            % Update the frequency vector according to the new f0
            if length(aks{ind})>1+frames(ind).Hmax
                aks{ind} = aks{ind}(1:1+frames(ind).Hmax);
            elseif length(aks{ind})<1+frames(ind).Hmax
                aks{ind} = [aks{ind}, zeros(1,(frames(ind).Hmax+1)-length(aks{ind}))];
            end

            % Interpolate parameters

            % Get the anchors indices covering the window
            % Take first the indices strictly inside the window ...
            is = find(T>=(T(ind)-(frames(ind).winlen-1)/(fs*2)) & T<=(T(ind)+(frames(ind).winlen-1)/(fs*2)));
            % ... then add one extra on each side to fully cover the window
            is = [is(1)-1; is; is(end)+1];
            is = is(find(is>=1 & is<=nt));

            U(ind,:)=0; U(ind,is)=1;

            Ho = floor(((fs/2)-max(f1(ids))/2)/max(f1(ids)));

            % pm = cumtrapz(fm); % Create Bad conditionned matrices !
            pm = p1(ids) - p1(ids((frames(ind).winlen-1)/2+1));
            pm = pm*(-Ho:Ho); % Purely harmonic
            Nbk = size(pm,2);

            % Build matrices to compute the LS solution of ak bk
            E = cos(pm)+1j*sin(pm); % dimension of E: (2N+1)x(Ho*2+1)
            E = [E repmat(winids',1,Nbk).*E]; % dimension of E:(2N+1)x(2Ho2+1)
            Ew = repmat(win,1,2*Nbk).*E;
            R = Ew'*Ew;
            frames(ind).RCN = rcond(R); % Compute the matrix condition number

            x = R\(Ew'*(win.*s)); % The LS solution
            ak = x(Ho+1:Nbk);     % amplitudes
            bk = x(end-Ho:end);   % slopes
            aks_new{ind} = ak.';

            % Reconstruct the signal with ak only (bk is never used in synth)
            % y = real(E*[ak;bk]);
            y = 2*real(E(:,Ho+1:Nbk)*ak); % The synthesis will be of that kind
            SNRs(ind,adapt) = mag2db(std(s)/std(s-y));

            % Compute frequencies bias
            df = fs/(2*pi)*((imag(bk).*real(ak) - imag(ak).*real(bk))./abs(ak).^2)';
            Ks = (0:length(ak)-1);

            % Check frequency corrections for update of f0
            % Consider only df which are consistent
            iscons = ones(size(df));
            iscons(1) = 0;
            iscons = iscons & ~isnan(df) & ~isinf(df);
            iscons(min(Ho,frames(ind).H)+1+1:end) = 0;

            % Ignore any frequency mismatch which makes the frequency track
            % diverging too much from the harmonic position.
            % Not so good according to quantitative evaluation. The estimated
            % df can be meaningfull even though abs(df)>f0/2
            % However, good for resynthesis because if the f0 curve is not
            % constrained enough it can degenerate in unvoiced segments
            % (can make strong low sin in fricatives or clics in transients).
            iscons = iscons & (abs(df) < f0s(ind)/2);

            % According to Pantazis "Adapt AM-FM [...]" |phi1|/2pi < B/3 (24)
            Bw = interp1(wls, Bws, frames(ind).winlen);
            iscons = iscons & (abs(df) < Bw/3);

            % Avoid unrealistic small frequencies (i.e. and thus negative freq)
            iscons = iscons & (Ks*f0s(ind)+df > opt.f0min);

            idx = find(iscons);
            if length(idx)==0;  df0s(ind) = 0;
            else                df0s(ind) = median(df(idx)./Ks(idx));     end

            % Take into account the correction
            f0s_new(ind) = f0s(ind) + df0s(ind);% '+' because df is a missmatch
            % and clip f0
            f0s_new(ind) = max(opt.f0min, min(opt.f0max,f0s_new(ind)));

            % Compute step size from window bandwidth and current error
            % Update current harmonic level according to df
            if frames(ind).H<frames(ind).Hmax

                % Update the number of harmonics in the frame
                frames(ind).Hmax = floor((0.5*fs-0.5*f0s_new(ind))/f0s_new(ind));

                % Update the level where harmonics are updated
                Hprev = frames(ind).H;
                if adapt==opt.maxit-1
                    % If not finished, force it
                    frames(ind).H = frames(ind).Hmax;
                else
                    Bw = interp1(wls, Bws, frames(ind).winlen);
                    Bw = min(Bw/3,0.5*f0s_new(ind));
                    Hnew = floor(Bw/abs(df0s(ind)));
                    frames(ind).H = max(frames(ind).H,Hnew); % Force increment
                    frames(ind).H = max(frames(ind).H,opt.Kimin);
                    frames(ind).H = min(frames(ind).H,frames(ind).Hmax);
                end
            end

            if opt.debug>0; pb = progressbar(pb,ind); end
        end
        if opt.debug>0; pb = progressbar(pb,nt); end

        f0s = f0s_new;
        aks = aks_new;
        
        if numel(aks{1})==0; keyboard; end

        % Stopping conditions
        isfinished = [frames.H]'>=[frames.Hmax]';
        hasconverged = abs(df0s).*[frames.Hmax]' < opt.f0_prec*f0s;
        if adapt==1;    dSNR = SNRs(:,adapt);
        else            dSNR = SNRs(:,adapt)-SNRs(:,adapt-1); end
        isimproving = dSNR>opt.snr_prec;
        isdiverging = dSNR<-opt.snr_prec;

        disp('');
        disp(['        max|df0|=' num2str(max(abs(df0s))) ' min(RCN)=' num2str(min(log10([frames.RCN]))) ' max(dSNR)=' num2str(max(dSNR))]);

        disp(['        ' num2str(nt) ' frames: ' num2str(sum(isfinished)) ' finished, ' num2str(sum(hasconverged)) ' converged, ' num2str(sum(isimproving)) ' improving (' num2str(sum(isdiverging)) ' diverging)']);

        % Print current adaptation
        if opt.debug>1
            hold off;
            plot((0:length(wav)-1)/fs, wav, 'k');
            hold on;
            plot(T, log2(f0sin(:,2)), 'k');
            plot(T, log2(f0s), 'r');
            plot([T(1), T(end)], log2(opt.f0min)*[1 1], ':r');
            plot([T(1), T(end)], log2(opt.f0max)*[1 1], ':r');
            plot(T, [frames(ind).H]'./[frames(ind).Hmax]', 'k');
            plot(T, (abs(df0s).*[frames.Hmax]')./(opt.f0_prec*f0s)+0.05, 'b');
            plot(T, SNRs(:,adapt)./10, 'g');
            plot(T, dSNR, 'r');
%              plot(T, log10(rcns), ':r');

            stem(T(toups), -1*ones(size(T(toups))), 'or');

            ylim([-5 4]);

            legend({'wav','f0 in','f0 out','fullband','|df0|','SNR','dSNR','rcn'});

%              pause
            keyboard
        end

        % Next turn
        adapt = adapt + 1;

        SNRs(:,adapt) = SNRs(:,adapt-1); % because not all frames will be updated

        % Update only the frames which need to be updated
        % see which frame changed ...
        haschanged = isimproving | isdiverging;
        % ... and update all the frame depending on the changed frames + those which are not finished
        toups = find(~isfinished | sum(U(:,find(haschanged)),2)>0 |  ~hasconverged);
    end
    adapt = adapt - 1;

    f0sout = [T, f0s];

    % Drop components close to Nyquist which have zero amplitude
    for n=1:numel(aks)
        while aks{n}(end)==0 && numel(aks{n})>1
            aks{n} = aks{n}(1:end-1);
        end
    end

    % Fill the frames structures with f0s and sinusoidal parameters
    for ind=1:numel(frames)
        frames(ind).f0 = f0s(ind);
        Ks = 0:length(aks{ind})-1;
        frames(ind).sin = [f0s(ind)*Ks; abs(aks{ind}); angle(aks{ind}); Ks];
    end

    if opt.do_final_ahm_step
        optfinstep = sin_analysis();
        optfinstep.win_durnbper  = opt.win_durnbper;
        optfinstep.win_ls_marg= opt.win_ls_marg;
        optfinstep.win_fn     = opt.win_fn;
        optfinstep.use_f0time = true;
        optfinstep.harmonic   = true;
        optfinstep.fquasi     = false;
        optfinstep.fadapted   = true;
        optfinstep.use_ls     = true;
        frames = sin_analysis(wav, fs, f0sout, optfinstep);
    end

    if opt.debug>1
        disp('Ready to leave');
%          times = (0:length(wav)-1)/fs;
%          plot(times, wav, 'k');
%          hold on;
        keyboard
    end

return


