% SEDREAMS is a method for Glottal Closure Instant (GCI) determination.
%
% Octave compatible
%
% Description
%  The Speech Event Detection based on the Residual Excitation And
%  Mean-based Signal (SEDREAMS) is described in [1] and [2]. It acts in two
%  successive steps. First short intervals where GCIs are expected to occur
%  are derived from the mean-based signal. This ensures good reliability
%  performance. In the second step, the precise location of each GCI is
%  refined by inspecting the residual excitation signal. This leads to good
%  accuracy performance.
%
% Inputs
%  wave             : [samples] [Nx1] input signal (speech signal)
%  fs               : [Hz]      [1x1] sampling frequency
%  f0mean           : [Hz] [1x1] an estimate of the averaged F0 used by
%                     the speaker. This can be estimated using the "pitch_srh"
%                     function included in this toolbox (or using any other pitch
%                     tracker).
%  polarity         : [1x1] the polarity of the speech signal (-1 or +1).
%                     If unknown, you can use the "polarity_reskew" function
%                     available in this toolbox.
%
% Outputs
%  gci              : [s] vector containing the location of the estimated GCIs.
%                     Note that GCIs are estimated even during unvoiced segments.
%                     Proper voicing boundaries should then be used as a mask
%                     (you can use the "pitch_srh" function included in this 
%                     toolbox for this).
%  MeanBasedSignal  : [samples] [Nx1] mean-based signal as used by the
%                     SEDREAMS algorithm
%  res              : [samples] [Nx1] the residual excitation signal
%
% Example
%  And see the HOWTO_glottalsource.m example file.
%
% References
% [1] T.Drugman, T.Dutoit, Glottal Closure and Opening Instant Detection
%     from Speech Signals, Interspeech09, Brighton, U.K, 2009
%
% [2] T.Drugman, M.Thomas, J.Gudnason, P.Naylor, T.Dutoit, Detection of
%     Glottal Closure Instants from Speech Signals: a Quantitative Review,
%     IEEE Transactions on Audio, Speech and Language Processing, vol. 20,
%     Issue 3, pp. 994-1006, March 2012

%  Publications available at the following link:
%  http://tcts.fpms.ac.be/~drugman/files/TASLP-GCI.pdf
%  http://tcts.fpms.ac.be/~drugman/files/IS09-GCIdetection.pdf
%
% Copyright (c) 2011 University of Mons, FNRS
%
% License
%  This code is a part of the GLOAT toolbox with the following
%  licence:
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Thomas Drugman thomas.drugman@umons.ac.be


function [gci, MeanBasedSignal, res] = gci_sedreams(wave, fs, f0mean, ...
                                                    polarity, opt)
    

    if nargin<5
        % Options
        opt.use_maxlpresidual = false; % Use the traditional LP
                                       % residual by default, as in
                                       % [1,2]
    end
    if nargin==0; gci=opt; return; end

    wave=polarity*wave;

    if ~opt.use_maxlpresidual
        res = lpcresidual(wave,round(25/1000*fs),round(5/1000*fs), ...
                          round(fs/1000)+2);
        
    else
        res = maxlpresidual(wave,fs,round(fs/1000)+2);
    end
    res(isnan(res))=0;

    % Calculation of the mean-based signal
    T0mean=round(fs/f0mean);
    halfL=round((1.7*T0mean)/2);
    Blackwin=blackman(2*halfL+1);

    % filter wave with blackwin and take mean
    MeanBasedSignal = filter(Blackwin,numel(Blackwin),wave)';
    % shift output since MATLAB's filter 'center' is the first element
    MeanBasedSignal(halfL+1:end-halfL) = MeanBasedSignal(1+2*halfL:end);
    % and pad begin and end with zeros
    MeanBasedSignal([1:halfL end-halfL+1:end]) = 0;
    
    % Remove the low-frequency contents of the mean-based signal
    Ws = 30/(fs/2);
    Wp = 50/(fs/2);
    Rp = 3; Rs = 60;
    [n,Wp] = ellipord(Wp,Ws,Rp,Rs);
    [b,a] = ellip(real(n),Rp,Rs,Wp,'high');

    MeanBasedSignal=filtfilt(b,a,MeanBasedSignal);
    MeanBasedSignal=MeanBasedSignal/max(abs(MeanBasedSignal));

    % Detect the minima and maxima of the mean-based signal
    [~,PotMaxis]=findpeaks(MeanBasedSignal);
    [~,PotMinis]=findpeaks(-MeanBasedSignal);

    while PotMaxis(1)<PotMinis(1)
        PotMaxis(1)=[];
    end

    while PotMinis(end)>PotMaxis(end)
        PotMinis(end)=[];
    end

    Minis=PotMinis;
    Maxis=PotMaxis;


    % Determine the median position of GCIs within the cycle
    res=res/max(abs(res));
    Posis=find(res>0.4);
    RelPosis=zeros(1,length(Posis));
    for k=1:length(Posis)

        Dists=abs(Minis-Posis(k));
        [~,pos]=min(Dists);
        interv=Maxis(pos)-Minis(pos);

        RelPosis(k)=(Posis(k)-Minis(pos))/interv;
    end

    if isempty(RelPosis)==0
        RatioGCI=median(RelPosis);
    else RatioGCI=0;
    end


    % Detect GCIs from the residual signal using the presence intervals derived
    % from the mean-based signal
    gci=zeros(1,length(Minis));

    Ind=1;
    for k=1:length(Minis)
        interv=Maxis(k)-Minis(k);    
        alpha=RatioGCI-0.35;
        start=Minis(k)+round(alpha*interv);    
        alpha=RatioGCI+0.35;
        stop=Minis(k)+round(alpha*interv);

        if start<1
            start=1;
        elseif start > length(res)
            break
        end
        if stop>length(res)
            stop=length(res);
        end
        if (stop > 1)
            vec=res(start:stop);
            [~,posi]=max(vec);
            gci(Ind)=start+posi(1)-1;
            Ind=Ind+1;
        end
    end


    gci = (gci-1)/fs;

end

