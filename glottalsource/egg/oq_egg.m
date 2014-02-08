% Open quotient measurement on EGG signal
%
% Octave compatible
% 
% Description
%  Measures the open quotient Oq on a EGG signal, knowing the sampling
%  frequency fs and the desired method (see below).
%  If freq_range is precised, freq_range = [f0min f0max],
%  otherwise freq_range = [80 1500].
%  If freq_range = f0, the fundamental frequency is not estimated in the function.
%
%  A minimum of two cycles (three gci) is required for the analysis.
%
% Inputs
%  s         : The input EGG signal.
%  fs        : [Hz] The sampling frequency of the EGG signal.
% 
%  method    : Select the method to compute Oq. It can be:
%              '35'      - 35% threshold method on EGG signal
%              '50'      - 50% threshold method on EGG signal
%              'howard'  - on DEGG and EGG with a threshold of 3/7
%
%  f0        : [1x2 Hz] The search limits of f0, in a row vector.
%              Default value is: f0=[80Hz, 1500Hz].
%              If a single value is used (e.g. f0=120), it replaces the
%              measured f0 and estimate Oq using this f0 value.
%
% Outputs
%  Oqeggmes  : Measured Open Quotient (Oq)
%  f0mes     : Measured Fundamental frequency
%
% Example
%  See the HOWTO_egg.m example file.
%
% Reference
% [1] N. Henrich, C. d'Alessandro, B. Doval, and M. Castellengo, "On the use
%     of the derivative of electroglottographic signals for characterization
%     of nonpathological phonation," Journal of the Acoustical Society of
%     America, vol.115, no.3, p.1321-1332, 2004.
%     Original publication of the code: http://voiceresearch.free.fr/egg
%
% Copyright (c) 2005 CNRS, Laboratoire d'Acoustique Musicale (Paris)
%                          and Institut de la Communication Parlée (Grenoble)
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
%  Nathalie Henrich <nathalie.henrich@gipsa-lab.fr> (May 20th, 2005)
%                                                   Revision : 2006-04-04
%

function [Oqeggmes,f0mes] = oq_egg(s, fs, f0, method)

if nargin <3 || isempty(f0)
    f0min = 80;
    f0max = 1500;
    f0 = [NaN;NaN];
elseif length(f0) == 2
    f0min = f0(1);
    f0max = f0(2);
end

% data length
L = length(s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNDAMENTAL FREQUENCY ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(f0) == 1
    f0mes = f0;
else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculuation of normalized (and biased) autocorrelation function
    % (non circular autocorrelation)
    %
    % equivalent method : (circular autocorrelation)
    %
    % autoc = real(ifft(abs(fft(s,2*L)).^2)) / L;
    % autoc = autoc(1:L);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    autoc = xcorr(s,'coeff');
    autoc = autoc(L:2*L-1);

    % peak detection corresponding to periods
    % threshold of 0.5 (voiced / unvoiced criterion)
    [maxa,indmaxa]=PicDetect(autoc,0.5);

    if isempty(maxa)
        f0mes=0;
        Oqeggmes=0;
    elseif autoc(max(1,indmaxa(1)-4))>=autoc(indmaxa(1))
        f0mes=0;
        Oqeggmes=0;
    else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fundamental frequency estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % rough estimation
        f0approx = fs/(indmaxa(1)-1);

        % precision of 0.5 Hz
        deltat = 0.5 / f0approx^2;

        % cubic interpolation between 5 points
        % centered around the detected peak
        indsel = [max(1,indmaxa(1)-2):indmaxa(1) indmaxa(1)+1 indmaxa(1)+2];
        tsel = (indsel-1)/fs;
        autocsel = autoc(indsel);

        % resampling
        ti = min(tsel):deltat:max(tsel);

        % interpolation
        autoci = interp1(tsel,autocsel,ti,'spline');
        % maximum
        [ymax1,indmax1] = max(autoci);

        % f0 estimation de f0 with a precision of 0.5 Hz
        t1 = ti(indmax1);
        f0mes = 1/t1;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % frequency range check
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ((f0mes>=f0max) || (f0mes<=f0min))
            f0mes=0;
            Oqeggmes=0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN QUOTIENT ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if f0mes~=0

    % zero-mean signal
    s = s-mean(s);
    % fundamental frequency in samples
    n0 = fix(1/f0mes*fs);
    % min and max detection
    [val_max,ind_max] = PicDetect(s,n0);
    [val_min,ind_min] = PicDetect(-s,n0);

    if length(ind_max)<=2 || length(ind_min)<=2
        Oqeggmes = 0;
    else

        % start on a minimum (open phase)
        while ind_max(1) <= ind_min(1)
            ind_max(1) = [];
            val_max(1) = [];
        end

        switch method
            case {'35';'50'}
                nlevel = str2num(method)/100;
            case 'howard'
                nlevel = 3/7;
            otherwise
                error(['The method ''' method ''' is not handled by oq_egg ...'])
        end

        % open quotient measurement between two min
        for k=1:min(length(ind_min)-1,length(ind_max))
            if ind_max(k) <= ind_min(k)
                warning('ind_max <= ind_min')
                Oqegg(k) = 0;
            else

                s_sel = s(ind_min(k):ind_min(k+1));
                threshold =  min(s_sel) + nlevel * (max(s_sel)-min(s_sel));

                switch method
                    case {'35';'50'}
                        indf = min(find(s_sel>=threshold));

                        % refining the measurement
                        % y = a*t+b with b = y2, a = (y2-y1)/(t2-t1)
                        % tx = t2 + (yx-y2)/a
                        if indf == 1
                            tf = 1/fs;
                        else
                            t2 = indf/fs; y2 = s_sel(indf);
                            t1 = (indf-1)/fs; y1 = s_sel(indf-1);
                            a = (y2-y1)/(t2-t1);
                            tf = t2 + (threshold-y2)/a;
                        end

                    case 'howard'
                        [val,indf] = max(diff(s_sel));
                        indf = indf+1; % décalage de 1 du à la dérivation
                        tf = indf/fs;

                    otherwise
                        error(['The method ''' method ''' is not handled by oq_egg ...'])
                end


                indo = max(find(s_sel>=threshold));
                % refining the measurement
                % y = a*t+b with b = y2, a = (y2-y1)/(t2-t1)
                % tx = t2 + (yx-y2)/a
                if indo == length(s_sel)
                    to = length(s_sel)/fs;
                else

                    t2 = indo/fs; y2 = s_sel(indo);
                    t1 = (indo+1)/fs; y1 = s_sel(indo+1);
                    a = (y2-y1)/(t2-t1);
                    to = t2 + (threshold-y2)/a;
                end

                % Oqegg = open phase / fundamental period = 1 - closed phase*f0
                Oqegg(k) = 1-(to-tf)*f0mes;
                t_f(k) = (ind_min(k)+indf-1)/fs;
                %                 plot([1:length(s_sel)],s_sel,indf,s_sel(indf),'*k',indo,s_sel(indo),'or')
                %                 title(num2str(Oqegg(k)))
                %                 drawnow
                %                 waitforbuttonpress
            end
        end

        ind = find(Oqegg>1|Oqegg<=0);
%         if length(ind)>=length(Oqegg)/4
%             f0mes = 0;
%             Oqeggmes = 0;
%             return;
%         end
        if ~isempty(ind)
            Oqegg(ind) = [];
        end

        if ~isempty(Oqegg)
            Oqeggmes = mean(Oqegg);
        else
            Oqeggmes = 0;
        end
    end
end


