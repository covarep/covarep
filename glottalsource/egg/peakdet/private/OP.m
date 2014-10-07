function [Oq,Oqval,DEOPA,goodperiods] = OP(simppeak,SIG,dSIG,valid,doublepeaks,COEF,resampC,Ts,method)
% OP: Opening Peaks detection.
% Input files:
% - <simppeak>: contains a list of closing-peak times in ms, and amplitudes
% (<simppeak> is a 2-column matrix)
% - <dSIG>: Derivative of EGG signal
% - <valid> : vector containing information on peaks that have several summits
% - <doublepeaks>: vector containing information on peaks that cluster together
% - <COEF>: contains, among other coefficients, the smoothing step chosen by the
% user. It is not used in analyzing the positive part of the signal, only in the
% negative part, for the calculation of Oq.
% - <resampC>: coefficient for reinterpolation
% - <Ts>: the inverse of the sampling rate
% - <method>: indication on the method which is to be used in the handling of 
% double / multiple peaks: selecting either the first or the highest. 
% (The fact that there is not one single peak will also be passed on
% to the following functions, so this choice is not of much consequence.)
%
% Output files: 
% <Oq>: vector containing open quotient values detected using the maxima
% <Oqval>: open quotient values computed by peak detection (single peaks only)
% <DEOPA>: vector containing DEOPA values
% <SM>, for Success of Method: variable indicating the number of Oq values that
% could be calculated using peak detection.
% <goodperiods>: beginning and end of cycles for which closing peaks are unique.
% <validGOI>: ratio of highest to second highest peak.
%
% No use of autocorrelation or comparison across cycles, otherwise the method 
% would not give results up to the end of voicing: border phenomena.
% Passages where voice quality changes quickly are not found only at the offset
% of voicing: they can be word-internal, as in the case of consonantal glottal
% stops, ejectives, implosives, or glottalized tones, e.g. in Vietnamese tone
% C2, which has medial glottal constriction.
%
% Method: detecting maxima; checking the plausibility of results; if results are
% implausible, detecting maxima anew with increased smoothing step

% Setting the output vector Oqval at []; otherwise, in case no single period has 
% single closing and opening peaks, this output argument is not assigned during
% call to OP.
Oqval = [];

% creating a matrix containing the beginning and end of periods for which Oq
% will be calculated. Values are converted to indices.
goodperiods = [];

if ~isempty(simppeak)
    if length(simppeak(:,1)) > 1
        for i = 1: length(simppeak(:,1)) - 1
            goodperiods(i,1) = round(simppeak(i,1) * (1 / Ts));
            goodperiods(i,2) = round(simppeak(i + 1,1) * (1 / Ts));
        end

        % Previous calculation method: eliminating periods for which closing peaks are
        % not unique:
        % %%%%%%%%%%%%%%%%%%%%%% Determining periods for which closings are unique
        % % Open quotient measurements are attempted only when closings are unique. For
        % % each glottal cycle, there are two condition: 
        % % 1) that the values associated to the peak in <valid> be
        % % above 0.6, i.e. within the peak there is no other summit that is higher than
        % % 60% of the highest summit
        % % 2) that none of the neighbouring peaks be higher than 60% of the candidate peak.
        % 
        % goodperiods = []; 
        % for i = 1:length(valid) - 1
        %     if (doublepeaks(i) == 0) & (valid(i) < 0.6) & (doublepeaks(i + 1) == 0) & (valid(i) < 0.6)
        %         goodperiodsnb = goodperiodsnb + 1;
        %         goodperiods(goodperiodsnb,1) = Tgci(i,1);
        %         goodperiods(goodperiodsnb,2) = Tgci(i + 1,1);
        %     end
        % end


        %%%%%%%%%%%%%%%%%%%%%% First measurement of openings: done by detecting
        %%%%%%%%%%%%%%%%%%%%%% local minima, without any condition on the values obtained.
        texis = exist('goodperiods');
        if texis == 1
            if length(goodperiods(:,1)) > 0
                for i = 1:length(goodperiods(:,1))
                    % A rounding of the indices is necessary: as there has been
                    % reinterpolation, the values are not integers.
                    [DEOPA(i),GOI(i)] = min(    dSIG( round(goodperiods(i,1)):round(goodperiods(i,2) ))   );
                    % correction: adding the index corresponding to the first point of the
                    % excerpt
                    GOI(i) = GOI(i) + goodperiods(i,1);
                    %% Calculation of Oq, in %
                    Oq(i) = 100 * ((  goodperiods(i,2) - GOI(i)  ) / (  goodperiods(i,2) - goodperiods(i,1)  ));
                end
        %     else
        %         DEOPA = 0;
        %         GOI = 0;
        %         Oq = 0;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%% Measurement by peak detection and barycentre calculation
        texis = exist('goodperiods');
        if texis == 1
            if length(goodperiods(:,1)) > 0
                for i = 1:length(goodperiods(:,1))
                    % if you want to see information displayed as analysis is being
                    % conducted, uncomment the two lines below:
                    %% disp(['Treating period number ',num2str(i)])
                    %% disp(goodperiods(i,:))
                    % placing the indices of beginning and end of the relevant portion of
                    % the signal in variables
                    % A rounding of the indices is necessary: as there has been
                    % reinterpolation, the values are not integers. In the present state of
                    % the programme, this rounding is performed earlier, when the
                    % <goodperiods> matrix is built up.
                    BE = goodperiods(i,1);           % BEginning of period
                    EN = goodperiods(i,2);           % ENd of period
                    % placing the extract in a vector, inverting it because CRO operates on
                    % a positive signal. An extra precaution must be taken: excluding the
                    % first and last eight points of the signal (by setting them at zero):
                    % it is physically impossible that a clear opening will take place 
                    % at these points, so close to the closing; technically, a small part of
                    % the signal immediately before the closing peak that marks the
                    % beginning of the period (res. after the closing peak that marks its
                    % end) may be included and result in program crash.
                    dSIGextr = - dSIG(BE:EN);
                    L = length(dSIGextr);
                    if L > 16
                        dSIGextr(1:8) = 0;
                        dSIGextr(L-7:L) = 0;
                    end
                    % Also creating an extract for the EGG signal, which is used later on,
                    % by AMPOS
                    SIGextr = - SIG(BE:EN);

                    %%%%%%%%%%%%%%%%%% Calling the function CRO to detect crossings. 
                    % The <method> toggle (3rd value in <COEF>) must be set at 0, and the
                    % amplitude threshold (4th value in <COEF>) passed on to the function.
                    COEFOPEN = COEF;
                    COEFOPEN(3) = 0;

                    % setting the threshold for peak detection (empirically). Choice in
                    % Henrich 2001:123 : setting it at 70% of the highest point.
                    COEFOPEN(4) = 0.70 * max(dSIGextr);

                    % finding out the crossings with the DEGG signal. In the best cases, it is
                    % expected that there will be only two crossings, one up and one down.
                    % But this is not frequent: the signal is often dented.
                    [rimsOPEN] = CRO(dSIGextr,COEFOPEN);

                    if isempty(rimsOPEN)
                        % If there is no single detected threshold crossing: the user should
                        % changed the maximum F0 value and ask for the results to be
                        % calculated anew. The user is advised to do so inside the <rim>
                        % function.
                        Oqval(i) = 0;
                        DEOPA(i) = 0;
                    else
                        % finding out the amplitude and position of the peaks. No
                        % reinterpolation of the signal is useful, as the signal-to-noise ratio in
                        % the negative part of the DEGG signal is poor. Quite the opposite: the
                        % signal is smoothed (over 3 samples, i.e. 0.07 ms).

                        % The threshold for exclusion of small peaks is set at 0.7 in the case of 
                        % opening peaks, following Henrich 2001:123. This has in fact already
                        % been applied before, in detection of peak rims by RIM; but the
                        % function PEAKSHAPE called by AMPOS needs this input.
                        propthresh = 0.7;
                       % figure(1)
                       % clf
                       % plot(SIGextr)
                        [TGOI,TGOIFo,validGOI,validtimeGOI] = AMPOS(SIGextr,rimsOPEN,1,Ts,method,propthresh);

                        % detection of double / multiple closing peaks ("peak clusters"). The
                        % threshold passed to this function is: (the inverse of) one-fifth of
                        % the period, under the assumption that peaks separated by more than
                        % one-fifth of the period should not be considered as a bundle. In these
                        % cases, no Oq is calculated.
                        maxF = 1 / ( (goodperiods(i,2) - goodperiods(i,1)) / 5);
                        doublepeaks = detectmult(TGOI,maxF);

                        % choice of opening peak, leaving only one, following the method
                        % specified by the user. If there are 2 peaks that are too far apart for
                        % averaging to make sense, there will be 2 lines in simppeak: two "real
                        % peaks" according to the criteria used in the SIMP function.
                        [simppeak,BUND] = simp(TGOI,doublepeaks,method);

                        % calculating the open quotient when possible. The condition is: if
                        % there is one single value in <simppeak>, the matrix where detected
                        % peaks are stored. If there are two (or more) values, this indicates
                        % that there were two (or more) peaks that were too far apart for
                        % averaging to make sense, and it seems best not to give an Fo value;
                        % the Oq value is then set at zero. This is also what happens if it 
                        % remains unaffected: it comes out as a zero inside the vector; but if
                        % the last values are zero, then the length of the vector is less than
                        % the size of the Fo vector, resulting in assignment problems. So the
                        % value is set at zero in this programme.
                        % inside the vector).
                        if length(simppeak(:,1)) == 1
                            % transforming the time of peak to an index
                            indexsimppeak = simppeak(1,1) * (1 / Ts);
                            Oqval(i) = 100 * (       ((goodperiods(i,2) - goodperiods(i,1)) - indexsimppeak ) / ...
                                  (goodperiods(i,2) - goodperiods(i,1))       );
                        else
                              Oqval(i) = 0;
                        end
                    % end of "if" loop: condition on existence of rimsOPEN values.    
                    end
                end
            end
        end
    else
        goodperiods = []; Oq = []; Oqval = []; DEOPA = [];
    end
else
    goodperiods = []; Oq = []; Oqval = []; DEOPA = [];
end
        %%%%%%%%%%% The method for detecting opening peaks is fairly similar to
        %%%%%%%%%%% that for detecting closing peaks. Concerning the detection
        %%%%%%%%%%% of multiple peaks, however, in the case of closing peaks
        %%%%%%%%%%% the threshold for detection was set rather low, and the
        %%%%%%%%%%% peaks were later sorted according to size and closeness to
        %%%%%%%%%%% one another. In the case of opening peaks, the period is
        %%%%%%%%%%% known already, so the threshold is fixed straight away,
        %%%%%%%%%%% on the basis of the highest point in this part of the
        %%%%%%%%%%% signal: if the second highest peak is 70% as high as the
        %%%%%%%%%%% highest peak, it is decided that the peak is not unique
        %%%%%%%%%%% (following Henrich 2001:123).

        %%%%%%%%%%% It is an issue whether this threshold must apply
        %%%%%%%%%%% within the peak region as well as in the comparison of
        %%%%%%%%%%% neighbouring peaks. The fine detail of the peak is often
        %%%%%%%%%%% slightly dented (see figures); this will result in
        %%%%%%%%%%% "secondary peaks" being detected very close to the maximum
        %%%%%%%%%%% (about 0.1 ms to its left and right). A time condition could
        %%%%%%%%%%% be added: if the secondary peak is very close
        %%%%%%%%%%% (threshold empirically set at a time value, e.g. 0.5 ms, or
        %%%%%%%%%%% as a proportion of the period), it would not count as
        %%%%%%%%%%% a doubled peak. In the present programme, the solution
        %%%%%%%%%%% chosen instead is to use smoothing of the signal.
        
%         % In an earlier version, the algorithm was: 
%         % 1) retrieve highest peak in TGOI; 2) check that it is unique
%         % internally; 3) check that the neighbouring peaks are less than 70% of its
%         % height. GOIC stands for Glottis-Opening-Instant Candidate.
%         % First step:
%         [GOICvalue,GOICindex] = max (TGOI(:,2));
%         % Second step: finding out whether the peak is unique internally. For
%         % this: look up the corresponding slot in vector <validtimeGOI>.
%         if validtimeGOI(GOICindex) < threshtime
%             % retrieving the index of maximum value in TGOI relative to the 
%             % beginning of the <dSIGextr> extract: it is the value
%             % in the first column of <TGOI> at the line indicated in
%             % <GOICindex>
%             GOICindex = TGOI(GOICindex,1);
%             % changing the relative index <GOICindex> into an index relative 
%             % to the <SIG> signal as a
%             % whole, for computation of the open quotient
%             GOICindex = GOICindex + goodperiods(i,1);
%             % calculation of open quotient
%             Oqval(i) = 100 * (       (goodperiods(i,2) - GOICindex ) / ...
%                 (goodperiods(i,2) - goodperiods(i,1))       );
%         end
            
            %%%%%%%%%%% Note:             
            %% In the first version of the programme, the threshold for detection 
            %% of peak regions was first set at 40% of the signal, and the
            %% condition on peak amplitude was added later, by recovering the 
            % second highest peak in TGOI and
            % comparing its amplitude to that of the highest peak. 
            %% In the present version of the function, only the peaks above 70% of the maximum
            %% within the dSIG extract <dSIGextr> are detected, so that
            %% this part of the function is not necessary anymore.
            %% The calculations were as follows:
            
%             COEFOPEN(4) = 0.7 * max(dSIGextr);
%             % Setting value for TGOI(GOICindex,2) at zero and re-calculating the
%             % maximum: 
%             TGOI(GOICindex,2) = 0;
%             [GOICvalue2,GOICindex2] = max (TGOI(:,2));
%             % condition on whether the proportion of the two peaks is below
%             % threshold
%             if GOICvalue2 / GOICvalue < threshGOI
%                 % retrieving the index of maximum value in TGOI: it is the value
%                 % in the first column of <TGOI> at the line indicated in
%                 % <GOICindex>
%                 GOICindex = TGOI(GOICindex,1);
%                 % changing the relative index <GOICindex> into an index relative 
%                 % to the <SIG> signal as a
%                 % whole, for computation of the open quotient
%                 GOICindex = GOICindex + goodperiods(i,1);
%                 % calculation of open quotient
%                 Oqval(i) = 100 * (       (goodperiods(i,2) - GOICindex ) / ...
%                     (goodperiods(i,2) - goodperiods(i,1))       );
%             end