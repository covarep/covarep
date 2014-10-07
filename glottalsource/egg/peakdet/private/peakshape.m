function [I,PA,valid,validtime,numberofpeaks] = peakshape(dREI,method,resampC,Ts,propthresh)
%%%%%%%%%%%%%%%% PEAKSHAPE: detects whether there is more than one peak
%%%%%%%%%%%%%%%%  in the portion of signal given as input, and if so,
%%%%%%%%%%%%%%%%  indicates the size of the second peak in proportion to the
%%%%%%%%%%%%%%%%  first, e.g. 0.5 (50%), and the time distance between the two.
% Rationale: It may happen that the part of the peak that is above the
% threshold is not a unique peak; if it is dented, this needs
% to be reported to the user.
% Input: 
% - portion of signal (derivative of re-interpolated signal)
% - method for handling double / multiple peaks: selecting either the first or
% the highest. (The fact that there is not one single peak will also be passed on
% to the following functions, so this choice is not of much consequence.)
% - <propthresh> : threshold for exclusion of small peaks; used only for
% glottis closures: for openings, the period is known already, and the threshold
% set at a certain height of the maximum of the signal within the period, earlier on
% in the calculations.
% Output: 
% 1) <I>, for Instant: the time of the point detected as peak
% 2) <PA>, for Peak Amplitude: the value of the signal at the point detected as peak
% (i.e. DECPA in the case of positive peaks, DEOPA for opening peaks).
% 3) <valid>, set at:
% -> if the peak is unique: 0
% -> otherwise: the ratio between the amplitude of the second and first
% peaks
% -> if the user specifies that double peaks must be excluded: <valid> is set at
% 2 in case of double peaks.
% 4) <validtime>: the distance in time between highest and second highest peaks
% 5) <numberofpeaks>: very simply, the number of peaks (or: "maxima") within the
% peak region.
%
% This function expects the peak to be positive, i.e. looks for positive
% values. If a negative peak is treated, it must have its sign inverted,
% i.e.: valid = peakshape(- dREI)

% setting the variables <valid> and <validtime> at zero (default value) to begin with
valid = 0;
validtime = 0;

% recovering the value and index (the latter being a relative indicator: 
% indication of position within <peakvect> vector) of the highest peak
[PA,I] = max(dREI);

% computing the second derivative of the signal
for i = 1: length(dREI) - 1
    ddREI(i) = dREI(i + 1) - dREI(i);
end

% retrieving the sign changes FROM POSITIVE TO NEGATIVE in the second
% derivative, i.e. the number of peak tips. If the portion of the signal contains
% one single peak, there should be only one such sign change to the negative.
numberofpeaks = 0; 
peakvect = [];
for i = 1: length(ddREI) - 1
    if (ddREI(i)>0) & (ddREI(i + 1)<0)
        % An additional condition is placed: one on peak amplitude. In the case
        % of opening peaks, each period is treated individually, and the
        % threshold is placed according to the local maximum; so no further
        % check needs to be placed at the present point; it is redundant. But
        % in the case of
        % closing peaks, the threshold is set for the whole sound extract
        % (typically, one syllable), so a check needs to be placed on the peaks
        % at the present stage. 
        if dREI(i + 1) / PA > propthresh
            numberofpeaks = numberofpeaks + 1;
            peakvect(numberofpeaks) = i + 1;
        end
    end
end

%%% if extra sign changes were detected: finding out which is the peak that is
%%% to be selected for the calculation of Fo and Oq
if numberofpeaks > 1
    %%%%%% calculating the amplitude ratio, and the distance, between 
    %%%%%% highest peak and each of the other peaks
    %%% The value and index of the highest summit are already stored in [PA,I].
    %%% But the index is relative to the sound extract as a whole, not to the
    %%% vector <peakvect>. The value is therefore calculated anew.
    [PA,Irel] = max(dREI(peakvect));

	% suppressing this maximum from <peakvect>
    TE = [];
	if Irel - 1 > 0
        for i = 1: (Irel - 1)
            % TE: for TEmporary variable
            TE (i) = peakvect(i);
        end
	end
	if Irel + 1 <= length(peakvect)
        for i = (Irel + 1) : length(peakvect)
            TE(i - 1) = peakvect(i);
        end
	end
	peakvect = [];
    peakvect = TE;
    
    % for calculation of Fo and Oq: yielding a single value, following the method
	% chosen by the user
    if method == 0
        %%%%% method "0": selecting the highest peak
        % This is what has already been done at the beginning of this function: 
        % recovering the value and index (the latter being a relative indicator: 
        % indication of position within <peakvect> vector) of the highest peak:
%         [PA,Irel] = max(dREI(peakvect));
        % converting Irel (index within <peakvect>) to an index within <dREI>
%         I = peakvect(Irel);
    elseif method == 1
        % method "1": selecting the first peak
        Irel = 1;
        I = peakvect(Irel);
        PA = dREI(I);
    elseif method == 2
        % method "2": selecting the last peak
        Irel = length(peakvect);
        I = peakvect(Irel);
        PA = dREI(I);
    elseif method == 3
        %% method "3": using barycentre method. Principle: if the first peak is
        %% at 0 in time and there is another peak at 100 in time, and the second is n % as
        %% high as the first, a time value in-between is computed: pondered mean of
        %% time values A and B. The present programme can deal with more than 2
        %% peaks, e.g. if there are 3 "peaks" at closure, as is sometimes
        %% observed. Variables: ICO for Index of COefficiented peak; PACO, for
        %% Peak Amplitude of COefficiented peak.
        
        % loop for peaks
        % for coefficiented mean: using a WEIGHT variable, set at 1 to
        % begin with, and then augmented as more peaks are added to the mean,
        % one by one.
        WEIGHT = 1;
        for i = 1:length(peakvect)
            % time distance between the peak being examined and the highest peak
            % (in index, relative to resampled signal)
            dist(i) = peakvect(i) - I;
            % This variable is
            % not an absolute value (abs(...)): it is negative if the peak being
            % examined comes before the main peak, positive if it comes after.
            
            % amplitude ratio between the peak being examined and the highest peak
            % (in %)
            ampratio(i) = dREI(peakvect(i)) / PA;
            
            % halving <dist(i)>, and adding a coefficient, to obtain the
            % distance to be subtracted from I
            dist(i) = (dist(i) / 2) * ampratio(i);
            % dividing by WEIGHT, to take into account all the peaks that
            % are averaged
            dist(i) = dist(i) / WEIGHT;
            % increasing WEIGHT
            WEIGHT = WEIGHT + ampratio(i);
            %%%%% correcting the peak position according to the coefficiented distance
            I = I + dist(i);
            % The amplitude value (called DECPA for closing peaks, DEOPA for opening peaks) 
            % is unchanged: it remains set at the peak amplitude of the highest maximum.
        % end of the loop for method 3    
        end

	%         % In previous versions, only the coordinates of the second highest peak (time and amplitude)
	%         % were retrieved, and stored in variables; either by a condition within the loop:       
	%         if dREI(peakvect(i)) > valuesecondpeak
	%             indexsecondpeak = peakvect(i);
	%             valuesecondpeak = dREI(indexsecondpeak);
	%         end
	%         % or by retrieving the maximum of <dist>, i.e. the distance (in seconds)
	%         % between the highest peak and the peak that is furthest from it in
	%         % time. 
	%         [validtime,POSITIO] = max(dist);
	%         % retrieving corresponding amplitude ratio
	%         valid = ampratio(POSITIO);
	
	
	        % Remember that
            % the coordinates of the second highest peak were recovered above, by
            % the following lines: 
	%         if (dist (i) > 5e-4) & ( dREI(peakvect(i)) > valuesecondpeak)
	%             indexsecondpeak = peakvect(i);
	%             valuesecondpeak = dREI(peakvect(i));
	%         end

    elseif method == 4
        % method "4": excluding all double peaks. The <valid> variable is set at
        % two, and the peak will be excluded from calculations.
        valid = 2;
    end

else
% case in which there is only one peak
    [PA,I] = max(dREI);
% end of the <if> loop
end