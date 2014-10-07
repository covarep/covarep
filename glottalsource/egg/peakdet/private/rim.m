function rims = rim(dSIG,thresh)
% Input: smoothed DEGG signal, and threshold for peak detection
% Output: matrix containing index of beginning and end of peaks: beginnings in
% left column, ends in right column.

% number of up-crossings (set at zero to begin with)
nbl = 0;
% number of right-crossings (set at zero to begin with)
nbr = 0;

% Loop: must start at point 2, and stop 1 point before the end, as the condition on the
% recognition of crossings starts 1 point before the i-th point and goes 1 point beyond.
for i = 2 : length(dSIG) - 1
    % The original programme detected any crossing of the threshold by the
    % following lines:
%     if and (dSIG(i) > thresh, dSIG(i + 1) <= thresh)
%        nbr = nbr + 1;
%        rims(nbr,1) = i ; 
%     end
%     if and (dSIG(i) < thresh, dSIG(i + 1) >= thresh)
%         nbl = nbl + 1;
%         rims(nbl,2) = i ; 
%     end
    % This method resulted in erroneous detections in the frequent cases where
    % there was a dent in the signal precisely at the threshold value. 
    % One solution envisaged consisted in adding a condition on the six
    % surrounding points (three to the left and right), to ensure that the point
    % in question is actually a crossing, not simply a dent in the signal. The
    % commands were as follow:
%     % Calculation of truth values for 3 points in succession, yielding a vector of three truth values
%     DECREASE_TRUTH = (dSIG(i - 2:i) >= thresh) & (dSIG(i + 1:i + 3) < thresh);
%     % Condition: if all three propositions are true (i.e. if their mean is 1),
%     % this is a down-crossing.
%     if mean(DECREASE_TRUTH) == 1
%         nbr = nbr + 1;
%         rims(nbr,2) = i ; 
%     end
%     % Same procedure for up-crossings.
%     INCREASE_TRUTH = (dSIG(i - 2:i) <= thresh) & (dSIG(i + 1:i + 3) > thresh);
%     if mean(INCREASE_TRUTH) == 1
%         nbl = nbl + 1;
%         rims(nbl,1) = i ; 
%     end
    % This condition, however, is too strong: at the point of crossing, there
    % may be dents in the signal, resulting in the crossing not being detected,
    % e.g. if there are "double crossings", so that no crossing was detected for
    % large stretches of the signal. The solution chosen is therefore to detect
    % all crossings, as in the first solution, only excluding "false crossings":
    % points where the limit is reached but not actually crossed, i.e. sequences
    % such as: 
    % threshold+3      threshold       threshold+1
    % or: 
    % threshold-3      threshold       threshold-1 :
    if (dSIG(i - 1) <= thresh) & (dSIG(i) <= thresh) & (dSIG(i + 1) > thresh)
        % (detection of an up-crossing)
        nbl = nbl + 1;
        rims(nbl,1) = i ; 
    end
    if (dSIG(i - 1) > thresh) & (dSIG(i) <= thresh) & (dSIG(i + 1) <= thresh)
       % (detection of a down-crossing)
       nbr = nbr + 1;
       rims(nbr,2) = i; 
    end
    % In case of successive values equal to the threshold: if they result from 
    % a decrease in the signal, a crossing is detected, both up and down; if
    % they result from an increase, no crossing is detected.
    % For example, if the threshold is 10: the
    % sequence 14 10 10 12 will be analyzed as: crossing-down at index 2,
    % crossing-up at index 3. The sequence 9 10 10 8 will be analyzed as: no
    % crossing at all (due to the > condition). This is of no consequence in the
    % results.
end

%%% Condition for cases where the first down-crossing precedes the first up-crossing (i.e. when the 
%%% excerpt cuts across a peak)
% Checking that the variable is not nonexistent (i.e. a value has been assigned to it)
EMPT = exist('rims');
if EMPT == 1
	if rims(1,1) > rims(1,2)
        for i = 2:length(rims(:,2))
            % suppression of final line
            rims(i - 1,2) = rims(i,2);
        end
        trimmedrims = rims(1:length(rims) - 1,:);
        rims = [];
        rims = trimmedrims;
	end

    % Similar condition for an up-crossing at the end of the extract.
    % In this case, a zero will have been added at the bottom right-hand corner of the matrix.
	if rims(length(rims(:,1)),2) == 0
        trimmedrims = rims(1:length(rims) - 1,:);
        rims = [];
        rims = trimmedrims;
    end
    
    % Cases where a down-crossing has been detected but no
    % up-crossing result in <rim> being equal to, e.g., [0 9].
    % In this case, <rim(1,1)> has to be modified, otherwise it cannot serve as
    % an index into a vector or matrix.
    if rims(1,1) == 0
        rims(1,1) = 1;
    end
else
    % case in which no single crossing was detected: an assignment must be made
    % to make it an empty matrix, and not a nonexistent matrix.
    disp(' ')
    disp(' ')
    error('No single peak has been detected in relevant portion of signal. See figure.')
    figure
    plot(dSIG)
    disp(' ')
    disp('Please verify the upper threshold set for fundamental frequency: ')
    disp('for instance, if it is set at 500 Hz, and the closing peaks do not stand out clearly')
    disp('in the derivative of the EGG signal, two periods may be detected instead of a single one; ')
    disp('one of the two will not have any opening peak, resulting in this kind of error.')
    disp('One solution is to set this threshold lower: for instance, 250 Hz for male voice.')
    disp('In the cases encountered in the development of this program, this solution was adequate.')
    disp(' ')
    rims = [];
end

