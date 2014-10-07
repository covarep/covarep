function doublepeaks = detectmult(TGEI,maxF)
%%%% Function for the detection of double / multiple peaks
%
% Input to the function: 
% - matrix TGEI, containing all the peaks detected by the AMPOS function; 
% - a toggle, set to 1 if a maximum plausible Fo value (maxF) is passed on 
% to the function, and set to 0 if the value maxF is left at its default
% value, which is 500 Hz
% - the value maxF if desired.
%
% 
% Output of the function: 
% for exclusion of the double / multiple peaks from the calculations: a vector
% of the same length as TGEI; for each peak, doublepeaks(i) is set at 0 if the
% peak is unique, and at 1, 2 or 3 if it is part of a peak cluster (= if it is not alone in 
% its glottal cycle): 1 if it is the beginning of a peak cluster; 2 if it is
% midway in a peak cluster; 3 if it is the end of a peak cluster.

%%%%%%%%%%% Calculation %%%%%%%%%%%%%%

% Setting all values in doublepeaks at 0
doublepeaks(length(TGEI(:,1))) = 0;

% condition to avoid crashes in case TGEI only has 1 single line
if length(TGEI(:,1)) > 1
	for z = 1: length(TGEI(:,1)) - 1
        % If the "period" has a frequency > maxF, both peaks are marked as
        % double / multiple. 
        % If this is the first peak in the matrix, or the previous peak has not been
        % marked as double, the peak of index z is the beginning of a peak cluster:
        % assignment of the value "1".
        if 1 / ( TGEI(z + 1,3) - TGEI(z,3) ) > maxF
            if z == 1
                % assigning value 1 to slot z
                doublepeaks(z) = 1;
                doublepeaks(z + 1) = 3;
            elseif ismember(doublepeaks(z - 1),[0 3])
                % assigning value 1 to slot z
                doublepeaks(z) = 1;
                doublepeaks(z + 1) = 3;
            elseif ismember(doublepeaks(z - 1),[1 2])
                doublepeaks(z) = 2;
                doublepeaks(z + 1) = 3;
            end
        % No "else" is needed, as the whole vector has been filled with zeros to
        % begin with: in case the value to be assigned is neither 1 nor 2 nor 3,
        % it will be 0 by default.
            
        % end of "if" condition
        end

        % end of the "for" loop
	end

% end of the condition on the size of TGEI    
end

% No special case is needed for the last peak: when the function deals with the
% point before last, it will consider the last period, and label the last point
% as "3" if this period is above the threshold.