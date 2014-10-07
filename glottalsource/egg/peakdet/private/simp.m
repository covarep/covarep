function [simppeak,BUND] = simp(TGEI,doublepeaks,method)
% On the basis of TGEI, creating 
% a matrix with only one maximum per cycle. 
%
% Input: 
% - TGEI matrix created by AMPOS
% - doublepeaks vector, indicating, for each peak in TGEI, whether it is a single
% peak (doublepeaks(i) == 0), the beginning of a peak cluster (==1), midway
% through a peak cluster (==2), or the end of a peak cluster (==3).
% - "method" variable: set at 0 if the highest of the peaks is to be selected,
% set at 1 if the first of the peaks is to be selected.
%
% Output: 
% - simppeak matrix, made up of two columns: one time (in seconds) and one
% amplitude, corresponding to the closing (resp. opening).
% Calculation of Fo in case of double closing peaks were perhaps
% best done by other methods, such as autocorrelation, when these can apply; the
% choice of the highest peak, the first peak, or a pondered value in-betweeen
% is specified in the input to
% the function: variable "method".
% - BUND vector, containing the number of peaks in each detected peak bundle:
% e.g. BUND = [2 2 2] indicates that three bundles of two peaks were detected.

%%% finding out peak bunches
% Variable to count Number of Real Peaks (one peak per cycle, excluding 
% double peaks)
NRP = 0; 
% Variable to count Number of Peak Bundles
NPB = 0;
BUND = []; simppeak = []; i = 1;
while i <= length(TGEI(:,1))
    %% disp(['Treating line number ',num2str(i),' of TGEI matrix.'])
    % In case the peak is unique:
    if doublepeaks(i) == 0
        % incrementing the count of real peaks
        NRP = NRP + 1;
        % placing the values in simppeak
        simppeak(NRP,:) = [TGEI(i,3) TGEI(i,2)];
        % incrementing i, so that the next item will be treated
        i = i + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% In case the peak is not unique: entering another loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    elseif doublepeaks(i) == 1
        INLOOP = 1;
        % counter of number of peaks in cluster:
        NBPEAKS = 0;
        % matrix containing all peaks in cluster: 
        PG = [];
        % incrementing counter of peak clusters
        NPB = NPB + 1;
        while INLOOP == 1
            %  If the value is 1 or 2: one is at the beginning or middle of a
            %  peak cluster.
            if ismember(doublepeaks(i),[1 2])
                % incrementing counter of peaks in cluster
                NBPEAKS = NBPEAKS + 1;
                % inserting peak in peak-cluster matrix. Variable PG: for Peaks
                % Group. The values are: time (in s), amplitude.
                PG(NBPEAKS,:) =  [TGEI(i,3) TGEI(i,2)];
                % incrementing i, to move on to the next peak
                i = i + 1;

            elseif doublepeaks(i) == 3
                % incrementing number of peaks in bundle
                NBPEAKS = NBPEAKS + 1;
                % changing INLOOP, so as to exit the loop
                INLOOP = 0;
                % inserting peak in peak-bundle matrix. Variable PG: for Peaks Group
                PG(NBPEAKS,:) =  [TGEI(i,3) TGEI(i,2)];
                % incrementing i, to pass on to the next peak
                i = i + 1;

            elseif doublepeaks(i) == 0
                % changing INLOOP, so as to exit the loop
                INLOOP = 0;

            end
        end
        
        % After the "while" loop has been exited, all the peaks in the peak
        % cluster are in the PG matrix. One single value must be arrived at, 
        % following the method specified by the user (in the 
        % "method" variable), and placed in the simppeak matrix.
        if method == 0
            % re-initializing variable <selectedpeak>
            selectedpeak = [];
            % method "0": selecting the highest peak
            % recovering the value and index of the highest peak
            [hiDECPA,hiindex] = max(PG(:,2));
            % assigning the DEEPA value to "selected peak"
            selectedpeak(2) = hiDECPA;
            % assigning the corresponding time point to "selected peak"
            selectedpeak(1) = PG(hiindex,1);
        elseif method == 1
            % method "1": selecting the first peak
            % re-initializing variable <selectedpeak>
            selectedpeak = [];
            % selecting the first peak:
            selectedpeak(1:2) = PG(1,:);
        elseif method == 2
            % method "2": selecting the last peak
            % re-initializing variable <selectedpeak>
            selectedpeak = [];
            % selecting the last peak:
            selectedpeak(1:2) = PG(length(PG),:);
        elseif method == 3
            % method "3": barycentre method
            %%%%%% calculation of the amplitude ratio, and the distance, between 
            %%%%%% highest peak and each of the other peaks
            %%% Retrieving the maxima: 
            [PA,indexMAX] = max(PG(:,2));
		    
            % retrieving time value <TIM> on the basis of <indexMAX>, which is simply
            % an index within the column in <PG>
            TIM = PG(indexMAX,1);

            % suppressing the maximum from <PG>
            TE = [];
			if  indexMAX - 1 > 0
                for j = 1: (indexMAX - 1)
                    % TE: for TEmporary variable
                    TE (j,1:2) = PG(j,1:2);
                end
			end
			if indexMAX + 1 <= length(PG(:,1))
                for j = (indexMAX + 1) : length(PG(:,1))
                    TE(j - 1,1:2) = PG(j,1:2);
                end
			end
			PG = [];
            PG = TE;        
		    
            % for coefficiented mean: using a WEIGHT variable, set at 1 to
            % begin with, and then augmented as more peaks are added to the mean,
            % one by one.
            WEIGHT = 1;
		
            for k = 1:length(PG(:,1))
                % time distance between the peak being examined and the highest peak
                dist(k) = PG(k,1) - TIM;
                % This variable is
                % not an absolute value (abs(...)): it is negative if the peak being
                % examined comes before the main peak, positive if it comes after.
                
                % amplitude ratio between the peak being examined and the highest peak
                % (in %)
                ampratio(k) = PG(k,2) / PA;
                
                % halving <dist(k)>, and adding a coefficient, to obtain the
                % distance to be subtracted from I
                dist(k) = (dist(k) / 2) * ampratio(k);
                % dividing by WEIGHT, to take into account all the peaks that
                % are averaged
                dist(k) = dist(k) / WEIGHT;
                % increasing WEIGHT
                WEIGHT = WEIGHT + ampratio(k);
                %%%%% correcting the peak position according to the coefficiented distance
                TIM = TIM + dist(k);
                % end of the loop for method 3    
            end
            % passing on the results to the <selectedpeak> variable, returned by
            % the function
            selectedpeak(1) = TIM;
            % The amplitude value (called DECPA for closing peaks, DEOPA for
            % opening peaks; common term for the two used here in variable
            % names: DEEPA, Derivative-Electroglottographic Event Peak Amplitude)
            % is unchanged: it remains set at the peak amplitude of the highest maximum.
            selectedpeak(2) = PA;
                    
                    
            %%%%%%%%% Note: 
            % No special condition if method == 4 : in that case, the double/
            % multiple peaks are simply excluded from calculations; no selection
            % of a compromised value.
		
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% End of the loop for non-unique peaks
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        % incrementing the number of "real" peaks
        NRP = NRP + 1;
        % assigning the selected peak in the results file
        simppeak(NRP,:) = selectedpeak;
        % storing the details on the peak bundle that has just been
        % treated
        BUND(NPB) = NBPEAKS;
    % end of the condition on peak unicity        
    end
% end of "while" loop
end

