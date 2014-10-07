function [TGEI,valid,validtime,numberofpeaks] = AMPOS(SIG,rims,resampC,Ts,method,propthresh)
%%%%%%%%% Input to AMPOS (= variables that are passed on to AMPOS):
% - <SIG>: the original EGG signal
% - <rims>: indices of beginning and end of detected peaks; created by function
% RIM
% - <resampC> : coefficient for resampling: e.g., if resampC is set at 100, the
% programme will interpolate 100 times more points than there were in the
% original signal, simulating a sampling frequency 100 times higher than that
% actually used: e.g. from 44,100 Hz to an equivalent of 4,410,000 Hz.
% - <Ts>: inverse of sampling frequency: sampling period.
% - <method>: indication on the method which is to be used in the handling of 
% double / multiple peaks: selecting the first or the highest or using a
% barycentre method, etc.
% - <propthresh> : threshold for exclusion of small peaks; used only for
% glottis closures: for openings, the period is known already, and the threshold
% set at a certain height of the maximum of the signal within the period, earlier on
% in the calculations.
%
%%%%%%%%% Output of AMPOS:
% <TGEI>: 3-column matrix containing: index, amplitude, and time (in ms) of
% detected closings, following the method specified by the user for the
% approximation of Fo: in case of double peaks, the user specified (in <method>)
% selection of the first peak or of the highest peak, etc. It is not called TGCI
% (Time of Glottis-Closure Instants) but only TGEI (Time of Glottal-Event
% Instant) because this function applies to opening instants as well.
% <valid>: created by PEAKSHAPE: provides an indication on the relative heights
% of the two highest summits within the peak.
% This function reinterpolates the original EGG signal to obtain a very accurate
% measurement of the amplitude of the peak. This follows advice given by Wolfgang Hess, 
% Univ. of Bonn, during discussions on the amplitude of the peak of the derivative
% at closure (DECPA) at the Speech Prosody 2004 conference, Nara, Japan. In the
% case of opening peaks, reinterpolation is not recommended, unless perhaps when visual
% observation shows exceptionally clear, unique opening peaks.

%%%%%%% finding out the amplitude and position of the peaks

% condition for cases where no single crossing was detected
EX = isempty(rims);
if EX == 1
    TGEI = [];
    valid = [];
    validtime = [];
    numberofpeaks = 0; 
else
	% Loop for individual peaks
	for z = 1: length(rims(:,1))
        peakSIG = []; timevector = []; REI = []; dREI = [];
        % retrieving corresponding portion of EGG signal
        peakSIG = SIG(rims(z,1):rims(z,2));
		% computing number of samples in the peak
		LP = rims(z,2) - rims(z,1) + 1;
        
        % if a coefficient for interpolation ( > 1) has been specified: 
        if resampC > 1
            % creating vector containing <resampC> (for: RESAMPLINGCoefficient) times more points, 
            % for spline interpolation
            timevector = [1:resampC * LP];
            timevector = timevector / resampC;
            % cubic spline interpolation of the (original) EGG signal
			REI = interp1(1:LP,peakSIG,timevector,'spline');
        else
            % otherwise: signal unmodified.
            REI = peakSIG;
        end
        
        % derivation of signal
        for k = 1:length(REI) - 1
            dREI(k) = REI(k + 1) - REI(k);
        end
	
        % multiplication of dREI by the <resampC> coefficient: otherwise the values will be much
        % lower than those in <SMOO_dSIG>, and cannot be plotted together with these
        % values in tests.
        dREI = resampC * dREI;
	
        % analysis of the shape of the peak: calling function PEAKSHAPE, which
        % returns the number of summits in peak, and the GCI and DECPA retrieved by
        % the method specified
        [TGEI(z,1),TGEI(z,2),valid(z),validtime(z),numberofpeaks(z)] = peakshape(dREI,method,resampC,Ts,propthresh);
        % converting the index into a relative time value: divide it by resampC, so
        % as to obtain a value that can be added to the index of the point where the
        % peak begins: 
        TGEI(z,1) = TGEI(z,1) / resampC;
        % The index is given relative to the portion of the vector that has been
        % passed on to the MAX function, e.g. 10 if the peak corresponds to the 10th
        % sample within dREI. To retrieve the index in
        % absolute terms (relative to the whole file), the index of the point 
        % where the peak begins is added: 
        TGEI(z,1) = TGEI(z,1) + rims(z,1);
        % translating this value into a time value, in seconds; this value is placed
        % in the 3rd column of the matrix <TGEI> :
        TGEI(z,3) = TGEI(z,1) * Ts;
	% end of loop for individual peaks
	end
end