function [rims] = CRO(dSIG,COEF);
%%%%%%%%%%%%%%%%%%% CRO : Detection of positive peaks in an electroglottographic
%%%%%%%%%%%%%%%%%%% signal, yielding the index of the points where the threshold
%%%%%%%%%%%%%%%%%%% (set at: absolute value of minimum in negative part of the
%%%%%%%%%%%%%%%%%%% DEGG signal) is crossed.
% Syntax : [rims] = CRO(dSIG,COEF)
%%%%%%%%%%%%%%%% Input to GCI : 
% - <dSIG>: the relevant portion of the derivative of the electroglottographic signal.
% - <COEF>: vector containing four values: 
% 1) sampling frequency ; 2) smoothing step requested by user ; 
% 3) toggle set by user : if == 1 : automatic setting of amplitude threshold;  
% if ==0: equal to the value indicated in the 4th column (4th value in vector).
%%%%%%%%%%%%%%%% Output of GCI :
% <rims>: vector containing the points where the threshold is crossed; see last
% line of this function, and function RIM.
%%%% Convention: positive values for closing peaks, negative values for opening peaks.

% Sampling period: (T for: period, s for: sampling) 
Ts = 1 / COEF(1);

% retrieving user choice for amplitude threshold setting method
toggle = COEF(3);

% finding out the intervals where closing peaks are likely to be found, 
% by detecting the left and right boundaries of peaks: 
% where the signal crosses the threshold.
% The threshold is chosen on the basis of: the highest negative peak (highest opening peak), 
% i.e. the minimum of vector. Its sign is opposite.
% thresh = - min(dSIG);
% A coefficient may be added, increasing the threshold by a certain coefficient to
% limit problems with "false peaks". Here, the threshold is not increased: 
% double peaks are spotted and treated later on.
% But conversely, if there are very small closing peaks, then this threshold may not be proper. 
if toggle == 1
	disp('Setting threshold for detection of closing peaks based on highest opening peak.')
    thresh =  - min(dSIG);
    % extra condition: it may happen (e.g. speaker F2 of the Naxi corpus) that
    % opening peaks are actually as high as closing peaks; in such cases,
    % setting the threshold at - min(dSIG) results in lack of detection of
    % numerous closing peaks. Solution: if the threshold is higher than a certain proportion 
	% of the highest closing peak, the threshold is calculated again.
    % This proportion was placed at one-third, then (realizing it was not enough) at one-fourth, 
    % and (in 2019) at one-ninth, then one-tenth, one-eleventh and one-twelfth. 
	
	% % For debugging:
	% disp('threshold: minimum of signal:')
	% disp(num2str(thresh))
	% disp('maximum positive peak:')
	% disp(num2str(max(dSIG)))
	% disp('divided by 5: ')
	% disp(num2str(max(dSIG) / 12))
    if thresh > (max(dSIG) / 12)
        thresh = max(dSIG) / 12;
    end
else
% if the method is manual choice: retrieving user choice for amplitude threshold 
    thresh = COEF(4);
end

% finding out the edges of peaks : yields a matrix of indices: 1st column: crossing the threshold, towards ceiling; 2nd column: 
% crossing the threshold towards the floor.
rims = rim (dSIG,thresh);