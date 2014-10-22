% Analysis of the electroglottographic signal for a continuously voiced
%  item of normal speech: typically, one syllable or a portion of a syllable.
%
% Description
%  This function was devised for SPOKEN materials, not for SUNG materials:
%  it was found that Nathalie Henrich's <decom.m> MatLab function, devised for
%  the singing voice, did not yield results for speech data
%  (specifically: Vietnamese syllables with strong glottalization) because
%  the EGG signal failed to meet the criterion for quasi-periodicity. Hence the
%  decision to develop a new function, which, like <decom.m>, detects closing
%  peaks and opening peaks on the derivative of the EGG signal, but by a
%  *threshold* method and not by autocorrelation [1,2].
%
%  Created in October 2004, with minor changes in July 2005 and July 2007, by
%  Alexis Michaud. Adapted to COVAREP by Nguyen Thi Lan.
%
% Inputs
%  SIG: The input EGG signal, in VFCA convention (vocal-fold contact area), 
%       i.e. EGG increases with glottal contact. The peak at closure is thus
%            positive in the EGG derivative.
%  FS : [Hz] The sampling frequency of the EGG signal.
%  F0 : The maximum possible F0 (e.g. F0=500Hz).
%  Method: Method chosen to handle double closing peaks in F0 calculation:
%                 if <method == 0>: selecting highest of the peaks
%                 if <method == 1>: selecting first peak
%                 if <method == 2>: selecting last peak
%                 if <method == 3>: using barycentre method
%                 if <method == 4>: exclude all double peaks
%       
%          The duration of glottal cycles is measured by detecting positive
%          peaks on the derivative of the electroglottographic signal (hereafter
%          DEGG). The inverse of this duration is still called
%          "fundamental frequency", but there is no check on periodicity, as
%          peaks are detected individually. This is different from correlation-
%          based methods (such as DECOM: Henrich et al. 2004). Thus this
%          algorithm provides indications on glottal cycles even for
%          laryngealized ("creaky") portions of signal.
%
% Outputs
%  results_matrix : Single matrix containing all the results.
%               Matrix content: 
%              - beginning and end of period : in 1st and 2nd columns
%              - F0 : 3rd column
%              - DECPA, Derivative-Electroglottographic Closure Peak Amplitude:
%                       4th column (on DECPA and DEOPA: see [3])
%              - Oq determined from raw maximum, and DEOPA : 5-6th col [1,2]
%              - Oq determined from maximum after smoothing : 7th col [1,2]
%              - Oq determined from peak detection : 8-9th col without
%                smoothing and with smoothing, respectively [1,2].
%
% Example
%  See the HOWTO_egg.m example file.
%
% Reference
% [1] Martine Mazaudon and Alexis Michaud, "Tonal Contrasts and Initial
%     Consonants: A Case Study of Tamang, a 'missing Link' in Tonogenesis",
%     Phonetica 65 (4): 231-56, 2008.
% [2] Alexis Michaud "Final Consonants and Glottalization: New Perspectives from
%     Hanoi Vietnamese", Phonetica 61 (2-3): 119-46, 2004.
% [3] Michaud, Alexis. "A Measurement from Electroglottography: DECPA, and its
%     Application in Prosody". In Bernard Bel & Isabelle Marlien (eds.), Proc.
%     Speech Prosody 2004, 633-636. Nara, Japan.
% [4] Guide available online at:
%      http://voiceresearch.free.fr/egg/softwares.htm#peakdet
% 
% Copyright (c) 2004 CNRS (Centre National de la Recherche Scientifique, France)
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
% Authors
%  Alexis Michaud <alexis.michaud@vjf.cnrs.fr> <michaud.cnrs@gmail.com>
%                  CNRS (Centre National de la Recherche Scientifique, France)
%  Nguyen Thi Lan <thi-lan.nguyen@mica.edu.vn>
%                 International Research Institute MICA, September 2014.
%                 Adaptation to COVAREP standards
%

function [ results_matrix, SdSIG, SIG ] = peakdet ( SIG, FS, F0, Method )

% Setting the resampling coefficient. The electroglottographic signal is
%  resampled (reinterpolated) at the closing and opening peaks for accurate
%  peak detection; the coefficient for reinterpolation is set at 100 by default;
%  this can be modified by the user, to define it as a function of the sampling
%  frequency of the original signal, or in relation to the fundamental frequency
%  of the sample under analysis.
resampC = 100;

% Setting the coefficient for recognition of "double peaks":
%  for closings, is 0.5, following Henrich N., d'Alessandro C., Castellengo M.
%  and Doval B., 2004, "On the use of the derivative of electroglottographic
%  signals for characterization of non-pathological voice phonation", Journal
%  of the Acoustical Society of America, 115(3), pp. 1321-1332.
propthresh = 0.5;

% Choice of method chosen to handle double closing peaks in Fo calculation
method = Method;

% Choosing maximum possible F0. To avoid absurd results in the case of double
% closing peaks, a threshold on maximum fundamental frequency is set by the
% user, and peaks that are so close that the corresponding F0 is above this
% threshold are considered as belonging in the same "peak cluster".
maxF = F0;
% In light of the great differences in F0 range across speakers and across the
% experimental tasks that they perform, it did not appear adequate to set the
% same threshold for all speakers (say, 500 Hz): the user must be allowed to
% modify the threshold.
% It was not found useful to set a lower threshold parallel to this upper
% threshold. Implausibly low values in the results point to one of the following
% situations :
%  - One closing peak has been detected before the onset of voicing or after the
%    offset of voicing, resulting in the detection of a ‘period’ the inverse of
%    which is under 20 Hz. These cases can be corrected by suppressing the first
%    or last period ; this option is offered by the program, at the stage where
%    the user is asked to check the results.
%  - Some closing peaks within a voiced portion of signal have gone undetected
%    because their amplitude is below the threshold. The user must then check on
%    the figure which amplitude threshold is to be chosen for all the peaks to
%    be detected, and set the amplitude threshold accordingly ; this option is
%    offered by the program when the user is asked to confirm the results.

% Choosing the smoothing step. Smoothing the DEGG signal turned out to be useful
% in the many cases where there is one single opening peak but its amplitude is
% very small and it tends to be drowned in noise. A smoothing step of 1 means
% smoothing 1 point to the left and right, i.e. each point in <SMOO_dSIG> is the
% average of 3 points in <dSIG> ; 2 means averaging over (2x2+1) = 5 points.
% As the programme computes the open quotient results by four methods two of
% which operate on the unsmoothed signal, it is avisable to choose a smoothing
% step (of 1) even if the user believes that this smoothing is unnecessary :
% this assumption can be verified by comparing the results with and without
% smoothing, which will show a complete fit if the original signal has very
% little background noise.
% Concerning the choice of a smoothing step for noisy signals : a step up to 5
% can be chosen ; visual comparison of the DEGG signals before and after
% smoothing is recommended to verify that this smoothing does not make
% neighbouring peaks coalesce. It must be remembered that in some cases the DEGG
% method simply does not apply, and (arguably) should not be forcibly applied:
% if taken to its limits, smoothing artificially creates a neat hump for opening
% and one for closing, but these humps fudge up the issue, as they do not
% correspond to any precise physiological reality anymore : the advantage of the
% DEGG method is that it is based on an established relationship between the
% DEGG signal and significant glottal events ; extreme smoothing blurs this
% relationship.
% In a nutshell : a smoothing step of 1 is adequate for high-quality signals
% (which already appear visually as very smooth), a smoothing step of 2 or 3
% increases correct peak detection in relatively noisy signals. Fudging-up of
% opening peaks was only observed with a smoothing step of 6 or more. 
smoothingstep = 3;

% assigning default values to <COEF> vector, used in FO: sampling frequency of
% the electroglottographic recording;
% smoothing step specified by user; 1, to indicate that the amplitude threshold
% for peak detection will be set automatically; the fourth value is the
% threshold value set manually; left at 0 to begin with, as the value will be
% set automatically and not manually.
COEF = [FS smoothingstep 1 0];

%%%%%%%%%%%%%% running main analysis programme
[Fo,Oq,Oqval,DEOPA,goodperiods,OqS,OqvalS,DEOPAS,goodperiodsS,simppeak,SIG,dSIG,SdSIG] = FO(COEF,method,propthresh,resampC,maxF,SIG,FS);	
        %%% Placing main results in a single matrix
        datafile = [];
        if isempty(Fo)
        else
            for k = 1:length(Fo)
                datafile(k,1) = simppeak(k,1);
                datafile(k,2) = simppeak(k + 1,1);
                datafile(k,3) = Fo(k);
                datafile(k,4) = simppeak(k,2);
                datafile(k,5) = Oq(k);
                datafile(k,6) = DEOPAS(k);
                datafile(k,7) = OqS(k);
                datafile(k,8) = Oqval(k);
                datafile(k,9) = OqvalS(k);
            end
        end
        %%%%%%%%%%%%%% placing results in matrices

          % checking that there is no doubling of the last line (this occasional
          % problem results in a bug that I have not identified, which causes
          % the last line to be written twice into the <datafile> matrix)
          ld = length(datafile(:,1));
          if ld > 1
              if datafile(ld,:) == datafile(ld - 1,:)
                  datafile = datafile(1:ld - 1,:);
              end
          end

          % calculating the number of periods (= nb of lines)
          period_nb = size(datafile,1);

          % calculating the number of columns
          nbcol = length(datafile(1,:));
        % end of the condition on non-emptiness of Fo variable 
   % assign returning result to returning matrix
   results_matrix=datafile;
   
end

