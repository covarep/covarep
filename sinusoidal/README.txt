At one instant, a set of sinusoids is represented by a matrix 4xK where K is
the number of components. The sinusoidal parameters are along the lines:
Freq. [Hz], amp. [linear], phase [rad], harmonic number [int], peakiness [bool]

Example:
         
sins = [0    , 130.7 , 261.4;      %Frequency [Hz]
        0.01 , 0.4   , 0.35;       %Amplitude, linear
        0    , 1.57  , -3.14;      %Phase [rad]
        0    , 1     , 2    ;      %Harmonic number [integer]
        false, true  , true   ]    %True if it comes from a peak
                                   %(False if it is likely to be from noise)

For the phase, the time reference is always the center of the window.

The window length is always odd such as the sample at the center of the window
correspond to the sample of the analysis instant.

From the start of the window the window's delay is (winlen-1)/2 in nb of samples
and thus the matlab index of the central sample is (winlen-1)/2+1;
