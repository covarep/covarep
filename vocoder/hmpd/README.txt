
               Harmonic Model + Phase Distortion (HMPD) vocoder

      Copyright (c) 2012 University of Crete - Computer Science Department \
               Foundation for Research and Technology-Hellas - \
                  Institute of Computer Science (FORTH-ICS)

                    Gilles Degottex <degottex@csd.uoc.gr>
                                October 2014


This file describes mainly the main options available in HMPD and provides some
remarks on its use.

A use case is given in the HOWTO_hmpd.m file in the howtos directory.


F0 ESTIMATION ------------------------------------------------------------------

For reason of fairness in evaluations, STRAIGHT [6] was used in [1,2] to compute
the f0 values in voiced segments. However, it is not possible to include
STRAIGHT in COVAREP for legal reasons. Consequently, the version of HMPD in
COVAREP uses the SRH [7] method instead.

You can obviously use your own f0 estimation curve given as argument of the 
hmpd_analysis function (or similarly hmpd_analysis_harmonic). For this
purpose, please see the documentation in these functions.

If no f0 curve is provided, it is necessary to set f0min and f0max options
appropriately with respect to the analyzed voice, as shown in HOWTO_hmpd.


COMPUTATION SPEED --------------------------------------------------------------

The analysis step is very slow, mainly because of the harmonic analysis which
estimates sinusoidal parameters at harmonic frequencies (it accounts for ~95%
of the computation time of the analysis).

It has been shown [5] that analysis/re-resynthesis of the aHM model with LS
solution (aHM-LS) has quasi-perfect reconstruction while providing among the
most accurate sinusoidal parameters [5]. Therefore, aHM-LS was used in the
original works of HMPD [1,2,3,4] in order to minimize the impact of the harmonic
analysis on the evaluation results. This ensures, as best as possible, that the
results presented in [1,2,3,4] are due to the phase processing techniques, the 
subject studied, and not due to the harmonic analysis method.
However, this doesn't mean that the results cannot be obtained with a simpler,
and quicker, harmonic analysis method.
Thus, the analysis step can be speed up by using a simple Peak Picking (PP)
method, without adaptation of the frequency basis:
    hmpdopt.sin.use_ls     = false;
    hmpdopt.sin.fadapted   = false;
(commented in the HOWTO_hmpd.m file)

Implications/Remarks:
    * It should drop the computation time of the harmonic analysis to ~10% of
      the original time.
    * WARNING: Keep in mind that using the peak picking method instead of the
      aHM-LS solution, the quality is degraded! (as shown in [5])
      (e.g. fricatives sound more buzzy using PP than using aHM-LS).
    * It could be a good practice to use the speed up options for testing a
      prototype, whereas final evaluations should be done using the LS solution.
    * For reason of replication of the published results, the slow settings are
      used in the HOWTO_hmpd.m file.

By looking at the content of hmpd_analysis.m, you will notice that this function
simply calls the two functions hmpd_analysis_harmonic and hmpd_analysis_features.
Thus, if you only want to play around with the computation of the features, you
can replace hmpd_analysis by it content in the HOWTO_hmpd file:
    % Estimate sinusoidal harmonic parameters
    frames = hmpd_analysis_harmonic(wav, fs, f0s, opt)
    % Compute amplitude envelope and phase statistics
    [f0s, AE, PDM, PDD] = hmpd_analysis_features(frames, fs, opt);
and save the intermediate frames structure instead of recomputing it each time
you run hmpd_analysis.


A second solution to speed up the computational time is to use the mex function
interp1ordered, which is a C implementation of the linear interpolation (without
any input checking and way faster than interp1q).
The slow interp1 function can be replaced by activating the option:
    hmpdopt.usemex = true;
(commented in the HOWTO_hmpd.m file)

Implications/Remarks:
    * It should drop the computation time of the synthesis to ~35% of the
      original time, making it real-time for a state-of-the-art machine (2012).
    * It almost does NOT affect the computation time of the analysis step.
    * This replaces some spline interpolations by linear interpolation.
      (e.g. RPS interpolation). Thus, creating steps in the frequency tracks
      of the harmonics. Though, to my best experience, I've never heard any
      difference.
    * Similarly to the previous option above, a good practice would be to
      develop prototypes with the speed up options, but the final results should
      be obtained without using any such options, as in [1,2,3,4].


COMPRESSION --------------------------------------------------------------------

The amplitude envelope and the short-term Phase Distortion statistics (PDD and
PDM) can be compressed using log scales by using the following options:
    hmpdopt.amp_enc_method=2; hmpdopt.amp_log=true; hmpdopt.amp_order=39;
    hmpdopt.pdd_log=true; hmpdopt.pdd_order=12; % 12 in [1,2], 24 in [3]
    hmpdopt.pdm_log=true; hmpdopt.pdm_order=24; % unused in [1,2], 24 in [3]
(commented in the HOWTO_hmpd.m file)

Note that the compressed PDM has never been used for synthesis, only for
classification [3].


PARAMETERS DIFFERENCES AMONG THE PUBLICATIONS ----------------------------------

Between the publications [1,2,3,4], the options used were not exactly the same.
The HMPD code in COVAREP with the default options implements the version used
in [1,2].

Roughly, for computing the Phase Distortion Deviation (PDD), the trend of the
phase of the pulse's shape is not removed in either [3] or [4]. In [3,4] PDD is
neither corrected, as described in [1,2].
For further details, please refer to the papers for the differences with the
other publications.


ERRATA -------------------------------------------------------------------------

A)  Unfortunately, the publications [1,2] do not mention a median filter and a
    hanning filter used when computing the smooth PD (the local trend).
    This has been forgotten during the writing of [1,2]. Please have a look at
    hmpd_phase_smooth.m for the implementation used.

B)  In [3], the description of the PDM compression is wrong. We tried the
    compression described in [3], which didn't work as expected,
    Later on, I replaced it by the one described in philin2philog.m, but I
    didn't replace the description in the paper we were writing.

Thanks to the COVAREP project, you can work with the implementations that have
been actually used for the experiments presented in [1,2] and not with the
article descriptions only.



REFERENCES ---------------------------------------------------------------------

[1] G. Degottex and D. Erro, "A uniform phase representation for the harmonic
    model in speech synthesis applications", EURASIP, Journal on Audio, Speech,
    and Music Processing - Special Issue: Models of Speech - In Search of Better
    Representations, 2014.
[2] G. Degottex and D. Erro, "A Measure of Phase Randomness for the Harmonic
    Model in Speech Synthesis", In Proc. Interspeech, Singapore. International
    Speech Communication Association (ISCA), September 2014.
[3] G. Degottex and N. Obin, "Phase Distortion Statistics as a Representation of
    the Glottal Source: Application to the Classification of Voice Qualities",
    In Proc. Interspeech, Singapore. International Speech Communication
    Association (ISCA), September 2014. 
[4] M. Koutsogiannaki, O. Simantiraki, G. Degottex and Y. Stylianou, "The
    Importance of Phase on Voice Quality Assessment", In Proc. Interspeech,
    Singapore. International Speech Communication Association (ISCA), September
    2014.
[5] G. Degottex and Y. Stylianou, "Analysis and Synthesis of Speech using an
    Adaptive Full-band Harmonic Model", IEEE Transactions on Acoustics, Speech
    and Language Processing, 21(10):2085-2095, 2013.
[6] H. Kawahara, I. Masuda-Katsuse and A. de Cheveigne, "Restructuring speech
    representations using a pitch-adaptative time-frequency smoothing and an
    instantaneous-frequency-based f0 extraction: Possible role of a repetitive
    structure in sounds", Speech Communication 27(3-4), 187-207, 1999.
[7] T. Drugman and A. Alwan, "Joint Robust Voicing Detection and Pitch
    Estimation Based on Residual Harmonics", Interspeech, Firenze, Italy, 2011.
