

                   Harmonic Model + Phase Distortion (HMPD)
      Copyright (c) 2012 University of Crete - Computer Science Department
                                October 2014


Author
 Gilles Degottex <degottex@csd.uoc.gr>


TODO


Daniel (default):
opt.pdm_nbper = 6;
opt.pdd_nbper = 2;
Nico:
opt.pdm_nbper = 3;
opt.pdd_nbper = 3;


--------------------------------------------------------------------------------
F0 estimator

Switched to SRH for legal reasons
DIFFERENCE WITH THE PUBLISHED WORKS
For reason of fairness in the evaluations in [TODO], STRAIGHT[TODO] was used to
compute the f0 values. However, it is not possible to include STRAIGHT
in COVAREP for license reasons. Consequently, in the present version of
HMPD, the SHR[TODO] is used.

--------------------------------------------------------------------------------
Compression

The amplitude, PDD and PDM features can be compressed on a log scale, using the
following options:
    hmpdopt.amp_enc_method=2; hmpdopt.amp_log=true; hmpdopt.amp_order=39;
    hmpdopt.pdd_log=true; hmpdopt.pdd_order=32;
    hmpdopt.pdm_log=true; hmpdopt.pdm_order=32;
(commented in the HOWTO_hmpd.m file)

Note that the compressed PDM has never been used for synthesis, only for
classification [TODO].


--------------------------------------------------------------------------------
Computational speed

The analysis step is very slow, mainly because of the harmonic analysis which
estimates sinusoidal parameters at harmonic frequencies.
It has been shown that analysis/re-resynthesis of the aHM model with LS solution
(aHM-LS) has quasi-perfect reconstruction while providing among the most
accurate sinusoidal parameters [TODO]. Therefore, aHM-LS was used in the
original works in order to minimize the impact of the harmonic analysis on the
evaluation results. This ensures, as best as possible, that the results
presented in [TODO,TODO] are due to the phase processing techniques, the studied
subject, and not the harmonic analysis method.
However, this doesn't mean that the results cannot be obtained with a simpler,
and quicker, harmonic analysis method.
Thus, the analysis step can be speed up by using a simple Peak Picking (PP)
method, without adaptation of the frequency basis:
    hmpdopt.sin.use_ls     = false;
    hmpdopt.sin.fadapted   = false;
(commented in the HOWTO_hmpd.m file)

Implications/Remarks:
    * It should drop the computation time to ~10% of the original
      time.
    * WARNING: Keep in mind that using the peak picking method instead of the
      aHM-LS solution, the quality is degraded! As shown in [TODO].
      (e.g. fricatives sound more bussy using PP than using aHM-LS).
    * It could be a good practice to use first the above options for testing a
      prototype whereas final evaluations should be done using the LS solution
      (i.e. without these speed up options).
    * For reason of reproduction of the published results, the slow settings are
      used in the HOWTO_hmpd.m file.


A second option to speed up the computational time is to use the mex function
interp1ordered, which is a C implementation of the linear interpolation (without
any input checking and way faster than interp1q).
The slow interp1 function can be replaced by activating the option:
    hmpdopt.usemex = true;
(commented in the HOWTO_hmpd.m file)

Implications/Remarks:
    * It should drop the computation time of the synthesis to ~35% of the
      original time, making it real-time for a state-of-the-art machine (2014).
    * This drop the computation time of the feature computation to only ~75%.
    * This replace some spline interpolations by linear interpolation.
      (e.g. RPS interpolation, thus, creating steps in the frequency tracks
            of the harmonics).
      Though, to my best experience, I've never heard any difference.
    * Similarly to the options above, a good practice would be to develope
      prototypes with the speed up optins of HMPD, but the final results should
      be obtained without using any options, as in [TODO].



--------------------------------------------------------------------------------
Subjects to investigate

* Obtain the quality obtained with aHM-LS by using a Peak Picking method,
  which is way faster.
* Improve the accuracy of PDD's estimate.
  (and study the tradeoff for opt.sin_nbat about quality and speed).
* Estimate an f0 curve which is convenient for HMPD.
  It should favor low frequency values in unvoiced segments.


