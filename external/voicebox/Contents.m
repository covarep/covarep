% Voicebox: Speech Processing Toolbox for MATLAB
%
% Audio File Input/Output
%   readwav       - Read a WAV file
%   writewav      - Write a WAV file
%   readhtk       - Read HTK waveform files
%   writehtk      - Write HTK waveform files
%   readsfs       - Read SFS files
%   readsph       - Read SPHERE/TIMIT waveform files
%   readaif       - Read AIFF Audio Interchange file format file
%   readcnx       - Raed BT Connex database files
%   readau        - Read AU files (from SUN)
%   readflac      - Read FLAC files
%
% Frequency Scales
%   frq2bark      - Convert Hz to the Bark frequency scale
%   frq2cent      - Convert Hertz to cents scale
%   frq2erb       - Convert Hertz to erb rate scale
%   frq2mel       - Convert Hertz to mel scale
%   frq2midi      - Convert Hertz to midi scale of semitones
%   bark2frq      - Convert the Bark frequency scale to Hz
%   cent2frq      - Convert cents scale to Hertz
%   erb2frq       - Convert erb rate scale to Hertz
%   mel2frq       - Convert mel scale to Hertz
%   midi2frq      - Convert midi scale of semitones to Hertz
%
% Fourier/DCT/Hartley Transforms
%   rfft          - FFT of real data
%   irfft         - Inverse of FFT of real data
%   rsfft         - FFT of real symmetric data
%   rdct          - DCT of real data
%   irdct         - Inverse of DCT of real data
%   rhartley      - Hartley transform of real data
%   zoomfft       - calculate the fft over a portion of the spectrum with any resolution
%   sphrharm      - calculate forward and inverse shperical harmonic transformations
%
% Probability Distributions
%   randvec       - Generate random vectors
%   randiscr      - Generate discrete random values with prescribed probabilities
%   rnsubset      - Select a random subset
%   randfilt      - Generate filtered random noise without transients
%   stdspectrum   - Generate standard audio and speech spectra
%   gausprod      - Calculate the product of multiple gaussians
%   maxgauss      - Calculate the mean and variance of max(x) where x is a gaussian vector
%   gaussmix      - Fit a gaussian mixture model to data values
%   gaussmixd     - Calculate marginal and conditional density distributions and perform inference
%   gaussmixk     - Estimate Kuleck-Leibler divergence between two GMMs
%   gaussmixg     - Calculate global mean, covariance and mode of a Gaussian mixture
%   gaussmixp     - Calculates and plots full and marginal probability density from a GMM
%   gmmlpdf       - Prob density function of a multivariate Gaussian mixture
%   lognmpdf      - Prob density function of a lognormal distribution
%   histndim      - N-dimensional histogram (+ plot 2-D histogram)
%   usasi         - Generate USASI noise (obsolete: use stdspectrum instead)
%
% Vector Distances
%   disteusq      - Calculate euclidean/mahanalobis distances between two sets of vectors
%   distchar      - COSH spectral distance between AR coefficient sets 
%   distitar      - Itakura spectral distance between AR coefficient sets 
%   distisar      - Itakura-Saito spectral distance between AR coefficient sets
%   distchpf      - COSH spectral distance between power spectra 
%   distitpf      - Itakura spectral distance between power spectra 
%   distispf      - Itakura-Saito spectral distance between power spectra 
%
% Speech Analysis
%   activlev      - Calculate the active level of speech (ITU-T P.56)
%   dypsa         - Estimate glottal closure instants from a speech waveform
%   enframe       - Divide a speech signal into frames for frame-based processing
%   correlogram   - calculate a 3-D correlogram
%   ewgrpdel      - Energy-weighted group delay waveform
%   fram2wav      - Interpolate frame-based values to a waveform
%   filtbankm     - Transformation matrix for a linear/mel/erb/bark-spaced filterbank from dft output 
%   fxpefac       - PEFAC pitch tracker
%   fxrapt        - RAPT pitch tracker
%   gammabank     - Calculate a bank of IIR gammatone filters
%   importsii     - Calculate the SII importance function (ANSI S3.5-1997)
%   modspect      - Caluclate the modulation specrogram
%   mos2pesq      - Convert MOS values to equivalent PESQ scores
%   overlapadd    - Reconstitute an output waveform after frame-based processing
%   pesq2mos      - Convert PESQ scores to equivalent MOS values
%   phon2sone     - Convert signal levels from phons to sones
%   psycdigit     - Experimental estimation of monotonic/unimodal psychometric function using TIDIGITS
%   psycest       - Experimental estimation of monotonic psychometric function
%   psycestu      - Experimental estimation of unimodal psychometric function 
%   psychofunc    - Psychometric functions
%   sigma         - Identify glottal closure and opening intstants from Lx or EGG waveform
%   snrseg        - Segmental SNR and Global SNR calculation
%   sone2phon     - Convert signal levels from sones to phons
%   soundspeed    - Returns the speed of sound in air as a function of temperature
%   spgrambw      - Spectrogram with many options
%   txalign       - Align two sets of time markers
%   vadsohn       - Voice activity detector
%   v_ppmvu       - Calculate the PPM, VU or EBU levels of a signal
%
% LPC Analysis of Speech
%   lpcauto       - LPC analysis: autocorrelation method
%   lpccovar      - LPC analysis: covariance method
%   lpc--2--      - Convert between alternative LPC representation
%   lpcrr2am      - Matrix with all LPC filters up to order p
%   lpcconv       - Arbitrary conversion between LPC representations
%   lpcbwexp      - Bandwidth expansion of LPC filter
%   ccwarpf       - warp complex cepstrum coefficients
%   lpcifilt      - inverse filter a speech signal
%   lpcrand       - create random stable filters
%
% Speech Synthesis
%   sapisynth     - Text-to-speech synthesis of a string or matrix 
%   glotros       - Rosenberg model of glottal waveform
%   glotlf        - Liljencrants-Fant model of glottal waveform
%
% Speech Enhancement
%   estnoiseg     - Estimate the noise spectrum from noisy speech using MMSE method
%   estnoisem     - Estimate the noise spectrum from noisy speech using minimum statistics
%   specsub       - Speech enhancement using spectral subtraction
%   ssubmmse      - Speech enhancement using MMSE estimate of spectral amplitude or log amplitude
%   specsubm      - (obsolete algorithm) Spectral subtraction 
%
% Speech Coding
%   lin2pcmu      - Convert linear PCM to mu-law PCM
%   pcma2lin      - Convert A-law PCM to linear PCM
%   pcmu2lin      - Convert mu-law PCM to linear PCM
%   lin2pcma      - Convert linear PCM to A-law PCM
%   kmeans        - Vector quantisation: k-means algorithm
%   kmeanlbg      - Vector quantisation: LBG algorithm
%   kmeanhar      - Vector quantization: K-harmonic means
%   potsband      - Create telephone bandwidth filter
%
% Speech Recognition
%   melbankm      - Mel filterbank transformation matrix
%   melcepst      - Mel cepstrum frontend for recogniser
%   cep2pow       - Convert mel cepstram means & variances to power domain
%   pow2cep       - Convert power domain means & variances to mel cepstrum
%   ldatrace      - constrained Linear Discriminant Analysis to maximize trace(W\B)
%
% Signal Processing
%   ditherq       - Add dither and quantize a signal
%   v_findpeaks   - Find peaks in a signal or spectrum (name changed to avoid conflict)
%   filterbank    - Apply a bank of IIR filters to a signal
%   maxfilt       - Running maximum filter
%   meansqtf      - Output power of a filter with white noise input
%   momfilt       - Generate running moments
%   schmitt       - Pass a signal through a schmitt trigger
%   sigalign      - Align a clean refeence with a noisy signal
%   teager        - Calculate the Teager energy waveform
%   windinfo      - Calculate window properties and figures of merit
%   windows       - Window function generation
%   zerocros      - Find interpolated zero crossings
%
% Information Theory
%   huffman       - Generate Huffman code
%   entropy       - Calculate entropy and conditional entropy
%
% Computer Vision
%   imagehomog    - Apply a homography transformation to an image with bilinear interpolation
%   polygonarea   - Calculate the area of a polygon
%   polygonwind   - Test if points are inside or outside a polygon
%   polygonxline  - Find where a line crosses a polygon
%   qrabs         - Absolute value of a real quaternion
%   qrdivide      - divide two real quaternions (or invert one)
%   qrdotdiv      - elmentwise division of two real quaternion arrays
%   qrdotmult     - elmentwise multiplication of two real quaternion arrays
%   qrmult        - multiply two real quaternion arrays
%   qrpermute     - permute the indices of a quaternion array
%   rectifyhomog  - Apply rectifing homographies to a set of cameras to make their optical axes parallel
%   rot--2--      - Convert between different representations of rotations
%   rotqrmean     - Find the average of several rotation quaternions
%   rotqrvec      - Apply a quaternion rotation to an array of 3D vectors
%   sphrharm      - forward and inverse spherical harmonic transform using uniform, Gaussian
%                   or arbitrary inclination (elevation) grids and a uniform azimuth grid.
%   upolyhedron   - Calculate the vertex coordinates and other characteristics of a uniform polyhedron
%
% Printing and Display functions
%   axisenlarge   - Selectively enlarge figure axis for clarity
%   xticksi       - Label x-axis tick marks using SI multipliers
%   yticksi       - Label y-axis tick marks using SI multipliers
%   xyzticksi     - Helper function for xticksi and yticksi
%   figbolden     - Make a figure bold for printing clearly
%   fig2emf       - Make a figure bold and save as a windows metafile
%   cblabel       - Add a label onto the colorbar
%   sprintsi      - Print a value with an SI multiplier
%   frac2bin      - Convert numbers to fixed   -point binary strings
%   v_colormap    - Set and plot colormap information
%
% Voicebox Parameters and System Interface
%   voicebox      - Global installation-dependent parameters
%   unixwhich     - Search the WINDOWS system path for an executable program (like UNIX which)
%   winenvar      - Obtain WINDOWS environment variables
%
% Utility Functions
%   atan2sc       - arctangent function that returns the sin and cos of the angle
%   bitsprec      - Rounds values to a precision of n bits
%   choosenk      - All choices of k elements out of 1:n without replacement
%   choosrnk      - All choices of k elements out of 1:n with replacement
%   dlyapsq       - Solve the discrete lyapunov equation
%   dualdiag      - Simultaneously diagonalise two hermitian matrices
%   finishat      - Estimate the finishing time of a long loop
%   fopenmkd      - like FOPEN() but creates any missing directories/folders
%   hostipinfo    - Get information about the computer name and internet connections
%   logsum        - Calculates log(sum(exp(x))) without overflow/underflow
%   minspane      - calculate the minimum (or shortest) spanning tree
%   mintrace      - find a row permutation to minimize the trace of a matrix
%   m2htmlpwd     - Create HTML documentation of matlab routines in the current directory
%   nearnonz      - Replace each zero element with the nearest non-zero element
%   permutes      - All n! permutations of 1:n
%   quadpeak      - Find quadratically-interpolated peak in a 2D array
%   rotation      - Generate rotation matrices
%   skew3d        - Generate 3x3 skew symmetric matrices
%   zerotrim      - Remove empty trailing rows and columns
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 1998-2013 Mike Brookes
%   Version: $Id: Contents.m 3758 2013-11-12 15:08:36Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

