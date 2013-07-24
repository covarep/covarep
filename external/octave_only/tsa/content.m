% Time Series Analysis (Ver 3.10)
% Schloegl A. (1996-2003,2008) Time Series Analysis - A Toolbox for the use with Matlab. 
% WWW: http://hci.tugraz.at/~schloegl/matlab/tsa/
%
%	$Id: content.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2003,2008 by Alois Schloegl <a.schloegl@ieee.org>
%
%  Time Series Analysis - a toolbox for the use with Matlab
%   aar		adaptive autoregressive estimator 
%   acovf       (*) Autocovariance function
%   acorf (acf)	(*) autocorrelation function	
%   pacf	(*) partial autocorrelation function, includes signifcance test and confidence interval
%   parcor	(*) partial autocorrelation function
%   biacovf	biautocovariance function (3rd order cumulant)
%   bispec	Bi-spectrum 
%   durlev      (*) solves Yule-Walker equation - converts ACOVF into AR parameters
%   lattice     (*) calcultes AR parameters with lattice method
%   lpc		(*) calculates the prediction coefficients form a given time series
%   invest0	(*) a prior investigation (used by invest1)
%   invest1	(*) investigates signal (useful for 1st evaluation of the data)
%   selmo	(*) Select Order of Autoregressive model using different criteria
%   histo	(*) histogram
%   hup     	(*) test Hurwitz polynomials
%   ucp     	(*) test Unit Circle Polynomials   
%   y2res	(*) computes mean, variance, skewness, kurtosis, entropy, etc. from data series 
%   ar_spa	(*) spectral analysis based on the autoregressive model
%   detrend 	(*) removes trend, can handle missing values, non-equidistant sampled data       
%   flix	floating index, interpolates data for non-interger indices
%   quantiles   calculates quantiles 
%
% Multivariate analysis (planned in future)
%   mvar	multivariate (vector) autoregressive estimation 
%   mvfilter	multivariate filter
%   arfit2	provides compatibility to ARFIT [Schneider and Neumaier, 2001]
       	
%  Conversions between Autocorrelation (AC), Autoregressive parameters (AR), 
%             	prediction polynom (POLY) and Reflection coefficient (RC)  
%   ac2poly 	(*) transforms autocorrelation into prediction polynom
%   ac2rc   	(*) transforms autocorrelation into reflexion coefficients
%   ar2rc	(*) transforms autoregressive parameters into reflection coefficients  
%   rc2ar	(*) transforms reflection coefficients into autoregressive parameters
%   poly2ac 	(*) transforms polynom to autocorrelation
%   poly2ar 	(*) transforms polynom to AR 
%   poly2rc 	(*) 
%   rc2ac 	(*) 
%   rc2poly 	(*) 
%   ar2poly 	(*) 
%   
% Utility functions 
%   sinvest1	shows the parameter calculated by INVEST1
%
% Test suites
%   tsademo		demonstrates INVEST1 on EEG data
%   invfdemo		demonstration of matched, inverse filtering
%   bisdemo		demonstrates bispectral estimation
%
% (*) indicates univariate analysis of multiple data series (each in a row) can be processed.
% (-) indicates that these functions will be removed in future 
%
% REFERENCES (sources):
%  http://www.itl.nist.gov/
%  http://mathworld.wolfram.com/
%  P.J. Brockwell and R.A. Davis "Time Series: Theory and Methods", 2nd ed. Springer, 1991.
%  O.   Foellinger "Lineare Abtastsysteme", Oldenburg Verlag, Muenchen, 1986.
%  F.   Gausch "Systemtechnik", Textbook, University of Technology Graz, 1993. 
%  M.S. Grewal and A.P. Andrews "Kalman Filtering" Prentice Hall, 1993. 
%  S.   Haykin "Adaptive Filter Theory" 3ed. Prentice Hall, 1996.
%  E.I. Jury "Theory and Application of the z-Transform Method", Robert E. Krieger Publishing Co., 1973. 
%  M.S. Kay "Modern Spectal Estimation" Prentice Hall, 1988. 
%  Ch.  Langraf and G. Schneider "Elemente der Regeltechnik", Springer Verlag, 1970.
%  S.L. Marple "Digital Spetral Analysis with Applications" Prentice Hall, 1987.
%  C.L. Nikias and A.P. Petropulu "Higher-Order Spectra Analysis" Prentice Hall, 1993.
%  M.B. Priestley "Spectral Analysis and Time Series" Academic Press, 1981. 
%  T. Schneider and A. Neumaier "Algorithm 808: ARFIT - a matlab package for the estimation of parameters and eigenmodes of multivariate autoregressive models" 
%               ACM Transactions on Mathematical software, 27(Mar), 58-65.
%  C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).
%  W.S. Wei "Time Series Analysis" Addison Wesley, 1990.
% 
% 
% REFERENCES (applications):
% [1] A. Schlögl, B. Kemp, T. Penzel, D. Kunz, S.-L. Himanen,A. Värri, G. Dorffner, G. Pfurtscheller.
%     Quality Control of polysomnographic Sleep Data by Histogram and Entropy Analysis. 
%     Clin. Neurophysiol. 1999, Dec; 110(12): 2165 - 2170.
% [2] Penzel T, Kemp B, Klösch G, Schlögl A, Hasan J, Varri A, Korhonen I.
%     Acquisition of biomedical signals databases
%     IEEE Engineering in Medicine and Biology Magazine 2001, 20(3): 25-32
%
% Features:
% - Multiple Signal Processing
% - Efficient algorithms 
% - Model order selection tools
% - higher (3rd) order analysis
% - Maximum entropy spectral estimation
% - can deal with missing values (NaN's)
