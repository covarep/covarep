% Cepstral Peak Prominence (CPP) for discriminating breathy from
% neutral/modal phonation.
%
% Octave compatible
%
% Description
% The Cepstral peak promience feature derived from the log magnitude
% cepstrum based on the paper by Hillenbrand & Houde (1996). The code
% facilitates time and quefrency smoothing variants, and regression 
% line (as per Hillenbrand), mean and no normalisation variants
%
% Inputs
%  x          : [samples]    [Nx1] Speech signal
%  fs         : [Hz]         [1x1] sampling frequency
%  smoothOpt  : [binary]     [1x1] flag for using smoothing
%  normOpt    : [string]     [1x1] 'line', 'mean' or 'nonorm'
%                                  for selecting normalisation type
% Outputs
%  CPP        : [s,samples] [Mx2] CPP with time values
%
% Example
%  Please see the HOWTO_glottalsource.m example file. 
%
% References
%  [1] Hillenbrand, J., Houde, R. A. (1996) ``Acoustic correlates of
%      breathy vocal quality: dysphonic voices and continuous
%      speech'', Journal of Speech and Hearing research 39:311-321
%
% This function is part of the Covarep project: 
% http://covarep.github.io/covarep
% 
% Author 
%  John Kane kanejo@tcd.ie
%
% TODO: Use vectorised regression line fitting 
%       Add band-pass filtering variant

function CPP = cpp( x, fs, smoothOpt, normOpt )

if nargin < 3
    smoothOpt = 0;
end
if nargin < 4
    normOpt = 'mean';
end

%% Settings
filterType = 'highpass';
HPfilt_b = [1 - 0.97];
frameLength = round( 0.04 * fs );

if smoothOpt
   frameShift = round( 0.002 * fs);
   timeSmoothLen = 10; 
   quefSmoothLen = 10; 
else 
   frameShift = round( 0.01 * fs);
end

halfLen = round( frameLength / 2 );
xLen = length( x );
frameLen = halfLen * 2 + 1;
NFFT = 2 ^ ( ceil ( log (frameLen) / log(2) ) );
quef = linspace( 0, frameLen / 1000, NFFT );
F0lim = [ 500, 50 ]; % Note that this differs from Hillenbrands
                     % settings of 60-300 Hz, to allow for a
                     % fuller range of potential F0 values
quefLim = fs ./ F0lim;
quefSeq = ( quefLim(1):quefLim(2) )';

time_samples = frameLength+1:frameShift:xLen-frameLength;
N = length(time_samples);
frameStart = time_samples-halfLen;
frameStop = time_samples+halfLen;

%% Apply filtering if requested
if strcmp( filterType, 'highpass' )
   x = filter( HPfilt_b, 1, x );
end

%% Create frame matrix
frameMat = zeros( NFFT, N );
for n=1:N
   frameMat(1:frameLen,n) = x( frameStart(n):frameStop(n) );
end

%% Apply Hann window function
win = hanning( frameLen);
winMat = repmat( win,1,N );
frameMat = frameMat(1:frameLen,:) .* winMat;

%% Compute magnitude spectrum
SpecMat = abs( fft( frameMat ) );
SpecdB = 10 * log10( SpecMat.^2 );

%% Compute log Cepstrum
Ceps = log( abs( fft( SpecdB ) ).^2 );

%% If selected smooth across time and quefrency
if smoothOpt
   timeFilter_b = ones( 1, timeSmoothLen ) / timeSmoothLen;
   quefFilter_b = ones( 1, quefSmoothLen ) / quefSmoothLen;
   Ceps = filter( timeFilter_b, 1, Ceps' )';
   Ceps = filter( quefFilter_b, 1, Ceps );
end

%% Take quefrency range and compute max
CepsLim = Ceps( quefSeq,: );
[CepsMax,maxIdx] = max( CepsLim, [], 1 );

%% Do normalisation
CepsNorm = zeros( 1, N );

if strcmp( normOpt, 'line' )
   for n=1:N
      p = polyfit( quefSeq, CepsLim(:,n), 1 );
      CepsNorm(n) = polyval( p, quefSeq( maxIdx(n) ) );
   end
elseif strcmp( normOpt, 'nonorm' )==0
   CepsMean = mean( CepsLim, 1 );
end
 
CPP = CepsMax - CepsNorm;
CPP = [CPP(:) time_samples(:)];



