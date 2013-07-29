% Function to estimate glottal closure instants (GCIs) using an extension
% of Thomas Drugman's SEDREAMS algorithm, optimised to better handle
% non-modal phonation types. This version is extended to handle speech with
% wide f0 excursions
%
% Octave compatible
%
% Description
%  Function to extract GCIs using an adapted version of the SEDREAMS
%  algorithm which is optimised for non-modal voice qualities (SE-VQ).
%  Ncand maximum peaks are selected from the LP-residual signal in the interval
%  defined by the mean-based signal. A dynamic programming algorithm is then
%  used to select the optimal path of GCI locations. Then a post-processing
%  method, using the output of a resonator applied to the residual signal, is
%  carried out to remove false positives occurring in creaky speech regions
%  Note that this method is slightly different from the standard SE-VQ
%  algorithm as the mean based signal is calculated using a variable window
%  length. This is set using an f0 contour interpolated over unvoiced
%  regions and heavily smoothed. This is particularly useful for speech
%  involving large f0 excursions (i.e. very expressive speech). Note also
%  that this algorithm has not yet been fully formally tested.
%
% Inputs
%  x               : [samples] [Nx1] Speech signal
%  fs              : [Hz]      [1x1] sampling frequency
%  f0              : [Hz]      [Px1] fundamental frequency
%  VUV             : [binary]  [Px1] Binary voicing decision
%  creak           : [binary]  [Nx1] Creak decision (OPTIONAL)
%
% Outputs
%  GCI             : [s] [Mx1] Glottal closure instants
%
% Example
%  GCI = se_vq_varf0(x,fs);
%
% References
%  [1] Kane, J., Gobl, C., (2013) `Evaluation of glottal closure instant 
%       detection in a range of voice qualities', Speech Communication 
%       55(2), pp. 295-314.
%
% Copyright (c) 2013 Trinity College Dublin
%
% License
%  This code is a part of the Voice Analysis Toolkit with the following
%  licence:
%  The software product (and any modifications thereof) is distributed under 
%  a dual licence: an open source license for individual, non-commercial 
%  purposes, and a commercial license. The opensource licence under which 
%  the product is distributed is GNU GPL v2. For individual users, this 
%  licence suits their use as these are not distributing proprietary 
%  modifications, additions to, or derivatives of the product and don't 
%  require legal protection of a commercial licence. For commercial users, 
%  where open source does not meet their requirements, we offer commercial 
%  licensing of the product. A commercial license permits customers to 
%  modify, add or produce derivative products without the obligation of 
%  making the subsequent code open source. For more information regarding 
%  our commercial licence, please contact john.whelan@tcd.ie
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  John Kane kanejo@tcd.ie
%

function GCI = se_vq_varf0(x,fs,f0,VUV,creak)

%% Settings
F0min=20;
F0max=500;

F0mean=median(f0(f0>F0min&f0<F0max&VUV==1));
F0max=max(medfilt1(f0(VUV==1),13));
if F0mean < 70
    disp('Utterance likely to contain creak')
    F0mean=80;
end

% Interpolate f0 over unvoiced regions and heavily smooth the contour
[TMP,f0_samp] = create_continuous_smooth_f0(f0,VUV,x);
T0mean = fs./f0_samp; % Rough period length for mean-based signal   

winLen = 25; % window length in ms 
winShift = 5; % window shift in ms 
LPC_ord = (fs/1000)+2; % LPC order
Ncand=5; % Number of candidate GCI residual peaks to be considered in the dynamic programming

trans_wgt=1; % Transition cost weight
relAmp_wgt=0.3; % Local cost weight

repNum=2; % Number of iterations for post-processing
removeThresh=0.4; % Threshold for removing false GCIs
search_reg=1.3/1000*fs; % Search region for post processing

%% Calculate LP-residual and extract N maxima per mean-based signal determined intervals
res = lpcresidual(x,winLen/1000*fs,winShift/1000*fs,LPC_ord); % Get LP residual
rep = RCVD_reson_GCI(res,fs,F0mean); % Get resonator output
MBS = get_MBS(x,fs,T0mean); % Extract mean based signal
interval = get_MBS_GCI_intervals(MBS,fs,T0mean,F0max); % Define search intervals
[GCI_N,GCI_relAmp] = search_res_interval_peaks(res,interval,Ncand); % Find residual peaks
GCI = RESON_dyProg_mat(GCI_relAmp',GCI_N',F0mean,x,fs,trans_wgt,relAmp_wgt); % Do dynamic programming

%% Remove false alarms as weak peaks in resonator output
if nargin > 4 && length(creak)==length(x)
    disp('Doing post-processing in detected creaky voice regions')
    GCI = GCI_creak_postproc(GCI,creak,search_reg,rep,removeThresh,repNum);
end

GCI = (GCI-1)/fs;
