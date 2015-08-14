% The Maxima Dispersion Quotient (MDQ) quantifies how impulse-like the
% glottal excitation is through wavelet analysis of the Linear Prediction
% (LP) residual
%
% Octave compatible
%
% Description
%  Function to calculate the Maxima Dispersion Quotient (MDQ) from the
%  inputted LP-residual of the speech signal. The method involves carrying 
%  out zero-phase wavelet based filtering on the LP-residual of the input signal. 
%  The dispersion of peaks across the different frequency bands are measured
%  in relation to the glottal closure instant, are averaged and then
%  normalised to the local glottal period.
%
%  Performance appears to be optimal for F0 in the range [50, 200],
%  between [200, 300] it is satisfactory but after 300 Hz the
%  performance deteriorates significantly. One idea would be to
%  measure the strength of periodicity in the higher scales (lower
%  frequencies) to determine whether higher ones should be omitted
%  from the dispersion measurement.
%
% Inputs
%  res      : [samples] [Nx1] Linear Prediction residual of speech signal
%  fs       : [Hz]      [1x1] sampling frequency
%  GCI      : [s] [Mx1] Glottal closure instants 
%
% Outputs
%  m        : [s,samples] [Mx2] Maxima Dispersion Quotient (values aligned to GCI)
%
% Example
%  Please see the HOWTO_glottalsource.m example file.
%
% References
%  [1] Kane, J., Gobl, C., (2013)``Wavelet maxima dispersion for breathy to tense 
%      voice discrimination'', IEEE Trans. Audio Speech & Language
%      Processing, 21(6), pp. 1170-1179.
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

function m = mdq(res,fs,GCI)

GCI = round(GCI*fs)+1;

%% Initial settings
GCI=unique(GCI);

i=2:6; % Analysis scales [2000, 1000, 500, 250, 125]
s_num=length(i);
searchRate=0.2;

m=zeros(1,length(GCI));

%% Do wavelet-based decomposition
[~,y_n] = do_daless_decomp(res,fs,i);

%% Do processing
for n=1:length(GCI)
   
    % Get local Glottal period
    if n==1
        T0=GCI(n+1)-GCI(n);
    else T0=GCI(n)-GCI(n-1);
    end
    
    searchLen_cur=round(searchRate*T0); % Search region as a function of T0
    
    
    if GCI(n)-T0 > 0 && GCI(n)+T0 <= length(res)

        % Get refined frame start and end points
        start_ser=GCI(n)-searchLen_cur;
        finish_ser=GCI(n)+searchLen_cur;
        midpoint=searchLen_cur;
        
        % Current frame
        y_cur=y_n(:,start_ser:finish_ser);

        dist_cur=zeros(1,s_num); % Allocate space for dispersion measure

        % Measure dispersion at each frequency band
        for k=1:s_num
            y_curn = y_cur(k,:);
            [~,maxIdx] = max(y_curn);
            dist_cur(k) = abs(midpoint-maxIdx);
        end  

        % Normalise average disperion to the local glottal period
        m(n)= mean(dist_cur)/T0; 
       
     end
end

GCI(isinf(m))=0; % Remove any infinity values
m(isinf(m))=0; % Remove any infinity values
m=[(GCI(:)-1)/fs m(:)];
