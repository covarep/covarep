% get_MBS_GCI_intervals.m
% Function to detect intervals to be used for searching the LP-residual to
% determine GCI locations using the method described in Drugman et al
% (2012).
%
% Octave compatible
%
% Description
% Function to detect intervals to be used for searching the LP-residual to
% determine GCI locations using the method described in Drugman et al
% (2012).
%
% Inputs
%  MBS             : [samples] [Nx1] Mean-based signal
%  fs              : [Hz]      [1x1] sampling frequency
%  T0mean          : [samples] [1x1] Glottal period
%  F0max           : [Hz]      [1x1] Maximum fundamental frequency
%
% Outputs
%  interval        : [samples] [Mx2] Search intervals for GCI estimation
%
% Example
%  interval = get_MBS_GCI_intervals(MBS,fs,T0mean,F0max)
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
% $Id <info set by the versioning system> $

function interval = get_MBS_GCI_intervals(MBS,fs,T0mean,F0max)

% Function to detect intervals to be used for searching the LP-residual to
% determine GCI locations using the method described in Drugman et al
% (2012). 

%% Initial settings
if nargin < 4
    F0max=500;
end
F0max=F0max*2;
T0max=round(fs/F0max);
[idx,TMP]=v_findpeaks(MBS*-1,'minpeakdistance',T0max); % Find locations of negative peaks
N=length(idx);
search_rate=0.28;
search_left_rate=0.01;
interval=zeros(N,2);

%% Do processing
for n=1:N
   
    if length(T0mean)>1
        start=idx(n)-round(T0mean(idx(n))*search_left_rate);
        stop=idx(n)+round(T0mean(idx(n))*search_rate);
    else start=idx(n)-round(T0mean*search_left_rate);
        stop=idx(n)+round(T0mean*search_rate);
    end
    
    if start < 1
        start=1;
    end
    
    % Check start and end points of detected intervals
    if stop > length(MBS) && start < length(MBS)
        stop=length(MBS);
    elseif stop > length(MBS) && start >= length(MBS)
        break
    end
    
    interval(n,1)=start;
    interval(n,2)=stop;
end
    