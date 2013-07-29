% search_res_interval_peaks.m
% Function to search for Ncand peaks in the residual within each search
% interval
%
% Octave compatible
%
% Description
% Function to search for Ncand peaks in the residual within each search
% interval
%
% Inputs
%  res             : [samples] [Nx1] Linear Prediction residual
%  interval        : [samples] [Mx2] Search intervals for GCI estimation
%  Ncand           : [integer] [1x1] Number of candidates
%
% Outputs
%  GCI             : [samples] [MxNcand] Glottal closure instants candidates
%  rel_amp         : [samples] [MxNcand] Relative residual amplitudes  
%
% Example
%  [GCI,rel_amp] = search_res_interval_peaks(res,interval,Ncand);
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

function [GCI,rel_amp] = search_res_interval_peaks(res,interval,Ncand)

% Function to search for Ncand peaks in the residual within each search
% interval
N=size(interval,1);
GCI=zeros(N,Ncand);
rel_amp=zeros(N,Ncand);

for n=1:N
   
    start=round(interval(n,1));
    stop=round(interval(n,2));
    
    if stop <= start
        GCI(n,:)=NaN;
    elseif stop-start < Ncand
        [amp,idx]=max(res(start:stop));
        GCI_cur=GCI_cur+start-1;
        
        GCI(n,:)=GCI_cur;
        rel_amp(n,:)=0;
        
    else [amp,idx]=sort(res(start:stop),'descend');
        GCI_cur=idx(1:Ncand)+start-1;
        GCI(n,:)=GCI_cur(:)';
        rel_amp(n,:)=1-(amp(1:Ncand)/max(amp(1:Ncand)));
    end
end