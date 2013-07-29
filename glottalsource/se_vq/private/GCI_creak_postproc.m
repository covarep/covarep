% GCI_creak_postproc.m
% Function to do the post-processing step to remove false positive GCIs
% detected in creaky voice regions
%
% Octave compatible
%
% Description
% Function to do the post-processing step to remove false positive GCIs
% detected in creaky voice regions
%
% Inputs
%  GCI             : [samples] [Mx1] Glottal closure instants
%  creak           : [binary]  [Nx1] Creak decision 
%  search_reg      : [integer] [1x1] Search factor
%  rep             : [samples] [Nx1] Resonator output
%  removeThresh    : [integer] [1x1] Threshold for false GCI removal
%  repNum          : [integer] [1x1] number of repititions
%
% Outputs
%  GCI             : [samples] [Mx1] Glottal closure instants
%
% Example
%  GCI = GCI_creak_postproc(GCI,creak,search_reg,rep,removeThresh,repNum)
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

function GCI = GCI_creak_postproc(GCI,creak,search_reg,rep,removeThresh,repNum)

% Separate GCIs detected in creaky voice regions from other regions
creak_GCI=GCI(creak(GCI)==1);
GCI(creak(GCI)==1)=[];
    
for m=1:repNum
    n=2;
    while n < length(creak_GCI)
        cur_rep_max = abs(min(rep(creak_GCI(n)-round(search_reg):creak_GCI(n)+round(search_reg))));

        if mean([abs(rep(creak_GCI(n-1))) abs(rep(creak_GCI(n+1)))])*removeThresh > cur_rep_max
             creak_GCI(n)=NaN;
            n=n+2;
        else n=n+1;
        end
    end
    creak_GCI(isnan(creak_GCI))=[];
end
    
GCI=sort(unique([GCI(:)' creak_GCI(:)']));