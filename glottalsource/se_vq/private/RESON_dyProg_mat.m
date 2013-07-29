% RESON_dyProg_mat.m
% Function to carry out dynamic programming componenent of SE_VQ
%
% Octave compatible
%
% Description
% Function to carry out dynamic programming method described in Ney (1989)
% and used previously in the ESPS GCI detection algorithm. The method
% considers target costs and transition costs which are accumulated in
% order to select the `cheapest' path, considering previous context
%
% Inputs
%  rel_amp         : [samples] [MxNcand] Relative residual amplitudes 
%  GCI_N           : [samples] [MxNcand] matrix containing M by Ncandidate GCI locations 
%  F0mean          : [Hz]      [1x1] Mean fundamental frequency
%  x               : [samples] [Nx1] Speech signal
%  fs              : [Hz]      [1x1] sampling frequency
%  trans_wgt       : [integer] [1x1] Weight of transition cost
%  relAmp_wgt      : [integer] [1x1] Weight of target cost
%
% Outputs
%  GCI_opt         : [samples] [Mx1] Optimal GCI path
%
% Example
%  GCI_opt =
%  RESON_dyProg_mat(GCI_relAmp,GCI_N,F0_mean,x,fs,trans_wgt,relAmp_wgt)
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

function GCI_opt = RESON_dyProg_mat(GCI_relAmp,GCI_N,F0_mean,x,fs,trans_wgt,relAmp_wgt)

%% Initial settings
plots=0;
cost = (GCI_relAmp.*relAmp_wgt)';
ncands=size(GCI_N,1);
nframe=size(GCI_N,2);
GCI_N=GCI_N';
prev=zeros(nframe,ncands);      % traceback pointer
pulseLen = round(fs/F0_mean);
GCI_opt=zeros(1,nframe);

%% Do processing
for n=1:nframe
   
    if n>1
        costm=zeros(ncands,ncands);         % transition cost matrix: rows (previous), cols (current)
        
        for c=1:ncands
            % Transitions TO states in current frame
            start=GCI_N(n,c)-round(pulseLen/2);
            stop=GCI_N(n,c)+round(pulseLen/2);
            if stop > length(x)
                stop=length(x);
            end
            
            if start<1
                start=1;
            end
            pulse_cur = x(start:stop);
            
            for p=1:ncands
                % Transitions FROM states in previous frame
                start=GCI_N(n-1,p)-round(pulseLen/2);
                stop=GCI_N(n-1,p)+round(pulseLen/2);
                if start<1
                    start=1;
                end
                if stop > length(x)
                    stop=length(x);
                end
                if start<1
                    start=1;
                end
                
                pulse_prev = x(start:stop);
                
                if isempty(pulse_cur) ||isnan(pulse_cur(1)) || ...
                        isnan(pulse_prev(1))
                    costm(p,c)=0;
                else
                    if length(pulse_cur)~=length(pulse_prev)
                        cor_cur=0;
                    else
                        cor_cur=corrcoef(pulse_cur(:),pulse_prev(:));
			if length(cor_cur)==1
			   cor_cur=0;
			else cor_cur=cor_cur(2);
			end
                    end
                    costm(p,c) = (1-abs(cor_cur))*trans_wgt; % transition cost
                end
            end
        end
        
        costm=costm+repmat(cost(n-1,1:ncands)',1,ncands);  % add in cumulative costs
        [costi,previ]=min(costm,[],1);
        cost(n,1:ncands)=cost(n,1:ncands)+costi;
        prev(n,1:ncands)=previ;
    end
    
end

%% Do traceback
best=zeros(n,1);
[cbest,best(n)]=min(cost(n,1:ncands));
for i=n:-1:2
     best(i-1)=prev(i,best(i));
end

for n=1:nframe
    GCI_opt(n) = GCI_N(n,best(n));
end

%% Do plots
if plots
    GCI_norm=zeros(nframe,ncands);
    GCI_opt_norm=zeros(nframe,ncands);
    for n=1:nframe
        GCI_norm(n,:) = GCI_N(n,:)-GCI_N(n,1);
        GCI_opt_norm(n) = GCI_opt(n)-GCI_N(n,1);
    end
    subplot(211), plot(x), hold on, stem(GCI_N(:,1),ones(1,length(GCI_N(:,1)))*-.1,'r')
        stem(GCI_opt,ones(1,length(GCI_opt))*-.1,'k')
    subplot(212), plot(GCI_opt,GCI_norm,'rx'), ylim([-20 20]),
        hold on, plot(GCI_opt,GCI_opt_norm)
end