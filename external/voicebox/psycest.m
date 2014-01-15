function [xx,ii,m,v]=psycest(iq,x,r,xp)
% Estimate multiple psychometric functions
%
% Usage: [xx,ii,m,v]=psycest(-n,p,q,xp) % initialize n models
%        [xx,ii,m,v]=psycest(i,x,r)     % supply a trial result to psycest
%                    psycest(i)         % plot pdf of model i
%              [p,q]=psychest(0)        % output model parameters (or print them if no outputs)
%
% Inputs (see usage examples for argument positions):
%         -n        minus the number of models
%          p(:,n)   parameters for each model
%                      1  thresh [0.75]
%                      2  miss [0.04]
%                      3  guess [0.1]
%                      4  SNR min [-20]
%                      5  SNR max [20]
%                      6  Slope min [0]
%                      7  slope max [0.5]
%          q(:)     parameters common to all models (vector or struct)
%                      1  q.nx  number of SNR values in pdf [40]
%                      2  q.ns  number of slope values in pdf [21]
%                      3  q.nh  number of probe SNR values to evaluate [30]
%                      4  q.cs  weighting of slope relative to threshold in cost function [1]
%                      5  q.dh  minimum step size in dB for probe SNRs [0.2]
%                      6  q.sl  min slope at threshold [0.02]
%                      7  q.kp  number of std deviations of the pdfs to keep [4]
%                      8  q.hg  amount to grow expected gains in ni trials [1.3]
%                      9  q.cf  cost function: 1=variance, 2=entropy [2]
%                     10  q.pm  psychometric model: 1=logistic, 2=cumulative gaussian [1]
%                     11  q.lg  use log slope in pdf: 0=no, 1=yes [1]
%          xp{n}(:) list of available probe SNR values
%          i        model probed
%          x        probe SNR value used
%          r        response: 0=wrong, 1=right.
%
% Outputs:
%          xx       recommended probe SNR
%          ii       recommended model to probe next
%          m(2,n,3) estimated srt and slope of all models
%                   m(:,:,1:3) are respectively the mean, mode (MAP) and marginal mode estimates
%          v(3,n)   estimated covariance matrix entries:
%                   [var(srt) cov(srt,slope) var(slope)]'
%
% Algorithm parameters:
%
% The algorithm estimates the psychometric function for a number of models simultaneously.
% A psychometric function specifies the probability of correct recognition as a function of
% SNR and is completely specified by two parameters: the SNR (in dB) and the slope (in 1/dB)
% at a specified target recognition probability (e.g. 0.5 or 0.75). The p(:,n) parameters specify
% some of the details of the psychometric function and can be different for each model:
%   p(1,n) gives the target recognition probability
%   p(2,n) gives the probability of incorrect recognition at very good SNRs (the "miss" probability).
%          If this value is made too low, a single unwarrented recognition error will have a large
%          effect of the estimated parameters and greatly lengthen the time to convergence. The default
%          value is 4%.
%   p(3,n) gives the probability of correct recognition at very poor SNRs (the "guess" probabiity).
%          This should be set to 1/N where N is the number of possible responses to a stimulus.
% p(4:5,n) These give the initial min and max SNR at the target recognition probability. They will
%           be adjusted by the algorithm if necessary but wrong values will delay convergence.
% p(6:7,n) These given the initial min and max slope (in probability per dB) at the target
%          recognition probability.
% The remaining parameters are shared between all the models and control how the algorithm operates.
%   q(1:2) These parameters determine the sampling density of the joint pdf of SNR and Slope.
%   q(3)   This parameter specifies how many potential values of probe SNR to evaluate in order to
%          determine which is likely to give the most improvement in the parameter estimates.
%   q(4)   This parameter gives the relative weight of SNR and Slope in the cost function. Increasing
%          its value will improve the slope estimate (or log(slope) estimate) at the expense of the
%          SNR estimate. Actually its value is not all that critical.
%   q(5)   At each iteration psycest evaluates several potential probe SNR levels to see which is
%          likely to give the most improvement to the parameter estimates. This parameter fixes the
%          minimum spacing between these potential probe values. It only has an effect when the variances
%          of the parameter estimates are very small.
%   q(6)   To prevent the routine getting lost, this parameter specifies the smallest reasonable value
%          of the Slope parameter at the target recognition probability. Under normal circumstances, it
%          has no effect.
%   q(7)   For each model, the routine maintains a sampled version of the joint pdf of the SNR and slope.
%          The sampling grid is reclculated after each iteration to encompass this number of standard
%          deviations in each dimension.
%   q(8)   At each iteration, psycest advises which model to probe next on the basis of which
%          gives the greatest expected reduction in cost function. To ensure that all models
%          are periodically sampled, this expected reduction of each unprobed model is multiplied
%          by this parameter at each iteration. Thus a model will eventually be probed even if
%          its expected cost factor improvement is small.
%   q(9)   This determines whether the recommended probe SNR to use at the next iteration is chosen to
%          minimize the expected variance or the entropy of the estimated SNR and slope distributions.
%          My experience is that entropy (the default) gives faster convergence.
%  q(10)   This selects whether the underlying model is a logistic function (1) or a cumulative
%          gaussian (2). The choice makes little difference unless the "miss" or "guess" probabilities
%          are very very small.
%  q(11)   This selects whether the internal pdf samples "slope" or "log(slope)". It is recommended
%          that you stick to the default of log(slope) since this results in more symmetrical
%          distributions.

%      Copyright (C) Mike Brookes 2009-2010
%      Version: $Id: psycest.m 1933 2012-06-07 09:41:30Z dmb $
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bugs/Suggestions:
% (1) allow lookahead by 2 tests rather than 1
% (2)use a structure for input parameters
% (3)should only resample the pdfs when necessary
% (4)add a forgetting factor for old measurements
% (8) use quadratic interpolation to find cost function minimum
% (10) optionally output the whole pdf + its axis values
% (11) optionally output all model probe values
% (13) Should perhaps estimate and return the mean and slope compensated for the guess rate
%      i.e. if you want p, you actually test at guess+ p*(1-guess-miss)/(1-miss)
% (17) remember probe snrs, models and results and recalculate the entire
%      array when changing precision
% (18) could force probe SNRs to be a multiple of minimum step size
% (19) use a non-uniform prior e.e. 1/(1+x^2)
% (20) possibly use different parametrization (e.g. 1/slope or a second SRT threshold)
% (22) save inputs so that pdfs can be recalculated when rescaling
% (24) Parametrize slope as x=(100s^2-1)/20s to make the distribution more uniform
%      inverse is s=(sqrt(x^2+1)-x)/10; alternatively use log(slope)
% (25) optionally have no prior (i.e. maximum likelihood)
% (26) Check why the estimated variances are too small
% (27) Adapt range of potential probes if optimum is at the limit
% (28) Resample the pdf on the basis of the marginal distributions
% (29) Select the probe range on the basis of marginal distributions
% (30) Resample the pdfs by a factor of 2 always
% (31) use uneven samples in pdf concentrated in the middle
% (32) Use a parametric mixture including one-sided constants for the pdf
% (33) Selection of subset of prescribed probe SNRs is not optimum
% (34) use quadratic interpolation to select the probe SNR unless using fixed values
% (35) expand the pdf grid based on the effective number of samples (like particle filters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent variables
% sizes: nx=SNR values in pdf
%        ns=slope values in pdf (1/std dev of cumulative gaussian)
%        nh=potential test SNRs
%        ni=simultaneous models being estimated
% wq(ns*nx,ni) = prob of model; each column is a vectorized ns*nx matrix
% nr(1,10)      = parameter values: [SNR-pdf slope-pdf SNR-probe ...]
% pr(7,ni)     = input model parameters
% qr(4,ni)     = derived model parameters
% xq(nr(1),ni) = SNR values in pdf
% sq(nr(2),ni) = slope values in pdf
% mq(2,ni)     = estimated srt and slope of all models
% vq(3,ni)     = estimated covariance matrix entries:[var(srt) cov(srt,slope) var(slope)]'
% xn(1,ni)     = next probe value to use
% hn(1,ni)     = expected decrease in cost function after next probe

persistent wq xq sq nr pr qr mq vq xn hn hfact xz

% algorithm parametes that could be programmable
pfloor=0.001;    % pdf floor relative to uniform distribution
trynsig=2; % number of standard deviations to explore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iq<0  % initialization
    ni=-iq;                                                 % number of models
    pr=repmat([0.75 0.04 0.1 -20 20 0 0.5]',1,ni);          % default parameters
    if nargin>1
        if size(x,2)>1
            if size(x,2)>ni
                error('initialization parameter argument has too many columns');
            end
            pr(1:size(x,1),1:size(x,2))=x;
        else
            pr(1:size(x,1),:)=repmat(x,1,ni);
        end
    end
    nr=[40 21 30 1 0.2 0.02 4 1.3 2 1 1]';                      % default parameter values
    nrf={'nx','ns','nh','cs','dh','sl','kp','hg','cf','pm','lg'};     % parameter field names
    numnr=length(nr);
    if nargin>2
        if isstruct(r)
            fnn=fieldnames(r);
            for i=1:length(fnn)
                mk=strcmp(fnn{i},nrf);
                if any(mk)
                    nr(mk)=r.(fnn{i});
                end
            end
        else
            nr(1:min(numnr,numel(r)))=r(:);
        end
        nr(1:3)=round(nr(1:3));
    end
    pr(6,:)=max(pr(6,:),nr(6));                         % low limit of slope in prob/dB
    nsxq=nr(1)*nr(2);
    xq=(0:nr(1)-1)'*(pr(5,:)-pr(4,:))/(nr(1)-1)+repmat(pr(4,:),nr(1),1);
    if nr(11)  % if log slope
        sq=(1-nr(2):0)'*(log(pr(7,:))-log(max(pr(6,:),nr(6))))/(nr(2)-1)+repmat(log(pr(7,:)),nr(2),1);
    else
        sq=(0:nr(2)-1)'*(pr(7,:)-pr(6,:))/(nr(2)-1)+repmat(pr(6,:),nr(2),1);
    end
    wq=repmat(1/nsxq,nsxq,ni);                          % initialize to uniform pdf
    qr=zeros(5,ni);
    qr(1,:)=1-pr(2,:)-pr(3,:);                          % prob range covered by cumulative gaussian
    qr(2,:)=(pr(1,:)-pr(3,:))./qr(1,:);                 % cumulative gaussian probability at threshold
    switch nr(10)
        case 1
            qr(3,:)=log(qr(2,:)./(1-qr(2,:)));                      % x position of target in std measure
            qr(4,:)=qr(1,:).*qr(2,:).*(1-qr(2,:));                  % slope*"stddev" at threshold
        case 2
            qr(3,:)=norminv(qr(2,:));                           % x position of target in std measure
            qr(4,:)=qr(1,:).*normpdf(qr(3,:));                  % slope*stddev at threshold
        otherwise
            error('Unrecognised psychometric model selection');
    end
    mq=repmat([mean(xq,1); mean(sq,1)],[1,1,3]);        % initial means
    vq=[var(xq,1,1); zeros(1,ni); var(sq,1,1)];         % initial variances
    %     hn=[1-nr(4) 0 nr(4)]*vq;  % artificially high value of cost function ensures all models are probed early on
    hn=repmat(Inf,1,ni);    % very high initial cost function
    hfact=nr(8)^(1/ni);     % growth factor to ensure that no model is neglected for too long
    xn=mq(1,:);             % start at mean value
    if nargin>=4 && ~isempty(xp)
        if iscell(xp)
            xz=xp;
        else
            xz=repmat(num2cell(xp(:)',2),ni,1);        % else replicate for each model
        end
        for i=1:ni
            [dummy,j]=min(abs(xz{i}-mq(1,i)));      % select the closest available probe to the mean
            xn(i)=xz{i}(j);
        end
    else
        xz=cell(ni,1);          % make empty cells
    end
elseif iq>0 && nargin==3        % update pdfs with a new probe result
    nxq=nr(1);
    nsq=nr(2);
    nxh=nr(3);
    nsxq=nxq*nsq;
    %     thresh=pr(1,iq);            % target threshold
    guess=pr(3,iq);             % guess rate (1/choices)
    pscale=qr(1,iq);            % prob range left after subtracting miss and guess probs
    xtstd=qr(3,iq);             % x position of target in std measure
    sqi=sq(:,iq);               % slope (or log slope) values in PDF
    if nr(11)  % if log slope
        sqis=exp(sqi)/qr(4,iq);     % inverse std dev of gaussian (proportional to slope)
    else
        sqis=sqi/qr(4,iq);          % inverse std dev of gaussian (proportional to slope)
    end
    xqi=xq(:,iq);
    wqi=wq(:,iq);
    %
    % update probabilities with the previous test result
    % If r==1, we multiply by a horizontally reflected version of the psychometric function that equals
    % the threshold (e.g. 0.75) at the probe SNR, x, for all slopes.
    % If r==0, we multiply by 1- this reflected version which therefore equals (1-thresh) at x.
    %
    r0=r==0;
    switch nr(10)
        case 1
            wqi=wqi.*(r0+(1-2*r0)*(guess+pscale*((1+exp(reshape(sqis*xqi'-xtstd,nsxq,1)-repmat(sqis,nxq,1)*x)).^(-1)))); %  P(l | r,x)
        case 2
            wqi=wqi.*(r0+(1-2*r0)*(guess+pscale*normcdf(repmat(sqis,nxq,1)*x-reshape(sqis*xqi'-xtstd,nsxq,1)))); %  P(l | r,x)
        otherwise
            error('Unrecognised psychometric model selection');
    end
    wqi=wqi./sum(wqi,1);                        % normalize
    %     wq(:,iq)=wqi;               % save updated probabilities (not needed if we always interpolate them)

    % **** this section copied to graph plotting ****
    % Calculate mean and covariance and entropy

    wqsx=reshape(wqi,nsq,nxq);
    px=sum(wqsx,1);                             % p(x0)
    ps=sum(wqsx,2);                             % p(s0)
    xe=px*xqi;                                  % E(x0)
    se=ps'*sqi;                                 % E(s0)

    [pxpk,xm]=max(px);                          % marginal mode in x
    if xm>1 && xm<nxq                           % use quadratic interpolation if possible
        [dum,xm2]=quadpeak(px(xm-1:xm+1)');
        xm=xm+xm2-2;
    end
    xm=(2-xm)*xqi(1)+(xm-1)*xqi(2);             % marginal mode(x0)

    [pspk,sm]=max(ps);
    if sm>1 && sm<nsq                           % use quadratic interpolation if possible
        [dum,sm2]=quadpeak(ps(sm-1:sm+1));
        sm=sm+sm2-2;
    end
    sm=(2-sm)*sqi(1)+(sm-1)*sqi(2);             % marginal mode(s0)

    [wqpk,j]=max(wqi);
    i=1+floor((j-1)/nsq);
    j=j-nsq*(i-1);
    if i>1 && i<nxq && j>1 && j<nsq             % use quadratic interpolation if possible
        [dum,ji]=quadpeak(wqsx(j-1:j+1,i-1:i+1));
        i=i+ji(2)-2;
        j=j+ji(1)-2;
    end
    xj=(2-i)*xqi(1)+(i-1)*xqi(2);               % joint mode  x
    sj=(2-j)*sqi(1)+(j-1)*sqi(2);               % joint mode: s
    xv=px*(xqi.^2)-xe^2;                        % Var(x0)
    sv=ps'*(sqi.^2)-se^2;                       % Var(s0)
    sxv=wqi'*(repmat(sqi,nxq,1).*reshape(repmat(xqi',nsq,1),nsxq,1))-xe*se; % Cov(s0*x0)

    % ************   end of copied section

    mq(:,iq,1)=[xe; se];                        % save means
    mq(:,iq,2)=[xj; sj];                        % save means
    mq(:,iq,3)=[xm; sm];                        % save means
    vq(:,iq)=[xv; sxv; sv];                     % save covariance matrix
    xh=(px*log(px)')*(xqi(1)-xqi(2));           % h(x0)
    sh=(ps'*log(ps))*(sqi(1)-sqi(2));          	% h(s0)

    % now estimate the next probe SNR $$$$

    if ~numel(xz{iq})                               % if no list of probe SNRs was specified
        ytry=exp(0.25i*pi*(0:7))';   % points around the circle
        ytry=[real(ytry) imag(ytry)];
        [vtry,dtry]=eig([xv sxv; sxv sv]);  % eigendecomposition of covariance matrix
        tryxs=repmat([xe,se],8,1)+trynsig*ytry*sqrt(dtry)*vtry';
        pmin=0.05;                             % target probe success probability
        if nr(11)
            tryxs(:,2)=qr(4,iq)*exp(-tryxs(:,2));        % convert log(slope) to std dev
        else
            tryxs(:,2)=qr(4,iq)*tryxs(:,2).^(-1);             % convert slope to std dev
        end
        switch nr(10)
            case 1 % logistic
                qmax=max(tryxs(:,1)+(log((1-pmin)/pmin)-xtstd)*tryxs(:,2));
                qmin=min(tryxs(:,1)+(log(pmin/(1-pmin))-xtstd)*tryxs(:,2));
            case 2 % cumulative gaussian
                qmax=max(tryxs(:,1)+(norminv(1-pmin)-xtstd)*tryxs(:,2));
                qmin=min(tryxs(:,1)+(norminv(pmin)-xtstd)*tryxs(:,2));
        end
        dxt=max(nr(5),(qmax-qmin)/nxh);      % minimum step size of nr(5) [0.2 dB]
        xt=(qmin+qmax)/2+((1:nxh)-(1+nxh)/2)*dxt;
    else                                            % if a specific list of probe SNRs exists
        xzi=xz{iq};                                 % xzi is the list of available probe SNRs
        if numel(xzi)<=nxh                          % use all available probe SNRs if there are not too many
            xt=xzi;
        else
            [xt,ixt]=min(abs(xzi-xe));                  % find the closest one to xe ** not necessarily optimum ***
            ixt=max(1,min(1+numel(xzi)-nxh,ixt-floor((1+nxh)/2))); % arrange symmetrically around xt
            xt=xzi(ixt:min(ixt+nxh-1,numel(xzi)));
        end
    end
    nxhp=length(xt);  % xt are the potential probe SNRs
    switch nr(10)
        case 1
            prt=guess+pscale*((1+exp(repmat(reshape(sqis*xqi'-xtstd,nsxq,1),1,nxhp)-repmat(sqis,nxq,1)*xt)).^(-1)); %  P(r | l,x)
        case 2
            prt=guess+pscale*normcdf(repmat(sqis,nxq,1)*xt-repmat(reshape(sqis*xqi'-xtstd,nsxq,1),1,nxhp)); %  P(r | l,x)
    end
    wqt=repmat(wqi,1,nxhp);
    pl1=prt.*wqt;                       % posterior prob given success = p(l | x,r=1) unnormalized
    pl0=(wqt-pl1);                      % posterior prob given failure = p(l | x,r=0) unnormalized
    prh=sum(pl1,1);                     % p(r | x)=Sum{P(r | l,x)*P(l)} [note wqt is normalized] (row vector)
    pl1=pl1./repmat(prh,nsxq,1);        % normalize
    pl0=pl0./repmat(1-prh,nsxq,1);      % posterior prob given failure = p(l | x,r=0) normalized
    px1=squeeze(sum(reshape(pl1,nsq,nxq,[]),1));    % p(x0 | x,r=1)
    px0=squeeze(sum(reshape(pl0,nsq,nxq,[]),1));    % p(x0 | x,r=0)
    ps1=squeeze(sum(reshape(pl1,nsq,nxq,[]),2));    % p(s0 | x,r=1)
    ps0=squeeze(sum(reshape(pl0,nsq,nxq,[]),2));    % p(s0 | x,r=0)
    xet1=xqi'*px1;                                  % E(x0 | x,r=1)
    xvt1=(xqi.^2)'*px1-xet1.^2;                     % Var(x0 | x,r=1)
    xet0=xqi'*px0;                                  % E(x0 | x,r=0)
    xvt0=(xqi.^2)'*px0-xet0.^2;                     % Var(x0 | x,r=0)
    xvt=xvt1.*prh+xvt0.*(1-prh);                    % E(Var(x0 | x ))
    set1=sqi'*ps1;                                  % E(s0 | x,r=1)
    svt1=(sqi.^2)'*ps1-set1.^2;                     % Var(s0 | x,r=1)
    set0=sqi'*ps0;                                  % E(s0 | x,r=0)
    svt0=(sqi.^2)'*ps0-set0.^2;                     % Var(s0 | x,r=0)
    svt=svt1.*prh+svt0.*(1-prh);                    % E(Var(s0 | x ))
    xht1=sum(log(px1).*px1,1);                      % -H(x0 | x,r=1)
    xht0=sum(log(px0).*px0,1);                      % -H(x0 | x,r=0)
    xht=(xht1.*prh+xht0.*(1-prh))*(xqi(1)-xqi(2));	% h(x0 | x)
    sht1=sum(log(ps1).*ps1,1);                      % -H(s0 | x,r=1)
    sht0=sum(log(ps0).*ps0,1);                      % -H(s0 | x,r=0)
    sht=(sht1.*prh+sht0.*(1-prh))*(sqi(1)-sqi(2));	% h(s0 | x)
    switch nr(9)
        case 1
            hx=(xvt + nr(4)*svt)/(1+nr(4));           	% cost function for each possible test SNR
            [hxmin,ix]=min(hx);                         % find the minimum of cost function
            hn(iq)=(xv + nr(4)*sv)/(1+nr(4))-hxmin;    	% expected decrease in cost function
        case 2
            hx=(xht + nr(4)*sht)/(1+nr(4));            	% cost function for each possible test SNR
            [hxmin,ix]=min(hx);                        	% find the minimum of cost function
            hn(iq)=(xh + nr(4)*sh)/(1+nr(4))-hxmin;    	% expected decrease in cost function
        otherwise
            error('Unrecognised cost function option');
    end
    xn(iq)=xt(ix);                              % next probe value for this model
    %         fprintf('Probe range: %.3g %.3g; choose %.3g\n',xt(1),xt(end),xt(ix));


    % rescale the pdfs if necessary

    ssd=sqrt(sv);
    xsd=sqrt(xv);
    if nr(11)
        sq2=linspace(max(log(nr(6)),se-nr(7)*ssd),se+nr(7)*ssd,nsq)';
    else
        sq2=linspace(max(nr(6),se-nr(7)*ssd),se+nr(7)*ssd,nsq)';
    end
    xq2=linspace(xe-nr(7)*xsd,xe+nr(7)*xsd,nxq)';   % new x axis values
    % do linear interpolation in the x direction
    wqi=reshape(wqi,nsq,nxq);           % turn into a matrix for easy interpolation
    xqf=(xq2-xq(1,iq))/(xq(2,iq)-xq(1,iq));
    xqj=ceil(xqf);
    xqf=xqj-xqf;
    xqg=1-xqf;
    xqf((xqj<=0) | (xqj>nxq))=0;
    xqg((xqj<0) | (xqj>=nxq))=0;
    wq2=wqi(:,min(max(xqj,1),nxq)).*repmat(xqf',nsq,1)+wqi(:,min(max(xqj+1,1),nxq)).*repmat(xqg',nsq,1);
    % do linear interpolation in the s direction
    sqf=(sq2-sq(1,iq))/(sq(2,iq)-sq(1,iq));
    sqj=ceil(sqf);
    sqf=sqj-sqf;
    sqg=1-sqf;
    sqf((sqj<=0) | (sqj>nsq))=0;
    sqg((sqj<0) | (sqj>=nsq))=0;
    wq2=wq2(min(max(sqj,1),nsq),:).*repmat(sqf,1,nxq)+wq2(min(max(sqj+1,1),nsq),:).*repmat(sqg,1,nxq);
    % now normalize and apply a floor
    wq2=wq2(:);     % turn back into a vector
    wq2=wq2/sum(wq2,1);  % normalize
    wq2=max(wq2,pfloor/nsxq); %impose a floor
    wq2=wq2/sum(wq2,1);  % normalize again
    sq(:,iq)=sq2;
    xq(:,iq)=xq2;
    wq(:,iq)=wq2;
elseif iq==0        % output model parameters
    xx=pr;
    ii=nr;
    if ~nargout
        pdesc={'Threshold','Miss prob','Guess Prob','Min SNR','Max SNR','Min Slope','Max Slope'};
        fprintf('\n*********\nModel-specific Parameters\n');
        for i=1:7
            fprintf('%4d) %s: ',i,pdesc{i});
            if size(pr,2)>1
                fprintf('%.5g, ',pr(i,1:size(pr,2)-1));
            end
            fprintf('%.5g\n',pr(i,size(pr,2)));
        end
        qdesc={'nx  SNR values in PDF', ...
            'ns  Slope values in PDF', ...
            'nh  Probe SNR values to evaluate', ...
            'cs  Weighting of slope relative to threshold in cost function', ...
            'dh  Min step size in dB for probe SNRs', ...
            'sl  Min slope at threshold', ...
            'kp  Std deviations of the pdfs to keep', ...
            'hg  Amount to grow expected gains in ni trials', ...
            'cf  Cost function', ...
            'pm  Psychometric model'};
        qoptf=[9 10]; % fields with options
        qoptval={'variance','entropy'; ...
            'logistic','cumulative gaussian'};
        fprintf('\nShared Parameters\n');
        for i=1:10
            fprintf('%4d) %s: ',i,qdesc{i});
            j=find(qoptf==i,1);
            if numel(j)
                fprintf('%d=%s\n',nr(i),qoptval{j,nr(i)});
            else
                fprintf('%.5g\n',nr(i));
            end
        end
    end
end

% now select the appropriate model to probe next

if iq~=0
    [hnmin,ii]=max(hn);         % chose model with the biggest expected decrease
    hn=hn*hfact;         % increase values to ensure they all get a chance
    xx=xn(ii);
    m=mq;
    v=vq;
    if nr(11)
        m(2,:,:)=exp(m(2,:,:));  % convert to real slope
        v(2,:)=v(2,:).*m(2,:,1); % correct the covariance
        v(3,:)=v(3,:).*m(2,:,1).^2; % and the slope variance
    end
end

if ~nargout && iq>0
    sqi=sq(:,iq);
    nsq=length(sqi);
    xqi=xq(:,iq);
    nxq=length(xqi);
    wqi=wq(:,iq);
    nsxq=nxq*nsq;

    % **** this section copied from main processing ****
    % Calculate mean and covariance and entropy

    wqsx=reshape(wqi,nsq,nxq);
    px=sum(wqsx,1); % p(x0)
    ps=sum(wqsx,2); % p(s0)
    xe=px*xqi;                            % E(x0)
    se=ps'*sqi;                            % E(s0)
    [pxpk,xm]=max(px);
    if xm>1 && xm<nxq  % use quadratic interpolation if possible
        [dum,xm2]=quadpeak(px(xm-1:xm+1)');
        xm=xm+xm2-2;
    end
    xm=(2-xm)*xqi(1)+(xm-1)*xqi(2);          % marginal mode(x0)
    [pspk,sm]=max(ps);
    if sm>1 && sm<nsq  % use quadratic interpolation if possible
        [dum,sm2]=quadpeak(ps(sm-1:sm+1));
        sm=sm+sm2-2;
    end
    sm=(2-sm)*sqi(1)+(sm-1)*sqi(2);          % marginal mode(s0)
    [wqpk,j]=max(wqi);
    i=1+floor((j-1)/nsq);
    j=j-nsq*(i-1);
    if i>1 && i<nxq && j>1 && j<nsq % use quadratic interpolation if possible
        [dum,ji]=quadpeak(wqsx(j-1:j+1,i-1:i+1));
        i=i+ji(2)-2;
        j=j+ji(1)-2;
    end
    xj=(2-i)*xqi(1)+(i-1)*xqi(2);    % joint mode  x
    sj=(2-j)*sqi(1)+(j-1)*sqi(2);   % joint mode: s
    xv=px*(xqi.^2)-xe^2;               % Var(x0)
    sv=ps'*(sqi.^2)-se^2;               % Var(s0)
    sxv=wqi'*(repmat(sqi,nxq,1).*reshape(repmat(xqi',nsq,1),nsxq,1))-xe*se; % Cov(s0*x0)

    % ************   end of copied section

    % draw final 2-D distribution
    nxq=nr(1);
    nsq=nr(2);
    axm=[3 -2; 2 -1];
    xqi2=[axm*xqi(1:2);xqi];
    sqi2=[axm*sqi(1:2);sqi];
    imagesc(xqi2,sqi2,[zeros(1,2) px/pxpk; zeros(1,nxq+2);ps/pspk zeros(nsq,1) wqsx/wqpk]);
    hold on
    plot(xe,se,'wo',xm,sm,'w^',xj,sj,'w*');
    plot(xqi2(2),se,'wo',xqi2(2),sm,'w^',xqi2(2),sj,'w*');
    plot(xe,sqi2(2),'wo',xm,sqi2(2),'w^',xj,sqi2(2),'w*');
    t=linspace(0,2*pi,200);
    xcir=cos(t);
    scir=sin(t);
    vcir=sqrt((sv*xcir.^2+xv*scir.^2-2*sxv*xcir.*scir)/(xv*sv-sxv^2));
    plot(xe+xcir./vcir,se+scir./vcir,'w-');
    hold off
    axis 'xy';
    %     colorbar;
    %     cblabel('Relative Probability');
    xlabel(sprintf('SNR @ %d%%SRT (dB)',round(pr(1)*100)));
    if nr(11)
        ylabel('Log psychometric Slope at threshold (prob/dB)');
    else
        ylabel('Psychometric Slope at threshold (prob/dB)');
    end
    title('Joint pdf: o mean, * mode, \Delta marginal mode');
end
