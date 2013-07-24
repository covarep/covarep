function l=gmmlpdf(x,m,v,w)
%GMMPDF calculated the pdf of a mixture of gaussians p=(x,m,v,w)
%
% Inputs: n data values, k mixtures, p parameters
%
%     X(n,p)   Input data vectors, one per row.
%     M(k,p)   mixture means, one row per mixture.
%     V(k,p)   mixture variances, one row per mixture (singlton dimensions will be replicated as required)
%              or else V(p,p,k) for full mixture covariance matrixes           
%     W(k,1)   mixture weights, one per mixture. The weights will be normalized by their sum. [default: all equal]
%
% Outputs: (Note that M, V and W are omitted if L==0)
%
%     L(n,1)   log PDF values

%  Bugs/Suggestions
%     (1) Sort out full covariance maatrices
%     (2) Improve plotting
%     (3) Add an extra arument for plotting control

%      Copyright (C) Mike Brookes 2000-2006
%      Version: $Id: gmmlpdf.m,v 1.2 2007/05/04 07:01:38 dmb Exp $
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
[n,p]=size(x);
l=[];           % in case n=0
if nargout>0 || n>0
    x2=x.^2;            % need x^2 for variance calculation
    k=size(m,1);        % number of mixtures
    if size(m,2)~=p
        error('x and m must have the same number of columns');
    end
    if nargin<4
        if nargin<3
            v=1;
        end
        w=ones(k,1);
    end
    w=w/sum(w);         % normalize the weights
    sv=size(v);
    if length(sv)>2 || k==1 && p>1 && sv(1)==p     % full covariance matrices
        error('full covariance matrices not yet implemented');
    else                            % diagonal (or constant) covariance matrices
        if sv(1)==1
            v=v(ones(k,1),:);
        end
        if sv(2)==1
            v=v(:,ones(1,p));
        end
        
        % If data size is large then do calculations in chunks
        
        memsize=voicebox('memsize'); 
        nb=min(n,max(1,floor(memsize/(8*p*k))));    % chunk size for testing data points
        nl=ceil(n/nb);                  % number of chunks
        
        im=repmat(1:k,1,nb); im=im(:);
        
        lpx=zeros(n,1);             % log probability of each data point
        wk=ones(k,1);
        wnb=ones(1,nb);
        vi=v.^(-1);                 % calculate quantities that depend on the variances
        vm=sqrt(prod(vi,2)).*w;
        vi=-0.5*vi;
        
        % first do partial chunk
        
        jx=n-(nl-1)*nb;                % size of first chunk
        ii=1:jx;
        kk=repmat(ii,k,1);
        km=repmat(1:k,1,jx);
        py=reshape(sum((x(kk(:),:)-m(km(:),:)).^2.*vi(km(:),:),2),k,jx);
        mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
        px=exp(py-mx(wk,:)).*vm(:,ones(1,jx));  % find normalized probability of each mixture for each datapoint
        lpx(ii)=log(sum(px,1))+mx;
        ix=jx+1;
        
        for il=2:nl
            jx=jx+nb;        % increment upper limit
            ii=ix:jx;
            kk=repmat(ii,k,1);
            py=reshape(sum((x(kk(:),:)-m(im,:)).^2.*vi(im,:),2),k,nb);
            mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
            px=exp(py-mx(wk,:)).*vm(:,wnb);  % find normalized probability of each mixture for each datapoint
            lpx(ii)=log(sum(px,1))+mx;
            ix=jx+1;
        end
        l=lpx-0.5*p*log(2*pi);   % log of total probability of each data point
    end
end
if nargout==0                        % attempt to plot the result
    if p==1                            % one dimensional data          
        plot(x,l);
    end
end