function [y,k,y0]=maxfilt(x,f,n,d,x0)
%MAXFILT find max of an exponentially weighted sliding window  [Y,K,Y0]=(X,F,N,D,X0)
%
% Inputs:  X  Vector or matrix of input data
%          F  exponential forgetting factor in the range 0 (very forgetful) to 1 (no forgetting)
%             F=exp(-1/T)  gives a time constant of T samples
%          N  Length of sliding window [default = size(X,D)]
%          D  Dimension for work along [default = first non-singleton dimension]
%         X0  Initial values placed in front of the X data
%
% Outputs: Y  Output matrix - same size as X
%          K  Index array: Y=X(K). (Note that these value may be <=0 if input X0 is present)
%         Y0  Last N-1 values (used to initialize a subsequent call to maxfilt())
%
% This routine calaulates y(p)=max(f^r*x(p-r), r=0:n-1) where x(r)=-inf for r<1
% y=x(k) on output

% Example: find all peaks in x that are not exceeded within +-w samples
% w=4;m=100;x=rand(m,1);[y,k]=maxfilt(x,1,2*w+1);p=find(((1:m)-k)==w);plot(1:m,x,'-',p-w,x(p-w),'+')

%      Copyright (C) Mike Brookes 2003
%      Version: $Id: maxfilt.m,v 1.7 2008/06/02 08:11:45 dmb Exp $
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

if nargin>4                 % initial values specified
    y=shiftdim(cat(d,x0,x),d-1);    % concatenate x0 and x along d
    nx0=size(x0,d);
elseif nargin<4
    [y,d]=shiftdim(x);      % no dimension specified
    nx0=0;
else
    y=shiftdim(x,d-1);      % dimension specified
    nx0=0;
end
s=size(y);
l=s(1);
if nargin<3
    n=l;
end
if nargout>2
    ny0=min(l,n-1);
    sy0=s;
    sy0(1)=ny0;
    y0=reshape(y(1+l-ny0:end,:),sy0);
    if nargin<4
        y0=shiftdim(y0,-d);
    else
        y0=shiftdim(y0,ndims(x)-d+1);
    end
end
n=min(n,l);         % no point in having n>l
k=repmat((1:l)',[1 s(2:end)]);
j=1;
j2=1;
while j>0
    g=f^j;
    m=find(y(j+1:l,:)<=g*y(1:l-j,:));
    m=m+j*fix((m-1)/(l-j));
    y(m+j)=g*y(m);
    k(m+j)=k(m);
    j2=j2+j;
    j=min(j2,n-j2);
end
if nargout==0
    if nargin<3
        x=shiftdim(x);
    else
        x=shiftdim(x,d-1);
    end
    ss=prod(s(2:end));
    plot(1:l,reshape(y(nx0+1:end,:),l,ss),'r',1:l,reshape(x,l,ss),'b');
else
    if nargin<4
        y=shiftdim(y,-d);
        k=shiftdim(k,-d);
    elseif nargin>4
        s(1)=s(1)-nx0;
        y=shiftdim(reshape(y(nx0+1:end,:),s),ndims(x)-d+1);
        k=shiftdim(reshape(k(nx0+1:end,:),s),ndims(x)-d+1)-nx0;
    else
        y=shiftdim(y,ndims(x)-d+1);
        k=shiftdim(k,ndims(x)-d+1);
    end
end