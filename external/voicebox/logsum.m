function y=logsum(x,d,k)
%LOGSUM logsum(x,d)=log(sum(exp(x).*k,d))
%  d gives dimension to sum along

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: logsum.m,v 1.8 2010/05/10 15:07:34 dmb Exp $
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

if nargin==1 || ~numel(d)
    d=[find(size(x)-1) 1];
    d=d(1);
end
n=size(x,d);
if n<=1,
    if nargin<3
        y=x;
    else
        y=x+log(k);
    end
    return;
end
s=size(x);
p=[d:ndims(x) 1:d-1];
z=reshape(permute(x,p),n,prod(s)/n);
q=max(z,[],1);              % we subtract y from each row to avoid dynamic range problems
a=(q==Inf)|(q==-Inf);       % check for infinities
if nargin<3
    y=q+log(sum(exp(z-q(ones(n,1),:)),1));
else
    y=q+log(sum(exp(z-q(ones(n,1),:)).*reshape(permute(k,p),n,prod(s)/n),1));
end
y(a)=q(a);                  % correct any column whose max is +-Inf
s(d)=1;
y=ipermute(reshape(y,s(p)),p);

