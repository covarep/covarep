function [y,zf,u]=randfilt(pb,pa,ny,zi)
%RANDFILT Generate filtered gaussian noise without initial transient
%
% Generates noise by passing gaussian noise through a filter pb/pa
% output signal y is a column vector of length ny
%
% zf and zi are the output and optional input state as defined in filter()
% If zi is not specified, random numbers with the correct covariance will be used.
% output u*u' is the state covariance matrix for filter().

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: randfilt.m,v 1.4 2007/05/04 07:01:39 dmb Exp $
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

% check to see if we must generate zi

if pa(1)~=1
   pb=pb/pa(1); pa=pa/pa(1);
end
if nargin<4 | nargout>2
   lb=length(pb);
   la=length(pa);
   
   k=max(la,lb)-1;
   l=la-1;
   ii=k+1-l:k;
   
   % form controllability matrix
   
   q=zeros(k,k);
   [z,q(:,1)]=filter(pb,pa,1);
   for i=2:k [z,q(:,i)]=filter(pb,pa,0,q(:,i-1)); end
   
   % we generate m through the step-down procedure
   s=randn(k,1);
   if l
      m=zeros(l,l);
      g=pa;
      for i=1:l
         g=(g(1)*g(1:end-1)-g(end)*g(end:-1:2))/sqrt((g(1)-g(end))*(g(1)+g(end)));
         m(i,i:l)=g;
      end
      s(ii)=triu(toeplitz(pa(1:l)))*(m\s(ii));
      if nargout>2
         u=q;
         u(:,ii)=q(:,ii)*triu(toeplitz(pa(1:l)))/m;
      end
   else
      if nargout>2
         u=q;
      end
   end
   if nargin < 4
      if k
         zi=q*s;
      else
         zi=[];
      end
   end
end
if ny
   [y,zf]=filter(pb,pa,randn(ny,1),zi);
else
   zf=zi;
   y=[];
   end