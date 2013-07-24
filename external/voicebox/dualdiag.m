function [a,d]=dualdiag(w,b)
%DUALDIAG Simultaneous diagonalisation of two hermitian matrices [A,D]=(W,B)
% Given two hermitian matrices W and B with W positive definite, this routine
% calculates A such that A'*W*A=I and A'*B*A=diag(D). The D will be in descending order.
%
% Suppose we have several N-dimensional data row-vectors arising from each of C different classes of data.
% for each class, c, we can form the mean data vector m(c) and the within-class covariance matrix W(c)
% We can then form the between class covariance matrix B by taking the covariance of the mean vectors m(1), m(2), ...
% and also the averaged within-class covariance matrix W by averaging W(1), W(2), ...
% If we then take A=dualdiag(W,B) and postmultiply all our original data vectors by A, we obtain new
% data vectors for which the average within-class covariance matrix is the identity and for which
% the first few components contain most of the information that is useful in discriminating between classes.

% An alternative algorithm that is 20% faster but slightly less accurate is:
% n=size(w,1);
% [v,l]=eig(w\b);
% [s,i]=sort(-diag(l));
% s=-s;
% d=l(i*(n+1)-n);
% q=sqrt(diag(v'*w*v))'.^(-1);
% a=v(:,i).*q(ones(n,1),i);


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: dualdiag.m,v 1.4 2007/05/04 07:01:38 dmb Exp $
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

[y,l]=eig(w+w');
z=y*diag(sqrt(diag(l*0.5)).^(-1));
[u,s,v]=svd(z'*(b+b')*z);
d=diag(s)*0.5;
a=z*u;
