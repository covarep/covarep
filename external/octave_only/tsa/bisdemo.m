% BISDEMO (script) Shows BISPECTRUM of eeg8s.mat

%       $Id: bisdemo.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1997, 1998,2008 by Alois Schloegl <a.schloegl@ieee.org>
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

load eeg8s.mat;
[BISPEC,BICOV,ACF]=bispec(eeg8s,30);
[s1,s2]=size(BISPEC);
t1=(1:s1)/max(s1);
t2=(1:s2)/max(s2);
subplot(211);
mesh(t1,t2,abs(BISPEC));
title('Bispectrum - mesh plot');

subplot(212);
if exist('OCTAVE_VERSION')>5
        contour(abs(BISPEC),10,t1,t2);
else
        contour(t1,t2,abs(BISPEC),10);
end;        
title('Bispectrum - contour plot');
