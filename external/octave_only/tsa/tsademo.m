% TSADEMO	demonstrates INVEST1 on EEG data

%       $Id: tsademo.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1998-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>
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


if exist('OCTAVE_VERSION')>5;
    load -force eeg8s.mat 
elseif 1
    load eeg8s.mat 
else
    [FileName, PathName]=uigetfile('eeg8s.mat','load demo data EEG8S.MAT');
    load([PathName FileName],'eeg8s');
end;
s = eeg8s';
Pmax=100;
[AutoCov,AutoCorr,ARPMX,E,CRITERIA,MOPS]=invest1(s,Pmax,'s');

if size(ARPMX,2)==2*Pmax,
	%invest1(eeg8s,30,'s');
        AR=ARPMX(:,1:Pmax);
        RC=ARPMX(:,Pmax+1:2*Pmax);
else
	AR=ARPMX(:,Pmax/2*(Pmax-1)+(1:Pmax));
	RC=ARPMX(:,(1:Pmax).*(2:Pmax+1)/2);
end;
