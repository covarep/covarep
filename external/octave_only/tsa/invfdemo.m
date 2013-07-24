% invfdemo	demonstrates Inverse Filtering

%	$Id: invfdemo.m 5090 2008-06-05 08:12:04Z schloegl $%	Copyright (C) 1997-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>
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

load eeg8s.mat;	    % load signal
ly=length(eeg8s);
Fs=128;
%a=earpyw(eeg8s',11); % Calculates AR(11) parameters with Yule-Walker method
a=lattice(eeg8s',11); % Calculates AR(11) parameters with Yule-Walker method
		    % The AR parameters are the weight taps of IIR Filter	
isig=filter([1 -a],1,eeg8s);   % Inverse filtering


subplot(221);
plot((1:ly)/Fs,eeg8s);
title('Signal');
xlabel('t [sec]')

subplot(223);
plot((1:ly)/Fs,isig);
xlabel('t [sec]')
title('Inverse filtered process');

subplot(222);
H=abs(fft(eeg8s,128)/ly).^2;
plot(1:Fs,H);
%plot([H mean(H)*ones(Fs,1)]);
ylabel('S(f)')
xlabel('f [Hz]')
title('Spectrum of original signal');

subplot(224);
H=abs(fft(isig,128)/ly).^2;
plot(1:Fs,H);
%plot([H mean(H)*ones(Fs,1)]);
ylabel('S(f)')
xlabel('f [Hz]')
title('Spectrum of inverse filtered signal');
