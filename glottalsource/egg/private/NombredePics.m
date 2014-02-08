% Copyright (c) 2006 CNRS, Laboratoire d'Acoustique Musicale (Paris)
%                          and Institut de la Communication Parlée (Grenoble)
%
% License
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Nathalie Henrich <nathalie.henrich@gipsa-lab.fr> (April 4th, 2006)
%

function [npic] = NombredePics(s,f0,Fe,seuil,BW)
%NombredePics	Détection du nombre de pics d'un signal
%
%	[npic] = NombredePics(s,f0,Fe,seuil,BW) repère les pics
%	d'un signal situé au-dessus d'un certain seuil, par
%	intercorrélation avec une fonction peigne d'impulsion 
%	gaussienne de même largeur que le signal s, de fréquence f0
%	et de largeur de bande BW (BW exprimé en %. Plus BW
%	est petit, plus la largeur de bande sera grande
%   BW = 10 -> largeur à mi-hauteur de 8% de la période
%   BW = 20 -> largeur à mi-hauteur de 4% de la période)

%	Auteur : Nathalie Henrich, 04/08/00.
%	$Révision: 


d0 = Fe/f0;
ind = min(find(s == max(s)));  
indpos = fix((ind/d0-fix(ind/d0))*d0 - d0/2);

% génération d'un train d'impulsion gaussienne
% de largeur de bande BW : 
% t = 0:1/Fe:2/f0;
% ygauss = exp(-t.^2 * (pi/2)^2 / (-(1/BW/f0)^2*log(10^(-6/20)))).*cos(2*pi*f0*t);

T=0:1/Fe:length(s)/Fe;
D=indpos/Fe:1/f0:length(s)/Fe;

y=pulstran(T,D,'gauspuls',f0,BW);

% intercorrélation entre le train d'impulsion 
% et le signal s
N = length(s);
interc = xcorr(y,s); % DIVA
interc = interc(N : 2*N-1);

% sélection d'une période sur le signal
% d'intercorrélation
indmax = find(interc(1:min([fix(3*d0/2) length(interc)])) == max(interc(1:min([fix(3*d0/2) length(interc)]))));

if indmax<d0/2
   intercsel = interc(indmax+fix(d0/2):min(N,indmax+fix(3*d0/2)));
else
   intercsel = interc(max(1,indmax-fix(d0/2)):min(N,indmax+fix(d0/2)));
end

% détection des pics du signal
if fix(seuil*3/2) >= length(intercsel)
seuil
figure(3)
plot(intercsel)
waitforbuttonpress
    seuil = length(intercsel);
end

[valpic,indpic]	=	PicDetect(intercsel,seuil);

% nombre de pics détectés
npic = length(indpic);

