% Measure of f0 and Oq based on Derivated ElectroGlottoGraphic signal (DEGG)
%
% Octave compatible
% 
% Description
%  Using the time Derivative of an ElectroGlottoGraphic signal (DEGG), this
%  function measures the fundamental frequency (f0) (f0mes),
%  the Open Quotient (Oq) (Oqeggmes), the number of peaks at closure instant
%  (npicferm) and the number of peaks at opening instant (npicouver).
%
% Original French description:
%  DECOM  Mesure de f0 et Oq sur le signal DEGG
%  mesure la fréquence fondamentale f0mes, le quotient d'ouverture Oqeggmes,
%  le nombre de pics à la fermeture npicferm (vérification de l'unicité du pic),
%  le nombre de pics à l'ouverture npicouver, sur un signal DEGG,
%  connaissant les bornes en fréquence f0 =  [f0min,f0max] et la fréquence
%  d'échantillonnage fs. Par défaut, f0min = 80, f0max = 1500.
%  Si f0 ne donne qu'une valeur, elle est directement utilisée pour le
%  calcul.
%  La variable 'method' définit la méthode utilisée pour le
%  calcul de Oq : 'max', le pic de la fonction
%  d'intercorrélation ayant une amplitude maximale est
%  sélectionné, 'premier' le premier des pics de la fonction
%  d'intercorrélation est sélectionné, 'dernier' le dernier
%  des pics est sélectionné.
%
% Inputs
%  s       : The input EGG signal, in VFCA convention (vocal-fold contact area),
%             i.e. EGG increases with glottal contact. The peak at closure is
%             thus positive in the EGG derivative.
%            It is recommended to high-pass the signal (around 50Hz) in order
%            to remove the influence of the larynx global movements.
%  fs      : [Hz] The sampling frequency of the EGG signal.
%  f0      : [1x2 Hz] The search limits of f0, in a row vector.
%            Default value is: f0=[80Hz, 1500Hz].
%            If a single value is used (e.g. f0=120), it replaces the
%            measured f0 and estimate Oq using this f0 value.
% 
%  method    : Select the method to compute Oq. It can be:
%              'max'   : The peak of the intercorrelation with the maximum
%                        amplitude is used (default method).
%              'first' : The first peak of the intercorrelation is used.
%              'last'  : The last peak of the intercorrelation is used.
%
% Outputs
%  Oqeggmes  : Measured Open Quotient (Oq)
%  f0mes     : Measured Fundamental frequency
%  npicferm  : Number of peaks at glottal closing
%  npicouver : Number of peaks at glottal opening
%
% Example
%  See the HOWTO_egg.m example file.
%
% Reference
% [1] N. Henrich, C. d'Alessandro, M. Castellengo, B. Doval (2004) "On the use
%     of the derivative of electroglottographic signals for characterization
%     of nonpathological phonation," Journal of the Acoustical Society of
%     America, vol.115, no.3, p.1321-1332, 2004.
%     Original publication of the code: http://voiceresearch.free.fr/egg
%
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

function [Oqeggmes,f0mes,npicferm,npicouver] = decom(s, fs, f0, method)

if nargin <3 || isempty(f0)
    f0min = 80;
    f0max = 1500;
    f0 = [NaN;NaN];
end
if length(f0)==2
    f0min = f0(1);
    f0max = f0(2);
end
if nargin < 4 || isempty(method)
    method = 'max';
end

% initialisation des tableaux
indn=[];
indp=[];
ep=[];
en=[];
interc=[];
autoc=[];

% valeur par défaut
npicferm = -1;
npicouver = -1;


% separation partie positive / partie negative
indn = find(s < 0);
indp = find(s > 0);
ep = s;
en = - s;
ep(indn) =0;
en(indp) =0;

% longueur des vecteurs s, ep et en
L = length(s);


%%%%%%%%%%%%%%%
% CALCUL DE F0
%%%%%%%%%%%%%%%

if length(f0) == 1
    f0mes = f0;
    % précision de 0.5 Hz
    deltat = 0.5 / f0mes^2;

else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul de l'autocorrelation normalisée (et biaisée)
    % (autocorrélation non circulaire)
    %
    % autre méthode (autocorrélation circulaire)
    %
    % autoc = real(ifft(abs(fft(ep,2*L)).^2)) / L;
    % autoc = autoc(1:L);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    autoc = xcorr(ep,'coeff');
    autoc = autoc(L:2*L-1);

    % détection du pic correspondant à la période
    % seuil de 0.5 (critère voisé / non voisé)
    [maxa,indmaxa]=PicDetect(autoc,0.5);

    if isempty(maxa)
        f0mes=0;
        Oqeggmes=0;
    elseif autoc(max(1,indmaxa(1)-4))>=autoc(indmaxa(1))
        f0mes=0;
        Oqeggmes=0;
    else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % estimation de f0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % estimation grossière
        f0approx=fs/(indmaxa(1)-1);

        % précision de 0.5 Hz
        deltat = 0.5 / f0approx^2;

        % interpolation cubique entre 5 points
        % centrés autour du pic détecté

        % sélection temporelle
        indsel = [max(1,indmaxa(1)-2):indmaxa(1) indmaxa(1)+1 indmaxa(1)+2];
        tsel = (indsel-1)/fs;
        autocsel = autoc(indsel);
        % rééchantillonnage
        ti = min(tsel):deltat:max(tsel);
        % interpolation
        autoci = interp1(tsel,autocsel,ti,'spline');
        % détection du maximum
        [ymax1,indmax1] = max(autoci);

        % estimation de f0 à 0.5 Hz près
        t1 = ti(indmax1);
        f0mes = 1/t1;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vérification des bornes en fréquence
        % de f0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ((f0mes>=f0max) || (f0mes<=f0min))
            f0mes=0;
            Oqeggmes=0;
        end
    end
end



%%%%%%%%%%%%%%%
% CALCUL DE OQ
%%%%%%%%%%%%%%%

if f0mes == 0
    Oqeggmes = 0;

else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % vérification de l'unicité du pic
    % correspondant à l'ouverture
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    npicouver = NombredePics(en,f0mes,fs,0.7,40);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul de l'intercorrelation
    % entre le signal positif (fermeture)
    % et le signal négatif (ouverture)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interc = xcorr(ep,en,'coeff'); % DIVA
    interc = interc(L:2*L-1);

    % détection du (des) pic(s) correspondant
    % seuil à 0.8
    [maxid,indmaxid]=PicDetect(interc,0.75);

    [maxi,indi] = max(maxid);

    if isempty(indi)
        Oqeggmes = 0;
        return;
    end

    if npicouver == 1
        indmaxi = indmaxid(indi);

    else
        switch method
            case 'max'

                % sélection du pic ayant une amplitude maximale
                indmaxi = indmaxid(indi);

            case 'first'
                indmaxi = min(indmaxid);

            case 'last'
                indmaxi = max(indmaxid);
        end
    end


    % détection précise du maximum
    % sélection temporelle
    indsel = [max(1,indmaxi(1)-2):indmaxi(1) indmaxi(1)+1 indmaxi(1)+2];
    tsel = (indsel-1)/fs;
    intercsel = interc(indsel);

    % rééchantillonnage
    ti = min(tsel):deltat:max(tsel);
    % interpolation
    interci = interp1(tsel,intercsel,ti,'spline');
    % détection du maximum
    [ymax2,indmax2] = max(interci);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % valeur de Oq mesurée
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t2 = ti(indmax2);
    Oqeggmes= t2 * f0mes;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % vérification de l'unicité du pic
    % correspondant à la fermeture
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    npicferm = NombredePics(ep,f0mes,fs,0.5,40);


    if Oqeggmes >1
        Oqeggmes = 0;
    end


end



