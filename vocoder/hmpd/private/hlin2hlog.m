% Log-harmonic scale
%
% Inputs
%  Hb    : Below this harmonic limit, the scale is linear
%          (similar to the mel scale which is linear below 1000Hz)
%          (e.g. 12)
%          Based on observation of the LF model, the asymptotic behavior of the
%          spectrum starts around the 12th harmonic (for the most tense voice)
%  Hmax  : The maxmimum number of harmonic considered during synthesis
%          (e.g. 256)
%  order : The reduced number of phase coefficients (e.g. 24)
%
% Outputs
%  hsl   : The log-harmonic scale
%
% Copyright (c) 2013 University of Crete - Computer Science Department(UOC-CSD)/ 
%                    Foundation for Research and Technology-Hellas - Institute
%                    of Computer Science (FORTH-ICS)
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
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function hsl = hlin2hlog(Hb, Hmax, order)

    if 0
        % This one doesn't work so well because when the order is very small
        % (e.g. 12), the last compressed coefficient represent values
        % only up to a few kHz. Thus, random phase values are used above and
        % the voice sounds constantly noisy.
        hsl = 1:Hmax;
        hsl(Hb+1:end) = Hb*(1+log((Hb+1:Hmax)/Hb));
    else

        % This one correct the previous problem using:
        %   The last compressed coef (i.e. order) corresponds to a given Hmax

        % Use linear scale below Hb (higher coefs will be overwritten)
        hsl = 1:Hmax;

        % TODO DROP CASTELJAU
        % Build a Bezier curve to start with a linear scale and finish smoothly
        % at (Hmax,order)
        p(1,:) = [Hb, Hb];
        p(2,:) = [order, order];
        p(3,:) = [Hmax, order];
        t = 0:0.01:1;
        [X,Y,p_bez] = CASTELJAU(0,1,p,t);
        hsl(Hb+1:end) = interp1(p_bez(:,1), p_bez(:,2), (Hb+1):Hmax);

%      af=abs(frq);
%      mel = sign(frq).*log(1+af/700)*k;
%      mr=(700+af)/k;
    end

    if 0
        plot(hs, hsl, 'k');
        hold on;
        keyboard
    end

return


function [X,Y,val] = CASTELJAU(a,b,p,y)

% function val = CASTELJAU(a,b,p,y)
%
% INPUT:  a   Linke Intervallgrenze
%         b   Rechte Intervallgrenze
%         p   Stützstellen (nx2-Matrix)
%         y   Auswertungspunkte (Spaltenvektor)
%
% OUTPUT: val   Werte des Bezierpolynoms an y (mx2-Matrix)
%
% Date:   2007-11-05
% Author: Jonas Ballani

n = size(p,1);
m = length(y);
T = zeros(n,n);
val = zeros(m,2);
X(:,1) = p(:,1);
Y(:,1) = p(:,2);

for j = 1:m
    for i = 2:n
        X(i:n,i) = (b-y(j))/(b-a)*X(i-1:n-1,i-1) + (y(j)-a)/(b-a)*X(i:n,i-1);
        Y(i:n,i) = (b-y(j))/(b-a)*Y(i-1:n-1,i-1) + (y(j)-a)/(b-a)*Y(i:n,i-1);
    end
    val(j,1) = X(n,n);
    val(j,2) = Y(n,n);
end
  
return
