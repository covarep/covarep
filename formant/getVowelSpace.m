% This script estimates the vowel space ratio of a speaker with respect to
% his/her reference population.
%
%
% Description
% Based on the tracked formants F1 and F2 for the voiced regions of speech 
% we compute the vowel space for each recorded subject individually [1]. 
% We define the vowel space as the frequency region covered by the triangle 
% in the two dimensional frequency space spanned by F1 and F2 for the vowels 
% /i/ (as in heed), /a/ (as in hod), and /u/ (as in who'd) following [2].
% These three vowels represent the vowels with the most extreme positions 
% of the tongue and are therefore located in the extremes of this triangularly 
% shaped two-dimensional frequency space [3]. The method is developed
% following the recommendations in [4]. Reference frequencies of american
% english vowels are found in [5].
%
%
% [ratio, centroids] = getVowelSpace(formant_data, gender, mode, scaling, plot_formants)
%
% Inputs
%  formant_dat      : [samples] [Nx2] input formant observations F1 and F2
%  gender           : [nominal] [1x1] female = 0; male = 1; child = 2
%  mode             : [code]    [string] 'triangle' (default); 'polygon'
%                               (add /ae/ (as in had) to the vowel space)
%  scaling          : [boolean] [1x1] bark scaling on/off (default);
%                               mel_scaling optional in code.
%  plot_formants    : [boolean] [1x1] visualize vowel space on/off
%                               (default)
%
% Outputs
%  ratio            : [1x1] vowel space ratio
%  centroids        : [1x12] vector containing prototype vowel locations
%
% References
% [1] Scherer, S., Morency, L.-P., Gratch, J., and Pestian, J., EDUCED VOWEL SPACE
% IS A ROBUST INDICATOR OF PSYCHOLOGICAL DISTRESS: A CROSS-CORPUS ANALYSIS,
% Proceedings of ICASSP 2015.
% [2] H.-M. Liu, F.-M. Tsao, and P. K. Kuhl, ?The effect of reduced vowel
% working space on speech intelligibility in mandarin-speaking young
% adults with cerebral palsy,? Journal of the Acoustical Society of America,
% vol. 117, no. 6, pp. 3879?3889, 2005.
% [3] B. Lindblom, ?Explaining phonetic variation: A sketch of the h&h
% theory,? Speech Production and Speech Modeling, pp. 403?439, 1990.
% [4] S. Sandoval, V. Berisha, R. L. Utianski, J. M. Liss, and A. Spanias,
% ?Automatic assessment of vowel space area,? Journal of the Acoustical
% Society of America, vol. 134, no. 5, pp. 477?483, 2013.
% [5] J. Hillenbrand, L. A. Getty, M. J. Clark, and K. Wheeler, ?Acoustic
% characteristics of american english vowels,? Journal of the Acoustical
% Society of America, vol. 97, no. 5, pp. 3099?3111, 1995.
%
% Copyright (c) 2014 University of Southern California, Institute for
% Creative Technologies
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
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
% 
% Author 
%  Stefan Scherer scherer@ict.usc.edu
%
function [ratio, centroids] = getVowelSpace(formant_data, gender, mode, scaling, plot_formants)

%% initialization
mel_scaling = @(x)(2595.*log10(ones(size(x))+x./700));
bark_scaling = @(x)(13.*atan(0.00076.*x)+3.5.*atan((x./7500).^2));
min_observations = 1000;

switch nargin
    case 1
        gender = 1;
        warning('getVowelSpace.m : no gender specified. assumes male gender.');
        plot_formants = 0;
        scaling = false;
        mode = 'triangle';
    case 2
        if(~ismember(gender, [0,1,2]))
            error('getVowelSpace.m : wrong use of gender variable; permitted values: 1 .. male; 0 .. female; 2 .. child');
        end
        plot_formants = 0;
        scaling = false;
        mode = 'triangle';
    case 3
        if(~ismember(gender, [0,1,2]))
            error('getVowelSpace.m : wrong use of gender variable; permitted values: 1 .. male; 0 .. female; 2 .. child');
        end
        if(~ismember(mode, {'triangle','polygon'}))
            error('getVowelSpace.m : wrong use of gender variable; permitted values: 1 .. male; 0 .. female; 2 .. child');
        end
        scaling = false;
        plot_formants = 0;
    case 4
        if(~ismember(gender, [0,1,2]))
            error('getVowelSpace.m : wrong use of gender variable; permitted values: 1 .. male; 0 .. female; 2 .. child');
        end
        plot_formants = 0;
    case 5
        if(~ismember(gender, [0,1,2]))
            error('getVowelSpace.m : wrong use of gender variable; permitted values: 1 .. male; 0 .. female; 2 .. child');
        end
    otherwise
        error('getVowelSpace.m : wrong use of getVowelSpace.');
end

if(size(formant_data,2) < 2)
    error('getVowelSpace.m : wrong use of getVowelSpace. you must provide at least the first two formants!');
end

%% Formants for vowels
% list of formants; based on [4]
% this list is quite flexible.. you can adjust that to your dialect/language if you like!
men_f1 = [342, 427, 476, 580, 588, 768, 652, 497, 469, 378, 623, 474];
men_f2 = [2322, 2034, 2089, 1799, 1952, 1333, 997, 910, 1122, 997, 1200, 1379];

fem_f1 = [437, 483, 536, 731, 669, 936, 781, 555, 519, 459, 753, 523];
fem_f2 = [2761, 2365, 2530, 2058, 2349, 1551, 1136, 1035, 1225, 1105, 1426, 1588];

child_f1 = [452, 511, 564, 749, 717, 1002, 803, 597, 568, 494, 794, 586];
child_f2 = [3081, 2552, 2656, 2267, 2501, 1688, 1210, 1137, 1490, 1345, 1546, 1719];


if strcmp(mode, 'polygon')
    target_vowels = [1, 5, 6, 10];
else
    target_vowels = [1, 6, 10];
end

%% initialize data for processing
if(gender == 1)
    init_formants = [men_f1; men_f2]';
    formant_data = formant_data(formant_data(:,1) < 900 & formant_data(:,2) > 800 & formant_data(:,1) > 250 & formant_data(:,2) < 2450, :); % remove some outliers of formants
elseif(gender == 0)
    init_formants = [fem_f1; fem_f2]';
    formant_data = formant_data(formant_data(:,1) < 1000 & formant_data(:,2) > 1000 & formant_data(:,1) > 350 & formant_data(:,2) < 2800, :); % remove some outliers of formants
elseif(gender == 2)
    init_formants = [child_f1; child_f2]';
    formant_data = formant_data(formant_data(:,1) < 1000 & formant_data(:,2) > 1000 & formant_data(:,1) > 350 & formant_data(:,2) < 2800, :); % remove some outliers of formants
end

if(size(formant_data,1) < min_observations)
    warning('getVowelSpace.m : More formant observations are required for this method to work robustly.');
    ratio = 0;
    centroids = [];
    return
end

%% bark scaling
if scaling
    formant_data = bark_scaling(formant_data);
    init_formants = bark_scaling(init_formants);
    % optional:
    % formant_data = mel_scaling(formant_data);
    % init_formants = mel_scaling(init_formants);
end

%% vector quantization
warning off;
[~,centroids] = kmeans(formant_data(:,1:2),[],'start',init_formants, 'emptyaction', 'singleton');
warning on;


%% temp data for clustering
% distance measure can be changed!
pairwise_distances = pdist2(init_formants, centroids, 'seuclidean');
[~, idx] = min(pairwise_distances(target_vowels, :), [], 2);

%% triangle ratio
convHullPoints=centroids(idx,:);
vowelSpacePoints = init_formants(target_vowels,:);
[convOrder , cur_space] = convhull(convHullPoints);
[spaceOrder , vowel_space] = convhull(vowelSpacePoints);
ratio = cur_space./vowel_space;

%% plot vowel space
if(plot_formants)
    figure
    scatter(formant_data(:,1),formant_data(:,2), '.', 'LineWidth', 0.05);
    hold on;
    plot(centroids(:,1), centroids(:,2), 'ro')
    plot(centroids(idx,1), centroids(idx,2), 'go')
    xlabel('Formant 1');
    ylabel('Formant 2');
    plot(init_formants(:,1), init_formants(:,2), 'kx')
    plot(init_formants(target_vowels,1), init_formants(target_vowels,2), 'mx')
    plot(convHullPoints(convOrder,1),convHullPoints(convOrder,2),'g-')
    plot(vowelSpacePoints(spaceOrder,1),vowelSpacePoints(spaceOrder,2),'m-')
    set(gca,'XDir','reverse');
    set(gca,'YDir','reverse');
    set(gca,'XAxisLocation','top','YAxisLocation','right');
    title(sprintf('vowel space ratio %2.3f (#obs = %d)', ratio, size(formant_data, 1)));
    drawnow;
    hold off;
end
