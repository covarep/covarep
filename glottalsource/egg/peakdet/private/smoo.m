function SMOO_dSIG = smoo(dSIG,C_SMOO)
% Smoothing a digitalized signal.
%
% This is the second version of this function (2.0), written in 2019.
% The earlier version (call it version 1) was unsophisticated: a simple moving average (i.e. unweighted moving average), thus:
% 
% % filling first and last values with original values
% % % for i = 1 : C_SMOO
    % % % SMOO_dSIG(i) = dSIG(i);
    % % % SMOO_dSIG(length(dSIG) + 1 - i) = dSIG(length(dSIG) + 1 - i);
% % % end
% %
% % smoothing
% % % for i = 1 + C_SMOO : length(dSIG) - C_SMOO
    % % % SMOO_dSIG (i) = sum(dSIG(i - C_SMOO:i + C_SMOO)) / (2 * C_SMOO + 1);    
% % % end
%
% Version 2 follows the example of <praatdet> (at https://github.com/kirbyj/praatdet#smoothing): it uses a linearly weighted symmetric moving average. A quick tutorial on this topic is available here:
%https://fr.mathworks.com/videos/using-convolution-to-smooth-data-with-a-moving-average-in-matlab-97193.html

% Creating the mask
for i = 1:(C_SMOO + 1)
	mask (i) = i;
	mask (2 * (C_SMOO+1) - i) = i;
end
% Dividing the mask by <C_SMOO>, so the amplitude remains constant
mask = mask / C_SMOO;

% Applying the convolution. The parameter 'same' means that the length remains unaffected.
SMOO_dSIG = conv(dSIG, mask, 'same');

% Decreasing the amplitude of the signal: otherwise the smoothed signal is much higher than the non-smoothed signal, making comparison difficult
SMOO_dSIG = SMOO_dSIG / (2 * C_SMOO);

