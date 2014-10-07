function SMOO_dSIG = smoo(dSIG,C_SMOO)
% Smoothing a digitalized signal.

% filling first and last values with original values
for i = 1 : C_SMOO
    SMOO_dSIG(i) = dSIG(i);
    SMOO_dSIG(length(dSIG) + 1 - i) = dSIG(length(dSIG) + 1 - i);
end

% smoothing
for i = 1 + C_SMOO : length(dSIG) - C_SMOO
    SMOO_dSIG (i) = sum(dSIG(i - C_SMOO:i + C_SMOO)) / (2 * C_SMOO + 1);    
end