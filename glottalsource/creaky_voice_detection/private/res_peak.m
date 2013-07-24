function [peak_inter,peak_prom,peak_t,rep] = res_peak(x,fs,F0mean,res,Es)

% Function to generate residual peak prominence contour which makes up the
% second component of the creak detection algorithm.


% Different settings according to the speakers baseline pitch
if F0mean>100
    maxCreakF0=90;
elseif F0mean < 100 && F0mean>=85
    maxCreakF0=65;
elseif F0mean<85
    maxCreakF0=55;
else maxCreakF0=80;
end

% Set window length based on maximum possible creaky F0
winLen=round(fs/maxCreakF0)*2; 

% Resonator settings
Phi=2*pi*1*F0mean/fs;
Rho=0.8;
rep=filter([1 0 0],[1 -2*Rho*cos(Phi) Rho^2],res);

% Measure residual peak prominence
[peak_prom,peak_t] = get_res_peak_prom(rep,fs,winLen,x,Es);

% Interpolate
if length(peak_prom)>1
    peak_inter=interp1(peak_t,peak_prom,1:length(x));
else peak_inter=zeros(1,length(x));
end