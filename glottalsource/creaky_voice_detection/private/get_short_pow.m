function [pow,pow_std,pow_std_inter] = get_short_pow(x,fs)

% Function to calculate the very short term power contour as is done in
% Ishi et al. (2008)

% Get very short term power contour
veryShort_len = 4*(fs/1000); % 4ms frame length for "very short-term" analysis
veryShort_shift = 2*(fs/1000); % 2ms shift for "very short-term" analysis
veryShort_powCont = zeros(1,ceil((length(x)-veryShort_len)/veryShort_shift));
t_pow = zeros(1,ceil((length(x)-veryShort_len)/veryShort_shift));
start=1;
finish=start+veryShort_len-1;

n=1;
while finish <= length(x)
    veryShort_powCont(n) = mean(x(start:finish).^2);
    t_pow(n)=mean([start finish]);
    start = start + veryShort_shift;
    finish=start+veryShort_len-1;
    n=n+1;
end

pow = 20*log10(veryShort_powCont);
inf_idx=find(isinf(pow));
pow2=pow;
pow2(inf_idx)=Inf;
pow(inf_idx)=min(pow2);

pow_std=zeros(1,length(pow));
std_len=16;

for n=std_len+1:length(pow)-std_len
    
    pow_std(n) = std(pow(n-std_len:n+std_len));
end
pow_std=medfilt1(pow_std,13);


pow_std_inter=interp1(linspace(1,length(x),length(pow)),pow_std,1:length(x));