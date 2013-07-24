function [zeroXingRate,pos] = get_zero_x_rate(x,N,shift)

% Function to measure zero crossing rate to be used as a measure in creaky
% voice detection

%% Calculate Intraframe Periodicity Measure (IFP)
start=1;
finish=start+N-1;
n=1;

while finish <= length(x)
    frame = x(round(start):round(finish)); % select frame
    
    xing_idx=zerocros(frame,'b');
    zeroXingRate(n) = length(xing_idx)/length(frame);
    
    start=start+shift;
    finish=start+N-1;
    pos(n) = mean([start finish]);
    n=n+1;
end

%%%%%%%
% INTERNAL FUNCTIONS
%%%%%%%

 function [t,s]=zerocros(x,m)
%ZEROCROS finds the zeros crossings in a signal [T,S]=(X,M)% find zero
%crossings in a signal

if nargin<2
   m='b';
end
s=x>=0;
k=s(2:end)-s(1:end-1);
if any(m=='p')
   f=find(k>0);
elseif any(m=='n')
   f=find(k<0);
else
   f=find(k~=0);
end
s=x(f+1)-x(f);
t=f-x(f)./s;
if ~nargout
   n=length(x);
   plot(1:n,x,'-',t,zeros(length(t),1),'o');
end

