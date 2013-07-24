% Synthesis of a time signal given sinusoids

% sins: matrix of size(3,N) [freq;amp;phase]

function [s t] = sin2sig(sins, fs, winlen, deriv)
    if nargin<4; deriv=0; end

%      t = (0:winlen-1)/fs;
%      t = (0:winlen-1)/fs;
    % The time reference used for the phase in sins is in the window center
    t = (-(winlen-1)/2:(winlen-1)/2)/fs;

    if deriv==0
        A = cos(2*pi*t'*sins(1,:) + ones(length(t),1)*sins(3,:));
        A = ones(length(t),1)*sins(2,:).*A;
        s = sum(A,2);
    elseif deriv==1
        A = -sin(2*pi*t'*sins(1,:) + ones(length(t),1)*sins(3,:));
        A = ones(length(t),1)*(sins(2,:).*(2*pi*partials(1,:))).*A;
        s = sum(A,2);
    end

return
