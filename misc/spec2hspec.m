% Drop the negative frequencies of a spectrum

% Gilles

function hS = spec2hspec(S)

    S = S(:).';

    if mod(length(S),2)==0; hS=S(1:end/2+1);
    else;                   fhS=S(1:(end-1)/2+1);   end
    
return
