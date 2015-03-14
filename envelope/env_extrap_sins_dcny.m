% Fill DC and Nyquist frequency gaps with fake components.
%
% Empty frequency intervals cause most envelope estimates to degenerate,
% at best case around the interval, at worst through all frequencies.
% DC is zero in acoustic signals, and sinusoidal components are not always 
% estimated up to Nyquist frequencies. Thus, such an empty interval always 
% appear at DC and often before Nyquist.
%
% Copyright (c) 2012 University of Crete
%
% Authors
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function sins = env_extrap_sins_dcny(sins, fs)

    if isstruct(sins)
        for fi=1:numel(sins)            
            sins(fi).sins = env_extrap_sins_dcny(sins(fi).sins, fs);
        end

    else
        % Fix DC, which is zero in acoustic signals, but complicate any env estimate
        % If there is a DC, just replace it
        if sins(1,1)==0
            a0 = sins(2,2);
    %          a0 = db2lin(ld(sins(2,1))+(ld(sins(2,1))-ld(sins(2,2))));
            sins(1:2,1) = [0;a0];
        end

        % Add extra sinusoids up to Nyquist because the sin model might be
        % band-limited lower than Nyquist.
        if sins(1,end)~=fs/2
            mfd = median(diff(sins(1,:)));
            while sins(1,end)<fs/2-mfd
                if size(sins,1)==2
                    sins = [sins, [sins(1,end)+mfd; sins(2,end)]];
                elseif size(sins,1)==4
                    sins = [sins, [sins(1,end)+mfd; sins(2,end); sins(3,end); sins(4,end)+1]];
                end
            end
        end
    end

return
