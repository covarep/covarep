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

function sins = env_extrap_sins_dcny(sins, fs, method, mf0)
    if nargin<3; method=1; end
    if nargin<4; mf0=[]; end

    if isstruct(sins)
        for fi=1:numel(sins)            
            sins(fi).sins = env_extrap_sins_dcny(sins(fi).sins, fs, method, mf0);
        end

    else
        % Drop any weird frequencies
        while sins(1,end)<0;    sins=sins(:,2:end); end
        while sins(1,end)>fs/2; sins=sins(:,1:end-1); end

        sins = sins(1:2,:);
        if isempty(mf0); mf0=median(diff(sins(1,:)),2); end

        % Fix DC, which is zero in acoustic signals, but complicate any env estimate
        % If there is a DC, just replace it
        if sins(1,1)==0
            a0 = sins(2,2);
    %          a0 = db2lin(ld(sins(2,1))+(ld(sins(2,1))-ld(sins(2,2))));
            sins(1:2,1) = [0;a0];
        elseif sins(1,1)>mf0/4
            % If no DC, add it
            a0 = sins(2,1);
            sins = [[0;a0], sins(1:2,:)];
        end

        % Add extra sinusoids up to Nyquist because the sin model might be
        % band-limited lower than Nyquist.
        if method==0
            % Add according to mean freq distance, up to Nyquist (not included)
            if sins(1,end)~=fs/2
                while sins(1,end)<fs/2-mf0
                    if size(sins,1)==2
                        sins = [sins, [sins(1,end)+mf0; sins(2,end)]];
                    elseif size(sins,1)==4
                        sins = [sins, [sins(1,end)+mf0; sins(2,end); sins(3,end); sins(4,end)+1]];
                    end
                end
            end
        elseif method==1
            % Add points uniformly from last frequency up to Nyquist
            d = fs/2-sins(1,end);
            nbp = ceil(d/mf0);
            afs = [(sins(1,end)+d/nbp:d/nbp:fs/2-0.5*mf0), fs/2];
            sins = [sins, [afs; sins(2,end)*ones(size(afs))]];

%          if ~isempty(find(isnan(sins(1,:)))); keyboard; end

        end
    end

return
