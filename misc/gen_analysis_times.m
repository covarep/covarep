% Generate analysis time instants, optionaly using a given f0 curve
%
% Octave compatible
% 
% Inputs
%  wav            : Wavform vector
%  fs             : Sampling frequency
%  tmargin        : Margin at start and end
%  usesampletimes : Time instants are thoses of samples.
%                   Time instants can't be between two sample times,
%                   such as window centers correspond to samples times
%  method
%     0: Regular instants.
%        varargin{1}: step size [s]
%     1: Keep f0sin time instants and fill with interpolated values.
%        varargin{1}: An f0 curve
%        varargin{2}: number of instants per period
%     2: Generate new instants from the given f0 curve.
%        varargin{1}: An f0 curve
%        varargin{2}: number of instants per period
%        Skip nan, inf and zero values
%
% Outputs
%  f0sout : f0 feature with the generated times
%
% Copyright (c) 2011 University of Crete - Computer Science Department (UOC-CSD)
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
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%

function at = gen_analysis_times(wav, fs, tmargin, usesampletimes, method, varargin)

    if nargin<5;    usesampletimes = true; end
    if nargin<6;    method = 3; end

    wavdur = (length(wav)-1)/fs;

    if method==0
        step = varargin{1};
        disp('    Generate regular analysis time');
        if usesampletimes
            at = round(tmargin*fs)/fs:round(step*fs)/fs:round((wavdur-tmargin)*fs)/fs;
        else
            at = tmargin:step:(wavdur-tmargin);
        end
        at = at';

    elseif method==1
        f0sin = varargin{1};
        nbperper = varargin{2};
        if nbperper>1
            disp(['    Add extra ' num2str(nbperper) ' analysis instants between each f0 value']);
            f0sin = f0sin(find(~isnan(f0sin(:,2)) & ~isinf(f0sin(:,2)) & f0sin(:,2)~=0),:);
            at = [];
            for n=1:size(f0sin,1)-1
                for m=0:nbperper-1
                    t = f0sin(n,1) + (m/nbperper)*(f0sin(n+1,1)-f0sin(n,1));
                    at = [at; t];
                end
            end
        else
            at = f0sin(:,1);
        end

    elseif method==2
        f0sin = varargin{1};
        nbperper = varargin{2};
        disp(['    Adapt analysis instants to f0 curve with ' num2str(nbperper) ' analysis intants per period']);
        f0sin = f0sin(find(~isnan(f0sin(:,2)) & ~isinf(f0sin(:,2)) & f0sin(:,2)~=0),:);
        at = [];
        nt = tmargin;
        while nt < wavdur - tmargin
            at = [at; nt];
            nt = nt + (1/nbperper)/interp1td(f0sin, nt);
        end

    end
    idx = find(at>tmargin & at<wavdur-tmargin);
    at = at(idx);

    if usesampletimes; at = round(fs*at)/fs; end

return
