% Synthesis of stationary sinusoids using OverLap-Add (OLA)
% 
% Octave compatible
% 
% Description
%  This is a simple OverLap-Add (OLA) synthesis method. Thus, no frequency
%  matching (connections between sinusoids between neighbor frames) is done.
%
% Inputs
%  frames  : N structures containing sinusoidal parameters.
%            For each frame, the sinusoid parameters are in a matrix with format:
%               [4xK] for each column: the frequency [Hz], the linear amplitude,
%               the instantaneous phase [rad] and the harmonic number of each
%               sinusoidal component.
%               The DC is ALWAYS included at the beginning of the matrix.
%
%  [wavlen]: [samples] Length of the waveform
%  [fs]    : [Hz] Sampling frequency of the synthesized waveform
%  [opt]   : Additional options (see code below)
%  
% Outputs
%  syn     : The synthesized waveform
%  fs      : The sampling frequency of the synthesized waveform
%
% Example
%  Please se the HOWTO_sinusoidal example
%  
% Copyright (c) 2012 University of Crete - Computer Science Department
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

function [syn fs opt] = sin_synthesis_ola(frames, wavlen, fs, opt)

    if nargin<4
        % Options
        opt.win_fn     = @hanning; % Bounds have to be strictly >0

        % Options
        opt.syn_dc     = true;

        opt.debug      = 1;

        opt.fmapin{1} = {{'frames', 'wavlen', 'fs'}, 'sin'};
        opt.fmapout{1}  = {{'syn', 'fs'}, 'snd'};
    end
    if nargin==0; syn=opt; return; end

    disp('Sinusoidal Synthesis');

    if opt.debug>0; disp(opt); end

    syn  = zeros(wavlen,1);
    wins = zeros(wavlen,1);

    pb = progressbar(numel(frames));
    for ind=1:numel(frames)

        winlen = frames(ind).winlen;
        win = opt.win_fn(winlen);
        dftlen = frames(ind).dftlen;

        winids = -(winlen-1)/2:(winlen-1)/2;
        ids = round(frames(ind).t*fs)+1 + winids; % Indices of the window in the signal

        if ids(1)<1 || ids(end)>wavlen;
            error('Part of the frame is out of the signal');
            continue;
        end

        % Get the sinusoids parameters used for the synthesis
        sins = frames(ind).sins;
        if ~opt.syn_dc
            sins = sins(:,2:end); % If asked, drop the DC
        end

        % Synthesize the frame content
        y = sin2sig(sins, fs, winlen);
        y = 2*y; % the signal is assumed to be real,
                 % thus "add" the negative frequencies 

        syn(ids) = syn(ids) + y.*win; % Window the frame and add it to the sig
        wins(ids) = wins(ids) + win;  % Sum the windows to normalize after

        if 0 && opt.debug>1 && T(ind)>0.3
            keyboard
        end
        
        pb = progressbar(pb, ind);
    end
    pb = progressbar(pb, numel(frames));

    % Normalize the signal by the sum of the windows
    idx = find(wins>0);
    syn(idx) = syn(idx)./wins(idx);

    if opt.debug>1
        keyboard
    end

return

%===============================================================================
