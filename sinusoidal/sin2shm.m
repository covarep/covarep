% Build a spectrum-like harmonic model using a given sinusoids set
%
% Input
%  sins     : Nx3 (or Nx4), frequency[Hz], amplitude, phase, (harmnb)
%             Only positive frequencies, with DC.
%  opt      : Options (see below)
%
% Output
%  M        : Spectrum built from sins.
%             (With same format as a DFT result)
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
% TODO
%  Manage holes in the harmonic structure
%

function [M, opt] = sin2shm(sins, opt)

    if nargin<2
        opt.nbh             = Inf; % Number of harmonics to be use in the model
                                   % (excluding DC).
        opt.max_freq        = 8000;% Consider harmonics up to this limit.

        opt.prevpow2        = false;% if true, reduce the harmonic model to the
                                  % biggest possible power of 2 length
        opt.dc_extrap       = 1;  % 0: Use the current value
                                  % 1; Copy the amplitude of the first harm
                                  % 2: Use linear extrapolation on log scale
                                  % both 1 and 2 assume angle(DC)=0
    end
    if nargin==0; M=opt; return; end

    % TODO manage holes in the harmonic structure

    if ~isempty(opt.max_freq) && opt.max_freq>0
        inds = sins(1,:)<=opt.max_freq;
        sins = sins(:,inds);
    end
    sins = sins(:,1:min(size(sins,2), 1+opt.nbh));
    if opt.prevpow2
        % in order to speed up the computation of the min phase
        harmsize = 2^(nextpow2(size(sins,2))-1) + 1;
        sins = sins(:,1:min(size(sins,2), harmsize));
    end

    M = [sins(2,:).'.*exp(1i*sins(3,:).')];

    M(1) = sign(real(M(1)))*abs(M(1));        % set DC phase to 0 or pi

    M(end) = sign(real(M(end)))*abs(M(end));  % set the Nyquist phase to 0 or pi

    M = hspec2spec(M);

    % Extrapolate the model DC
    % Doesn't keep the DC phase !
    % Angle(DC) never makes sens in acoustic measurement
    % For technical reasons, we can therefore assume any constant (here zero).
    if opt.dc_extrap==1;     M(1) = abs(M(2));
    elseif opt.dc_extrap==2; M(1) = db2lin(lin2db(M(2))+(lin2db(M(2))-lin2db(M(3))));   end

return
