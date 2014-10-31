% Get an interpolated value from a time-data vector
%
% Description
%  Linear interpolation of a time-data vector (see below the corresponding
%  input). Works only with sorted times.
%
% Inputs
%  values   : [time,data] vector where time and data are column vectors
%  t        : instant(s) where the interpoaltion is obtained.
%  value_fn : Instead of a linear interpolation, a specific function can be used
%  value_ifn : The inversion function of the previous one.
%
% Outputs
%  value : interpolated value
%  previ : The index of the previous element used for interpolation
%  nexti : The index of the next element used for interpolation
%  g     : The weight used between the previous and the next value
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

function [value previ nexti g] = interp1td(values, t, value_fn, value_ifn)

    if nargin<3
        value_fn=[];
        value_ifn=[];
    end

    if length(t)>1
        value = ones(length(t),size(values,2));
        for ind=1:length(t)
            value(ind,1) = t(ind);
            value(ind,2:end) = interp1td(values, t(ind), value_fn, value_ifn);
        end
    else
        times = values(:,1);

        previ = find(times<=t, 1, 'last');
        if isempty(previ);
            value = values(1,2:end);
            previ = 1;
            nexti = 1;
            g = 1;
        elseif previ==size(values,1)
            value = values(end,2:end);
            nexti = previ;
            g = 0;
        elseif times(previ)==t
            value = values(previ,2:end);
            nexti = previ;
            g = 0;
        else
            nexti=previ+1;

            g = (t-times(previ))./(times(nexti)-times(previ));

            if g<0 || g>1; error('linearity factor unbound, g not in [0;1]'); end
            
            % linear interpolation
            if nargin<3 || isempty(value_fn)
                value = values(previ,2:end)*(1-g) + values(nexti,2:end)*g;
            else
                value = value_fn(values(previ,2:end))*(1-g) + value_fn(values(nexti,2:end))*g;
                value = value_ifn(value);
            end
        end
    end

return
