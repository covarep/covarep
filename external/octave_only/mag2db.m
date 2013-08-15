% Convert linear amplitude to dB
%
% Copyright (c) 2011 University of Crete - Computer Science Department
%
% License
%  This file is part of libphoni. libphoni is free software: you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General
%  Public License as published by the Free Software Foundation, either version 3
%  of the License, or (at your option) any later version. libphoni is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
%  details. You should have received a copy of the GNU Lesser General Public
%  License along with libphoni. If not, see <http://www.gnu.org/licenses/>.
%
% Author
%  Gilles Degottex <degottex@csd.uoc.gr>
%
% $Id: lin2db.m 134 2012-07-31 15:12:33Z degottex $

function y = mag2db(x)

    if strcmp(class(x),'single');   e = realmin('single');
    else                            e = realmin;            end

    y = 20*log10(abs(x) + e);

return
