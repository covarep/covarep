function FLAG = flag_implicit_skip_nan(i)
% FLAG_IMPLICIT_SKIP_NAN sets and gets default mode for handling NaNs
%	1 skips NaN's (the default mode if no mode is set)
% 	0 NaNs are propagated; input NaN's give NaN's at the output
% 
% FLAG = flag_implicit_skip_nan()
% 	gets current mode
%
% flag_implicit_skip_nan(FLAG)    % sets mode 
%
% prevFLAG = flag_implicit_skip_nan(nextFLAG)
%	gets previous set FLAG and sets FLAG for the future
% flag_implicit_skip_nan(prevFLAG)
%	resets FLAG to previous mode
%
% It is used in: 
%	SUMSKIPNAN, MEDIAN, QUANTILES, TRIMEAN
% and affects many other functions like: 
%	CENTER, KURTOSIS, MAD, MEAN, MOMENT, RMS, SEM, SKEWNESS, 
%	STATISTIC, STD, VAR, ZSCORE etc. 
%
% The mode is stored in the global variable FLAG_implicit_skip_nan
% It is recommended to use flag_implicit_skip_nan(1) as default and
% flag_implicit_skip_nan(0) should be used for exceptional cases only.
% This feature might disappear without further notice, so you should really not
% rely on it. 


%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

%	$Id: flag_implicit_skip_nan.m 8351 2011-06-24 17:35:07Z carandraug $
% 	Copyright (C) 2001-2003,2009 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


persistent FLAG_implicit_skip_nan;

%% if strcmp(version,'3.6'), FLAG_implicit_skip_nan=(1==1); end;	%% hack for the use with Freemat3.6

%%% set DEFAULT value of FLAG
if isempty(FLAG_implicit_skip_nan),
	FLAG_implicit_skip_nan = (1==1); %logical(1); % logical.m not available on 2.0.16
end;

FLAG = FLAG_implicit_skip_nan;
if nargin>0,
	FLAG_implicit_skip_nan = (i~=0); %logical(i); %logical.m not available in 2.0.16 
	if (~i)
		warning('flag_implicit_skipnan(0): You are warned!!! You have turned off skipping NaN in sumskipnan. This is not recommended. Make sure you really know what you do.')
	end;
end;

