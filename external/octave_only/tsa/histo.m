function [H,X]=histo(Y,Mode)
% HISTO calculates histogram for each column
% [H,X] = HISTO(Y,Mode)
% 	   
%   Mode 
%	'rows' : frequency of each row
%	'1x'   : single bin-values 
%	'nx'   : separate bin-values for each column
%   X  are the bin-values 
%   H  is the frequency of occurence of value X 
%
% HISTO(Y) with no output arguments:
%	plots the histogram bar(X,H)
%
% more histogram-based results can be obtained by HIST2RES2  
%
% see also:  HISTO, HISTO2, HISTO3, HISTO4
%
% REFERENCE(S):
%  C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).

%	$Id: histo.m 5148 2008-06-27 10:36:37Z schloegl $
%	Copyright (C) 1996-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the TSA-toolbox 
%	http://hci.tugraz.at/~schloegl/matlab/tsa/
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin<2,
        Mode='1x';
end;
Mode=lower(Mode);

if strcmp(Mode,'rows')
        R = histo4(Y);
        
elseif strcmp(Mode,'column')
        R = histo4(Y');
        R.X = R.X';
        
elseif strcmp(Mode,'1x')
        R = histo3(Y);
        
elseif strcmp(Mode,'nx')
        R = histo2(Y);
        
end;

H = R.H;
X = R.X;
if nargout == 0,
        if any(size(X)==1),
                if exist('OCTAVE_VERSION')<5,
                        bar(R.X,R.H,'stacked');
                else
                        bar(R.X,R.H);   
                end
        else
                warning('2-dim X-values not supported\n')
		%bar3(R.X,R.H);
        end;
end;
