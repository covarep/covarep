function fig2emf(h,s,p)
%FIG2EMF save a figure in windows metafile format (H,S,P)
% Usage:  (1) fig2emf      % save current figure in current folder
%         (2) fig2emf('',0)   % embolden the lines before saving
%         (3) fig2emf('',-2)   % embolden the lines and make width/height = 2
%         (4) fig2emf('C:/temp/<m>-<n>',0) embolden and save in C:\temp
%                                          with the m-file name and figure number
%
% Inputs: h   figure handle [use [] or omit for the current figure]
%         s   optional file name which can include <m> for the top level
%                 mfile name and <n> for figure number [if empty or missing: '<m>_<n>']
%                 '.' suppresses the save, if s ends in '/' or '\', then '<m>_<n>' is appended
%         p   argument for call to figbolden if required
%                 or 0 for default argument [ [] for no figbolden call]
%               [xmin ymin width height] gives the lower left corner position and the window size in pixels
%               [width height] leaves the lower left corner alone
%               [width] has a standard aspect ratio of 4:3
%               [-width/height] leaves the area unchanged but fixes the aspect ratio

%      Copyright (C) Mike Brookes 2012
%      Version: $Id: fig2emf.m 3103 2013-06-14 09:01:21Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 0
        h=[];
        s='';
        p=[];
    case 1
        if ischar(h) || ~numel(h)   % fig2emf(s)
            s=h;
            h=[];
        else
            s='';
        end
        p=[];
    case 2
        if ischar(h) || ~numel(h)   % fig2emf(s,p)
            p=s;
            s=h;
            h=[];
        else                        % fig2emf(h,s)
            p=[];
        end
end
if ~numel(h)
    h=gcf;
else
    figure(h);
end
if ~numel(s)
    s='<m>_<n>';
elseif s(end)=='/' || s(end)=='\'
    s=[s '<m>_<n>'];
end
[st,sti]=dbstack;
if numel(st)>1
    mfn=st(end).name;  % ancestor mfile name
else
    mfn='Figure';
end
fn=num2str(round(h));
ix=strfind(s,'<m>');
while numel(ix)>0
    s=[s(1:ix-1) mfn s(ix+3:end)];
    ix=strfind(s,'<m>');
end
ix=strfind(s,'<n>');
while numel(ix)>0
    s=[s(1:ix-1) fn s(ix+3:end)];
    ix=strfind(s,'<n>');
end
if numel(p)>0
    if numel(p)==1 && p==0
        figbolden;
    else
        figbolden(p)
    end
end
if ~strcmp(s,'.')
    eval(['print -dmeta ' s]);
end