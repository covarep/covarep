function [ih,xa,ya]=imagehomog(im,h,m,clip)
%IMAGEHOMOG Apply a homography transformation to an image with bilinear interpolation
%Inputs: im(ny,nx,nc)  input image (uint8)
%        h(3,3)        homography
%        m             mode string
%                         i  show image [default if no output arguments]
%                         m  matlab coordinates (1,1) is top left [default]
%                         c  central coordinates (0,0) is the centre = {(1+nx)/2,(1+ny)/2} in 'm'
%                         k  clip to original image dimensions
%                         x  extend zero-fill to the specified clipping rectangle
%        clip(4)       bounding box [xmin xmax ymin ymax]
% Outputs:
%        ih(my,mx,nc)  output image (uint8)
%        xa(mx)        x axis
%        ya(my)        y axis

% Bugs/Suggestions:
% (1) cope with non-uint8
% (2) cope with (a) multiple inputs, (b) multiple transformations
% (3) output a boundary mask as an alpha channel
% (4) do anti-aliasing along the boundary
% (5) check that origin shift is correct for central coordinates

%      Copyright (C) Mike Brookes 2010
%      Version: $Id: imagehomog.m,v 1.3 2010/05/10 15:07:59 dmb Exp $
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

maxby=1e7;  % maximum memory to use
[ny,nx,nc]=size(im);
if nargin<4
    clip=[];
    if nargin<3
        m='';
        if nargin<2
            h=eye(3);
        end
    end
end
imr=reshape(im,nx*ny,nc);
t=eye(3);
if any(m=='c')   % convert homography and clipping box to matlab coordinates
    t(7:8)=0.5+[nx ny]/2;  % shift origin to image centre
    h=t*h/t;  % change homography so input and output use MATLAB coordinates
    if numel(clip)
        clip=clip+t([7 7 8 8]);  % make clipping use MATLAB coordinates as well
    end
end
box=h*[1 1 nx nx; 1 ny 1 ny; 1 1 1 1];
box=box(1:2,:)./box([3 3],:);
box=[min(box(1,:)) max(box(1,:)) min(box(2,:)) max(box(2,:))];
if any(m=='k')
    clip=[1 nx 1 ny];
end
if ~numel(clip)
    clip=box;
end
clip=clip(:)';
clip(1:2:3)=floor(clip(1:2:3));
clip(2:2:4)=ceil(clip(2:2:4));
box(1:2:3)=floor(max(clip(1:2:3),box(1:2:3)));  % no point in calculating non-existant points
box(2:2:4)=ceil(min(clip(2:2:4),box(2:2:4)));
g=inv(h);
mx=box(2)-box(1)+1; % number of columns in destination
my=box(4)-box(3)+1; % number of rows in destination

ih=zeros(my*mx,nc,'uint8');
ncol=max(1,min(mx,floor(maxby/(my*nc*8)))); % number of columns to do in a chunk
nloop=ceil(mx/ncol);
cmax=1+rem(mx-1,ncol); % final column of first iteration
jxinc=ncol*my;      % increment target indices each loop
ginc=g(:,1)*ncol; % increment transformed targets each loop
wj=ones(1,jxinc);   % repeat index for ginc
jx=1+cmax*my-jxinc:cmax*my;  % initial target indices (some might be negative)
gk=g*[reshape(repmat(cmax-ncol+box(1):cmax+box(1)-1,my,1),1,jxinc); repmat(box(3):box(4),1,ncol); ones(1,jxinc)];
gn=gk(1:2,:)./gk([3 3],:);   % normalize source coordinates
mn=[zeros(1,jxinc-cmax*my) ones(1,cmax*my)]; % mask for initial iteration
mn=mn & (gn(1,:)>-0.5 & gn(2,:)>-0.5 & gn(1,:)<nx+0.5 & gn(2,:)<ny+0.5); % mask valid pixels
w3=ones(nc,1);
for i=1:nloop
    fn=max(floor(gn(:,mn)),1);
    fn1=min(max(fn(1,:)',1),nx-1);
    fn2=min(max(fn(2,:)',1),ny-1);
    dn=gn(:,mn)-[fn1 fn2]';
    dn1=min(max(dn(1,:)',0),1);
    dn2=min(max(dn(2,:)',0),1);
    dn1c=1-dn1;
    dn2c=1-dn2;
    ih(jx(mn),:)=uint8(dn1c(:,w3).*(dn2c(:,w3).*single(imr(fn2+ny*(fn1-1),:)) ...
        +dn2(:,w3).*single(imr(fn2+1+ny*(fn1-1),:))) ...
        +dn1(:,w3).*(dn2c(:,w3).*single(imr(fn2+ny*(fn1),:)) ...
        +dn2(:,w3).*single(imr(fn2+1+ny*(fn1),:))));
    jx=jx+jxinc;  % target indices
    gk=gk+ginc(:,wj);
    gn=gk(1:2,:)./gk([3 3],:);   % normalize source coordinates
    mn=gn(1,:)>-0.5 & gn(2,:)>-0.5 & gn(1,:)<nx+0.5 & gn(2,:)<ny+0.5; % mask valid pixels
end
ih=reshape(ih,[my,mx,nc]);
if any(m=='x')          % extend blank area to specified clipping rectangle
    ih=[zeros(box(3)-clip(3),clip(2)-clip(1)+1,nc,'uint8'); ...
        zeros(my,box(1)-clip(1),nc,'uint8') ih zeros(my,clip(2)-box(2),nc,'uint8');
        zeros(clip(4)-box(4),clip(2)-clip(1)+1,nc,'uint8')];
    xa=(clip(1):clip(2))-t(7);
    ya=(clip(3):clip(4))-t(8);
else
    xa=(box(1):box(2))-t(7);
    ya=(box(3):box(4))-t(8);
end

if ~nargout || any(m=='i')
    imagesc(xa,ya,ih);
    axis equal
end

