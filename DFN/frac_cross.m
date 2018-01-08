%  Generates two perpendicular intersecting fractures
%  ---------------------------------------------------------------------
%  Copyright (C) 2016 by the Thermaid authors
% 
%  This file is part of Thermaid.
% 
%  Thermaid is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  Thermaid is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with Thermaid.  If not, see <http://www.gnu.org/licenses/>.
%  ---------------------------------------------------------------------
% 
%  Authors: Gunnar Jansen, University of Neuchatel, 2016-2017

% global dxf Nf_i frac_angle XY1

udata.N_fractures = 2;

dx_f = 1.1*min(udata.dx);

xc = [];
yc = [];
xb = [];
yb = [];
xe = [];
ye = [];

Nf_i = [];

dxi = 0;
dyi = dx_f;

x = 0.5*udata.len(1); y=2/9*udata.len(2);
r = floor(5/9.*min(udata.len(1),udata.len(2))./dx_f);
ybi = y:dyi:y+r*dyi-dyi; 
yei = y+dyi:dyi:y+r*dyi;

xci = 0.5*udata.len(1)*ones(size(yei));


xb =[xb round(xci .* 1e4) ./ 1e4];
xe =[xe round(xci .* 1e4) ./ 1e4];
yb =[yb round(ybi .* 1e4) ./ 1e4];
ye =[ye round(yei .* 1e4) ./ 1e4];

udata.Nf_i(1) = length(xci);

dxi = dx_f;
dyi = 0;

x = 2/9*udata.len(1); y=0.5*udata.len(2);
xbi = x:dxi:x+r*dxi-dxi; 
xei = x+dxi:dxi:x+r*dxi;

yci = 0.5*udata.len(2)*ones(size(xci));

udata.Nf_i(2) = length(xei);

xc =[xc round(xci .* 1e4) ./ 1e4];
yc =[yc round(yci .* 1e4) ./ 1e4];
xb =[xb round(xbi .* 1e4) ./ 1e4];
xe =[xe round(xei .* 1e4) ./ 1e4];
yb =[yb round(yci .* 1e4) ./ 1e4];
ye =[ye round(yci .* 1e4) ./ 1e4];

ListOfVariables = who;
for k = 1:length(ListOfVariables)
   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
end 

XY1 = [xb' yb' xe' ye'];
XY2(:,1,1) = xb';
XY2(:,1,2) = yb';
XY2(:,2,1) = xe';
XY2(:,2,2) = ye';
XE = [xe; ye]';

udata.Nf_f = length(xe);
udata.frac_angle = atan2((ye-yb),(xe-xb))*180/pi;
udata.dxf = sqrt((xb(1)-xe(1))^2 + (yb(1)-ye(1))^2);
