%  Generates two perpendicular intersecting fractures rotated 45 degrees
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
%  Authors: Gunnar Jansen, University of Neuchatel, 2016

Nf_i = [];
dxf = 1.2*sqrt(udata.dx(1)*udata.dx(1)+udata.dx(2)*udata.dx(2));

xc = []; yc = [];
xb = []; yb = [];
xe = []; ye = [];

N = 2;

rot = [-45 45]; %Fracture rotation [°] (Must be in interval (0,90)

midx = 0.5*udata.len(1);
midy = 0.5*udata.len(2);

L = 0.75*min(udata.len(1),udata.len(2));

xstart = [midx - 0.5*L.*cos(rot*pi/180)];
ystart = [midy - 0.5*L.*sin(rot*pi/180)];

xend = [midx + 0.5*L.*cos(rot*pi/180)];
yend = [midy + 0.5*L.*sin(rot*pi/180)];

beta = atan2((yend-ystart),(xend-xstart))*180/pi;

L = sqrt((yend-ystart).^2+(xend-xstart).^2);

L = floor(L./dxf);

xy = [xstart; ystart]';

for i=1:N
    if L(i) > 8
        pass = 0;
        x = xy(i,1);
        y = xy(i,2);
        r = L(i);
        
        dyi = sin(beta(i)*pi/180)*dxf;
        dxi = sqrt(dxf*dxf-dyi*dyi);
        
        xbi = x:dxi:x+r*dxi-dxi;
        xei = x+dxi:dxi:x+r*dxi;
        ybi = y:dyi:y+r*dyi-dyi;
        yei = y+dyi:dyi:y+r*dyi;
        
        if (all(xbi <= udata.len(1)) && all(xei <= udata.len(1)) && ...
                all(ybi <= udata.len(2)) && all(yei <= udata.len(2)) && ...
                all(xbi >= 0) && all(xei >= 0) && ...
                all(ybi >= 0) && all(yei >= 0) && ...
                length(ybi) == length(xbi) && ...
                length(ybi) == length(yei) && ...
                length(yei) == length(xei) && ...
                length(yei) == length(xbi) )
            pass = 1;
            
            udata.Nf_i(i) = length(xei);
            
            xb =[xb round(xbi,4)];
            xe =[xe round(xei,4)];
            yb =[yb round(ybi,4)];
            ye =[ye round(yei,4)];
            
        end
        
    end
end

udata.Nf_f = length(xe);
udata.N_fractures = length(udata.Nf_i);

XY1 = [xb' yb' xe' ye'];
XE = [xe; ye]';

udata.frac_angle = atan2((ye-yb),(xe-xb))*180/pi;
udata.dxf = dxf;

% Comment in the following code to visualize the fracture network before
% the simulation is started
% -------------------------------------------------------------------------
% ListOfVariables = who;
% for k = 1:length(ListOfVariables)
%     assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
% end
% figure(222)
% line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','[0.8500    0.3250    0.0980]','LineWidth',1.5);
% hold on
% % line([XYp(:,1)';XYp(:,3)'],[XYp(:,2)';XYp(:,4)'],'Color','[0 0.4470 0.7410]','LineWidth',1.5);
% 
% xlim([0 udata.len(1)]);
% ylim([0 udata.len(2)]);
% xlabel('x [m]');
% ylabel('y [m]');
% axis equal