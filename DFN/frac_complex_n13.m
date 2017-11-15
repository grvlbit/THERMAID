%  Generates a complex fracture network with 13 fractures
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

udata.N_fractures = 13;
Nf_i = [];
dxf = 1.0*sqrt(udata.dx(1)*udata.dx(1)+udata.dx(2)*udata.dx(2));

xc = []; yc = [];
xb = []; yb = [];
xe = []; ye = [];

beta = [-41 23 18 33 -19 -35 15 63 -73 45 46 47 74]; % Fracture angle [°]
xy = [50 450;
      50 350;
      50 240;
      50 50;
      100 300;
      150 150;
      200 100;
      260 150;
      300 450;
      300 50;
      350 400;
      325 250;
      400 300;
      ];
xy = xy./500;

L = [500 130 350 200 300 225 150 170 160 200 60 200 150];

L = L./500;
L = floor(L*min(udata.len(1),udata.len(2))./dxf);

for i=1:udata.N_fractures 
    pass = 0;
    while pass < 1
        
        x = xy(i,1)*udata.len(1);
        y = xy(i,2)*udata.len(2);
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
            all(ybi >= 0) && all(yei >= 0))
            pass = 1;
        end
    end
    
    udata.Nf_i(i) = length(xei);
    
    xb =[xb round(xbi,4)];
    xe =[xe round(xei,4)];
    yb =[yb round(ybi,4)];
    ye =[ye round(yei,4)];
end

udata.Nf_f = length(xe);

XY1 = [xb' yb' xe' ye'];
XE = [xe; ye]';

udata.frac_angle = atan2((ye-yb),(xe-xb))*180/pi;
udata.dxf = dxf;

% Comment in the following code to visualize the fracture network before
% the simulation is started
% -------------------------------------------------------------------------
% ListOfVariables = who;
% for k = 1:length(ListOfVariables)
%   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
% end 
% figure(221)
% line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','[0 0.4470 0.7410]','LineWidth',3);
% xlim([0 500]);
% ylim([0 500]);
% xlabel('x [m]');
% ylabel('y [m]');
% pause()
