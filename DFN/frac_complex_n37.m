%  Generates a complex fracture network with 37 fractures
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

udata.N_fractures = 37;
Nf_i = [];
dxf = 1.1*sqrt(udata.dx(1)*udata.dx(1)+udata.dx(2)*udata.dx(2));

xc = []; yc = [];
xb = []; yb = [];
xe = []; ye = [];


% Data sorted by y begin
xstart = [5 7.5 7.5 6 6 6 8.5 8 7.5 ... 
          9.5 9 9 9 11.5 11.5 11.5 11.5 10.5 ...
          0 0 0 0 0 0 0 0 0 ...
          2 1.5 3 3.5 3.5 4.5 5 4.5 5.5 ...
          13.5 12.5];
ystart = [0 8.5 8 3.5 2 0 9.5 5.5 2 ... 
          7 6.5 5 2.5 8 6.5 3.5 2.5 0 ...
          9.5 8.5 8.0 5.5 4 3 2 1.5 0.5 ...
          2.5 0.5 8.5 4 1 9.5 5.5 4 6.5 ...
          8.5 8.0];

xend = [7 10.5 8.5 9.5 10 10 10 10.5 11.5 ... 
        12 14.5 14.5 12.5 13.5 15 13 15 14.5 ...
        2 5 3 2 1.5 1.5 5 5 5  ...
        3 4 4 7 4.5 8.5 8 7 7.5 ... 
        14.5 15];
yend = [2 9.5 6 8.5 4.5 1.5 6.5 2.5 0 ... 
        10 4 6.5 6.5 10 9 5 1.5 3.5 ...    
        9.25 7.5 6 8 6 4.5 7 2.5 1.5 ... 
        4.5 0 10 6 0.5 8.5 4.5 1.5 9.5 ... 
        10 7];

beta = atan2((yend-ystart),(xend-xstart))*180/pi;

L = sqrt((yend/15-ystart/15).^2+(xend/15-xstart/15).^2);

% L = floor(L*150./dxf);
L = floor(L*min(udata.len(1),udata.len(2))./dxf);

xy = [xstart/15; ystart/10]';

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
    
    xb =[xb round(xbi .* 1e4) ./ 1e4];
    xe =[xe round(xei .* 1e4) ./ 1e4];
    yb =[yb round(ybi .* 1e4) ./ 1e4];
    ye =[ye round(yei .* 1e4) ./ 1e4];
end

udata.Nf_f = length(xe);

XY1 = [xb' yb' xe' ye'];
XE = [xe; ye]';

udata.frac_angle = atan2((ye-yb),(xe-xb))*180/pi;
udata.dxf = dxf;

% % Comment in the following code to visualize the fracture network before
% % the simulation is started
% % -------------------------------------------------------------------------
% ListOfVariables = who;
% for k = 1:length(ListOfVariables)
%   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
% end 
% figure(221)
% line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','[0 0.4470 0.7410]','LineWidth',3);
% xlim([0 len(1)]);
% ylim([0 len(2)]);
% xlabel('x [m]');
% ylabel('y [m]');
% axis equal
% %pause()
