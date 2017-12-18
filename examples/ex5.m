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

addpath ../

close all; clear all;

THERMAID('Input_ex5',1)

global Nf

s = linspace(0,1,256);
rgb1=[0.230, 0.299, 0.754];
rgb2=[0.706, 0.016, 0.150];
cmap = diverging_map(s,rgb1,rgb2);


figure 
hold on;
pcolor(x,y,p'); shading interp
colormap(cmap)   
axis equal
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
hold off;

c=colorbar;
ylabel(c,'Pressure [Pa]')
xlabel('x [m]')
ylabel('x [m]')

figure 
hold on;
pcolor(x,y,tNew'); shading interp
colormap(cmap)
axis equal
[C,hfigc] = contour(x, y, tNew',10,'ShowText','off');
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
hold off;

c=colorbar;
ylabel(c,'Temperature [°C]')
xlabel('x [m]')
ylabel('x [m]')
