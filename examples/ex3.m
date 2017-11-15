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

global k_ratio

% Run with Kf/Km = 1e5
k_ratio = 1e5;

THERMAID('Input_ex3',1)

%% Plots
global Nf

figure(1) 
hold on;
pcolor(x,y,p'); shading interp
[C,hfigc] = contour(x, y, p',[0 1e6 2e6 3e6 4e6 5e6 6e6 7e6 8e6 9e6 10e6],'ShowText','off');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
hold off;

c=colorbar;
ylabel(c,'Pressure [Pa]')
xlabel('x [m]')
ylabel('x [m]')
colormap jet
axis equal
ng;


figure(2) 
hold on;
pcolor(x,y,tNew'); shading interp
%axis square
[C,hfigc] = contour(x, y, tNew',[50 70 90 110 130 150 170 ],'ShowText','on');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
hold off;

c=colorbar;
ylabel(c,'Temperature [°C]')
xlabel('x [m]')
ylabel('x [m]')
colormap jet
axis equal
ng;


rmpath ../
