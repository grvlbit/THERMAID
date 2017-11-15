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

THERMAID('Input_ex5',1)

global Nf

figure 
hold on;
pcolor(x,y,p'); shading interp
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

E = 2*shear_modulus*(1+poisson_ratio);
deltaT = tNew-T0.*ones(size(tNew));
dsigma_t = E/ (1 - 2 * poisson_ratio)*therm_exp_coeff*deltaT;

deltaTf = tNewf-T0f.*ones(size(tNewf));
dsigma_tf = E/ (1 - 2 * poisson_ratio)*therm_exp_coeff*deltaTf;

dsigma_t = (dsigma_t/1e6);
dsigma_tf =(dsigma_tf/1e6);

figure 
hold on;
pcolor(x,y,dsigma_t'); shading interp
axis equal
[C,hfigc] = contour(x, y, dsigma_t',10,'ShowText','on');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
scatter(xe,ye,20,dsigma_tf,'filled')
hold off;
c=colorbar;
ylabel(c,'Thermal stress [MPa]')
xlabel('x [m]')
ylabel('x [m]')

[vfx,vfy] = calcFractureVelocity(XY1,vf);
figure
quiver(x',y',vx(1:Nf(1),:)',vy(:,1:Nf(2))')
hold on; 
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color',[0.8 0.8 0.8],'LineWidth',4);
quiver(xb',yb',vfx,vfy,0.5,'LineWidth',1.5)
hold off
axis equal

rmpath ../
