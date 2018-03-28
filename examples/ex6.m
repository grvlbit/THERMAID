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

%% THERMAID

global k_ratio

k_ratio = 1e2;
THERMAID('Input_ex6',0)

figure(1); 
pdiag = zeros(udata.Nf(1),1);
dx_diag = sqrt(udata.dx(1)*udata.dx(1)+udata.dx(2)*udata.dx(2));
xdiag = 0:dx_diag:(udata.Nf(1)-1)*dx_diag;
% xdiag = dx_diag:dx_diag:(udata.Nf(1))*dx_diag;

for i = 1:udata.Nf(1)
    pdiag(i) = p(i,i);
end


%% Analytical
a = 45*pi/180; % Fracture angle °
b_max = 0.05;  % Maximum fracture aperture [m]
L = 2; 
K_m = 1e-12;
K_f = K_m*k_ratio;
mu = 1e-3;
q0 = 1e-4;

x = linspace(-5,5,301);
y = linspace(-5,5,301);
[X,Y] = meshgrid(x,y);
z = X+1i*Y;
z1 = (-1+0i)*exp(1i*a);
z2 = (+1+0i)*exp(1i*a);

Z = (z-0.5*(z1+z2))/(0.5*(z2-z1));
A = 0.5*K_f*b_max/(K_m*L+K_f*b_max)*q0*L*cos(a);

l = 1;
o1 = -A.*sign(real(Z)).*sqrt((Z-l).*(Z+l));
o2 = A.*Z;
o3 = -0.5*q0*L.*exp(1i*a).*Z;

omega = o1 + o2 + o3;
phi = real(omega)*mu/K_m;

phidiag = zeros(size(phi,1),1);
dx = X(1,2)-X(1,1);
dy = Y(2,1)-Y(1,1);
dx_phi = sqrt(dx*dx+dy*dy);
xphi = 0:dx_phi:(size(phi,1)-1)*dx_phi;

for i = 1:size(phi,1)
    phidiag(i) = phi(i,i);
end

%% Plot

p = plot(xdiag,pdiag/1e6,'ko',xphi,phidiag/1e6,'r','LineWidth',2.0,'MarkerSize',8);
nummarkers(p,55);
xlabel('Distance from bottom left corner [m]')
ylabel('Matrix Pressure [MPa]')
legend(p,'THERMAID','Analytical solution')

ListOfVariables = who;
for k = 1:length(ListOfVariables)
   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
end
