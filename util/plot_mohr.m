function fig = plot_mohr(sigma,name)
%  Plot Mohr circle for a given stress tensor
%  ---------------------------------------------------------------------
%  Copyright (C) 2017 by the Thermaid authors
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
%  Authors: Gunnar Jansen, University of Neuchatel, 2017 

% Plot Mohr's circles for given Stress ratio
if (nargin<2), 
    name = 'Mohr´s circles' ; end

th = linspace( pi/2, -pi/2, 100);

fig = figure();
axis equal;
hold on

smean_1 = (sigma(1)+sigma(3))/2;
R1 = abs(sigma(3)-sigma(1))/2;
y = R1*cos(th) + 0; % Center of circle is s_mean
x = R1*sin(th) + smean_1; % No offset in y
p = plot(x,y);
set(p,'LineWidth',2,'Color','k');

smean_2 = (sigma(2)+sigma(3))/2;
R2 = abs(sigma(3)-sigma(2))/2;
y = R2*cos(th) + 0; % Center of circle is s_mean
x = R2*sin(th) + smean_2; % No offset in y
p = plot(x,y);
set(p,'LineWidth',2,'Color','k');

smean_2 = (sigma(2)+sigma(1))/2;
R3 = abs(sigma(1)-sigma(2))/2;
y = R3*cos(th) + 0; % Center of circle is s_mean
x = R3*sin(th) + smean_2; % No offset in y
p = plot(x,y);
set(p,'LineWidth',2,'Color','k');

hold off

xlim([sigma(3)-0.1 sigma(1)+0.1])
ylim([0 max(max(R1,R2),R3)+0.1])

text(sigma(1)+0.01,0.01,'\sigma_1');
text(sigma(2)+0.01,0.01,'\sigma_2');
text(sigma(3)+0.01,0.01,'\sigma_3');

title(name)
