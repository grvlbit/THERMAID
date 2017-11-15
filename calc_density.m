function [density, density_f] = calc_density(T,Tf,p,pf)
%  Calculate the pressure and temperature dependent density of the fluid.
%  The equation of state follows the work of Hongbing et al. (2008)
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
%
%  calc_density(T, Tf, p, pf)
%
%  Input: 
%        T     (nx,ny)         current temperature distribution in the matrix
%        Tf    (nf,1)          current temperature distribution in the fractures
%        p     (nx,ny)         current pressure distribution in the matrix
%        pf    (nf,1)          current pressure distribution in the fractures
%
%  Output:
%        density      (nx,ny)  current fluid density distribution in the matrix
%        density_f    (nf,1)   current fluid density distribution in the fractures

%  Change the dimensions from Pa to MPa for this function 
p  = max(p./1e6, 0.1);
pf = max(pf./1e6, 0.1);

%  Calculate fluid density density with equation of state from Hongbing et al. (2008)
density = 9.992e2 + 9.539e-2.*T - 7.618e-3.*T.*T + 4.336e-1.*p + 1.762e-3.*p.*p;
density_f = 9.992e2 + 9.539e-2.*Tf - 7.618e-3.*Tf.*Tf + 4.336e-1.*pf + 1.762e-3.*pf.*pf;