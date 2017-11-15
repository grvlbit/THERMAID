function d_mean = calc_d_mean(dxx, dyy)
%  Calculate the mean distance <d> within a cell to a fracture segment
%  Please see the section 'Fracture Matrix Coupling' in the manual for more
%  detail.
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
%  calc_d_mean(dxx, dyy)
%
%  Input: 
%        dxx      x-grid spacing
%        dyy      y-grid spacing
%
%  Output:
%        d_mean   average distance between fracture and matrix

dx = [dxx dyy];
d_mean = (dx(1)*dx(2))/(3*sqrt(dx(1)*dx(1)+dx(2)*dx(2)));

if (3*sqrt(dx(1)*dx(1)+dx(2)*dx(2)) == 0)
    d_mean = 0;
end