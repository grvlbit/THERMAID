function [psif] = calc_interface_values_fracture(udata,psi_f)
%  Calculate harmonic mean interface values for the fracture interfaces
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
%  calc_interface_values_fracture(udata,psi_f)
%
%  Input: 
%         udata  [struct]    User data
%         psi    (Nf_f,1)    Source array
%  Output:
%         psif   (Nf_f+N_fractures,1)
%                            Harmonic mean interface values of the input
% 

persistent frac_mean_mat
if isempty(frac_mean_mat)
    N_fractures = length(udata.Nf_i);
    frac_mean_mat = calc_frac_mean_mat(N_fractures,udata.Nf_f, udata.Nf_i);
end

psif = 2./(frac_mean_mat*(1./psi_f));
