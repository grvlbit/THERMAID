function [frac_left_flux_mat, frac_right_flux_mat] = calc_frac_flux_mat(N_fractures,Nf_f, Nf_i)
%  Calculate the matrix operator for flux extraction on the fractures
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
%  calc_frac_flux_mat(N_fractures,Nf_f, Nf_i)
%
%  Input: 
%        N_fractures              number of fractures
%        Nf_f                     number of fracture segments
%        Nf_i(N_fractures,1)      fracture segments by fracture
%     
%  Output:
%        frac_left_flux_mat       matrix operator for left-side flux
%                                 calculation for the fractures
%        frac_right_flux_mat      matrix operator for right-side flux
%                                 calculation for the fractures

%% Extract left fluxes
Ni=zeros(size(Nf_i));
for i=1:N_fractures
    for j=i:-1:1
        Ni(i) = Ni(i) +Nf_i(j);
    end
end
Nii = 1:1:length(Ni);
Nf_ff = Nf_f +N_fractures;

T = ones(Nf_ff,1);
T(Ni-Nf_i+1) = 0;

frac_left_flux_mat = zeros(Nf_f,Nf_ff);
j=0;
k=1;
count =1;
for i=1:Nf_f
    frac_left_flux_mat(i,i+j) = T(i);
    if count==Nf_i(k) && k<N_fractures
        j= j+1; count=0; k=k+1;
    end
    count = count +1;
end
frac_left_flux_mat = sparse(frac_left_flux_mat);

%% Extract right fluxes
T = ones(Nf_ff,1);
T(Ni) = 0;

frac_right_flux_mat = zeros(Nf_f,Nf_ff);
j=0;
k=1;
count =1;
for i=1:Nf_f
    frac_right_flux_mat(i,i+j+1) = T(i);
    if count==Nf_i(k) && k<N_fractures
        j= j+1; count=0; k=k+1;
    end
    count = count +1;
end
frac_right_flux_mat = sparse(frac_right_flux_mat);



