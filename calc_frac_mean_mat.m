function frac_mean_mat = calc_frac_mean_mat(N_fractures,Nf_f, Nf_i)
%  Calculate the matrix operator for the mean on the fractures
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
%  calc_frac_mean_mat(N_fractures,Nf_f, Nf_i)
%
%  Input: 
%        N_fracture               number of fractures
%        Nf_f                     number of fracture segments
%        Nf_i(N_fractures,1)      fracture segments by fracture
%     
%  Output:
%        frac_mean_mat(Nf_f+N_fractures,Nf_f)
%                                 matrix operator for mean calculation
%                                 for the fractures
%  
%  Usage:       
%        - Arithmetic Mean
%                          arit_mean = 0.5*frac_mean_mat*input;
% 
%        - Harmonic Mean
%                          harm_mean = 2./(frac_mean_mat*(1./input));

Ni=zeros(size(Nf_i));
for i=1:N_fractures
    for j=i:-1:1
        Ni(i) = Ni(i) +Nf_i(j);
    end
end
Nii = 1:1:length(Ni);
Nf_ff = Nf_f +N_fractures;

%% Averaging Matrix
T = ones(Nf_ff,1);
T(Ni+Nii-Nf_i) = 0;
T(Ni+Nii) = 2;

U = ones(Nf_ff,1);
U(Ni+Nii-Nf_i)= 2;
U(Ni+Nii) = 0;

frac_mean_mat = zeros(Nf_ff,Nf_f);
j=0; k=1; count =1;
for i=1:Nf_ff
    if i<Nf_ff, frac_mean_mat(i,i-j)=U(i); end
    if i>1, frac_mean_mat(i,i-j-1)=T(i); end
    if count==Nf_i(k)+1 && k<N_fractures
        j= j+1; count=0; k=k+1;
    end
    count = count +1;
end
frac_mean_mat = sparse(frac_mean_mat);
