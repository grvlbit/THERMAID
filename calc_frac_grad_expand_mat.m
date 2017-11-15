function frac_grad_expand_mat = calc_frac_grad_expand_mat(N_fractures,Nf_f, Nf_i)
%  Calculate the matrix operator for the (expanded) gradient on the fractures.
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
%  calc_frac_grad_expand_mat(N_fractures,Nf_f, Nf_i)
%
%  Input: 
%        N_fracture               number of fractures
%        Nf_f                     number of fracture segments
%        Nf_i(N_fractures,1)      fracture segments by fracture
%     
%  Output:
%        frac_grad_expand_mat
%                                 matrix operator for expanded gradient 
%                                 calculation for the fractures
%  
%  Usage:       
%          grad_input = frac_grad_expand_mat*input;
% 
%  This is utilized for example in the velocity calculation in the fractures 
%    (vf = -Tf.*(frac_grad_expand_mat*pf-gf);)
%

Ni=zeros(size(Nf_i));
for i=1:N_fractures
    for j=i:-1:1
        Ni(i) = Ni(i) +Nf_i(j);
    end
end
Nii = 1:1:length(Ni);
Nf_ff = Nf_f +N_fractures;

%% Expanding gradient matrix
T = ones(Nf_ff,1);
T(Ni+Nii-Nf_i) = 0;
T(Ni+Nii) = 0;

U = ones(Nf_ff,1);
U(Ni+Nii-Nf_i)= 0;
U(Ni+Nii) = 0;

Ds      = [-T U];

frac_grad_expand_mat = zeros(Nf_ff,Nf_f);
j=0;
k=1;
count =1;
for i=1:Nf_ff
    if i<Nf_ff, frac_grad_expand_mat(i,i-j)=Ds(i,2); end
    if i>1, frac_grad_expand_mat(i,i-j-1)=Ds(i,1); end
    if count==Nf_i(k)+1 && k<N_fractures
        j= j+1; count=0; k=k+1;
    end
    count = count +1;
end
frac_grad_expand_mat = sparse(frac_grad_expand_mat);

