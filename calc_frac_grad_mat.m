function frac_grad_mat = calc_frac_grad_mat(N_fractures,Nf_f, Nf_i)
%  Calculate the matrix operator for the gradient on the fractures
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
%  calc_frac_grad_mat(N_fractures,Nf_f, Nf_i)
%
%  Input: 
%        N_fracture               number of fractures
%        Nf_f                     number of fracture segments
%        Nf_i(N_fractures,1)      fracture segments by fracture
%     
%  Output:
%        frac_grad_mat(Nf_f,Nf_f+N_fractures)
%                                 matrix operator for gradient calculation
%                                 for the fractures
%  
%  Usage:       
%          grad_input = frac_grad_mat*input;

Nf_ff = Nf_f +N_fractures;

%% Gradient Matrix
U = ones(Nf_f,1);
T = ones(Nf_f,1);
D = [-T U];

frac_grad_mat = zeros(Nf_f,Nf_ff);
j=0;
k=1;
count =1;
for i=1:Nf_f
    frac_grad_mat(i,i+j)=D(i,1);
    frac_grad_mat(i,i+j+1)=D(i,2);
    if count==Nf_i(k) && k < N_fractures
        j= j+1; count=0; k=k+1;
    end
    count = count +1;
end
frac_grad_mat = sparse(frac_grad_mat);
