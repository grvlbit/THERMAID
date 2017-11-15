function [vx,vy,vf,Vfm,Vmf,Vff] = calcVelocity(udata,p,pf,Tx,Ty,Tf,Tfm,Tff,g,gf)
%  Calculates the velocity field from the pressure
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
%           Ivan Lunati, Rouven Kuenze, University of Lausanne, 2012
%
%  Acknowledgement:  thanks are due to Manav Tyagi and Hadi Hajibeygi for
%                    contributing to the very early development of the code. 
%
%  calcVelocity(udata,p,pf,Tx,Ty,Tf,Tfm,Tff,g,gf)
%
%  Input: 
%        udata  [struct]        user data
%        p      (nx,ny)         matrix pressure
%        pf     (nf,1)          fracture pressure
%        Tx     (nx+1,ny)       trasmissibility in the x direction
%        Ty     (nx,ny+1)       trasmissibility in the y direction
%        Tf     (nf,1)          fracture transmissivity
%        Tfm    (nf,nx*ny)      coupling matrix between fracture and matrix
%        Tmf    (nx*ny,nf)      coupling matrix between matrix and fracture
%        Tff    (nf*nf,nf*nf)   fracture-fracture intersection
%                               transmissivity
%        g      (nx,ny)         gravity term with respect to the interfaces
%        gf     (nf,1)          gravity term with respect to the fracture interfaces
%  Output:
%        vx     (nx+1,ny)       matrix velocity in the x direction
%        vy     (nx,ny+1)       matrix velocity in the y direction
%        vf     (nf,1)          fracture velocity
%        Vfm    (nf,nx*ny)      velocity matrix between fracture and matrix
%        Vmf    (nx*ny,nf)      velocity matrix between matrix and fracture
%        Vff    (nf*nf,nf*nf)   fracture-fracture intersection velocity
 
persistent frac_grad_expand_mat
if isempty(frac_grad_expand_mat)
    N_fractures = length(udata.Nf_i);
    frac_grad_expand_mat = calc_frac_grad_expand_mat(N_fractures,udata.Nf_f, udata.Nf_i);
end

N = size(p);

%-------------------------------------------------------------------------%
%   Calculate velocity direction for alpha                                %
%-------------------------------------------------------------------------%

%-------------- velocity inside matrix cells -----------------------------%
vx = sparse(N(1)+1,N(2));
vy = sparse(N(1),N(2)+1);

vx(2:N(1),:) = -Tx(2:N(1),:).*(p(2:N(1),:)-p(1:N(1)-1,:));                 % velocity direction in x direction
vy(:,2:N(2)) = -Ty(:,2:N(2)).*(p(:,2:N(2))-p(:,1:N(2)-1)) - g(:,2:N(2));   % velocity direction in y direction

%--------------------- velocity on matrix boundaries ---------------------%                                 
  
vx(1,:) = -Tx(1,:).*(p(1,:)-udata.Fix(1:N(2))').*udata.ibcs(1:N(2))';                                                        % West Dirichlet
vx(1,:) = vx(1,:)+(udata.Fix(1:N(2)).*~udata.ibcs(1:N(2)))';                                                                 % West Neumann

vx(N(1)+1,:) = -Tx(N(1)+1,:).*(udata.Fix(N(2)+1:2*N(2))'-p(N(1),:)).*udata.ibcs(N(2)+1:2*N(2))';                             % East (adjust Neumann sign!)
vx(N(1)+1,:) = vx(N(1)+1,:)-(udata.Fix(N(2)+1:2*N(2)).*~udata.ibcs(N(2)+1:2*N(2)))';

vy(:,1) = (-Ty(:,1).*(p(:,1)-udata.Fix(2*N(2)+1:2*N(2)+N(1))) - g(:,1)).*udata.ibcs(2*N(2)+1:2*N(2)+N(1));                   % South
vy(:,1) = vy(:,1)+(udata.Fix(2*N(2)+1:2*N(2)+N(1)).*~udata.ibcs(2*N(2)+1:2*N(2)+N(1)));   

vy(:,N(2)+1) = (-Ty(:,N(2)+1).*(udata.Fix(2*N(2)+N(1)+1:2*N(2)+2*N(1))-p(:,N(2))) - g(:,N(2)+1)).*udata.ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1));   % North (adjust Neumann sign!)
vy(:,N(2)+1) = vy(:,N(2)+1)-(udata.Fix(2*N(2)+N(1)+1:2*N(2)+2*N(1)).*~udata.ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1)));

%--------------------- velocity in fractures -----------------------------%
vf = -Tf.*(frac_grad_expand_mat*pf-gf);

%--------------------- velocity fracture-fracture-------------------------%                                 
Vff = bsxfun(@times, pf, Tff) - bsxfun(@times,pf',Tff);

%--------------------- velocity fracture-matrix---------------------------% 
Vfm = bsxfun(@times, pf, Tfm) - bsxfun(@times, p(:)',Tfm);
%--------------------- velocity matrix-fracture---------------------------%
Vmf = -Vfm';

