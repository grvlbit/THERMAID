function [tNew] = transport_heat_System(udata,tN,tNf,vx,vy,vf,Vfm,Vmf,Vff,Dfm,Q,QT,density_l,density_lf)
%  Assembles and solves the transport system. It consists of the
%  subfunctions for advective and diffusive transport
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
%  Acknowledgement:  Thanks are due to Manav Tyagi and Hadi Hajibeygi for
%                    contributing to the very early development of the code.
%
%  transport_heat_System(udata,tN,tNf,vx,vy,vf,Vfm,Vmf,Vff,Dfm,Q,QT,density_l,density_lf)
%
%  Input: 
%        udata  [struct]        user data
%        tN     (nx,ny)         old heat transport solution on the matrix
%        tNf    (nf,1)          old heat transport solution in the fractures
%        vx     (nx+1,ny)       matrix velocity in the x direction
%        vy     (nx,ny+1)       matrix velocity in the y direction
%        vf     (nf,1)          fracture velocity
%        Vfm    (nf,nx*ny)      velocity matrix between fracture and matrix
%        Vmf    (nx*ny,nf)      velocity matrix between matrix and fracture
%        Vff    (nf*nf,nf*nf)   fracture-fracture intersection velocity
%        Dfm    (nf,nx*ny)      fracture-matrix diffusivity
%        Q      (nx,ny)         source term of the pressure problem
%        QT     (nx,ny)         source term of the temperature problem
%        density_l (nx,ny)      fluid density of matrix grid cells
%        densisty_lf (nf,1)     fluid density of the fracture segments
%
%  Output: 
%        tNew  (nx*ny+nf,1)     new heat transport solution vector

%-------------------------------------------------------------------------%
% Calculate matrix interface permeability & porosity (harmonic average)
%-------------------------------------------------------------------------%
[phix, phiy] = calc_interface_values(udata.phi);
[phif] = calc_interface_values_fracture(udata,udata.phi_f);
[Bf] = calc_interface_values_fracture(udata,udata.b0);

%-------------------------------------------------------------------------%
%    Build Transport System                                               %
%-------------------------------------------------------------------------%
[Up,r]  = transport_heat_Advection(udata,vx,vy,vf,Vfm,Vmf,Vff,Q,QT,tN,tNf,udata.cp_l.*density_l,udata.cp_l.*density_lf);        % Convective Upwind Matrix

if (udata.phi*udata.lambda_l+(1-udata.phi)*udata.lambda_s) == 0;
   Di = sparse(prod(udata.Nf)+udata.Nf_f,prod(udata.Nf)+udata.Nf_f);
   ri = sparse(prod(udata.Nf)+udata.Nf_f,1);
else
  [Di,ri] = transport_heat_Diffusion(udata,phix,phiy,phif,Bf,Dfm);         % Diffusive Matrix 
end

Ac      = sparse((udata.phi.*udata.cp_l.*density_l+(1-udata.phi).*udata.cp_s.*udata.density_s).*prod(udata.dx)/udata.dt);                                        % Accumulation term matrix
Acf = sparse((udata.phi_f.*udata.cp_l.*density_lf+(1-udata.phi_f).*udata.cp_s.*udata.density_sf).*udata.dxf.*udata.b0/udata.dt);                                % Accumulation term fracture

Ac = [Ac(:); Acf(:)];                                                      % Combine the accumulation terms
snnf = [tN(:); tNf(:)];                                                    % for matrix and fracture

A       = Up + diag(Ac(:)) + Di;                                           % Stiffness matrix
rhs     = r + ri + sparse(snnf(:).*Ac(:));                                 % Right hand side

%-------------------------------------------------------------------------%
%    Solve Transport System                                               %
%-------------------------------------------------------------------------%
tNew = A\rhs;
