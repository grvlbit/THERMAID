function [Tx,Ty,Tf,Tfm,Tff,DfmT,g,gf,density_l,density_lf] = initialize(udata,T,Tf,p,pf,CI,row,col)
%  Initialize updates the transmissivity based on pressure, concentraion and
%  temperature
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
%  initialize(udata,c,cf,T,Tf,p,pf,CI,row,col)
%
%  Input:
%        udata  [struct]        user data
%        T      (nx,ny)         matrix temperature
%        Tf     (nf,1)          fracture temperature
%        p      (nx,ny)         matrix pressure
%        pf     (nf,1)          fracture pressure
%        CI     (#intersectedCells,3) connectivity index (matrix-fracture)
%        row    (#ffIntersections,1) row indices of fracture-fracture
%                               intersections
%        col    (#ffIntersections,1) col indices of fracture-fracture
%                               intersections
%
%  Output:
%        Tx     (nx+1,ny)       trasmissibility in the x direction
%        Ty     (nx,ny+1)       trasmissibility in the y direction
%        Tf     (nf,1)          fracture transmissivity
%        Tfm    (nf,nx*ny)      fracture-matrix transmissivity
%        Tff    (nf,nf)         fracture-fracture transmissivity
%        DfmT   (nf,nx*ny)      fracture-matrix heat diffusivity
%        g      (nx,ny)         gravity term with respect to the interfaces
%        gf     (nf,1)          gravity term with respect to the fracture interfaces
%        density_l (nx,ny)      fluid density of matrix grid cells
%        densisty_lf (nf,1)     fluid density of the fracture segments

persistent frac_mean_m
if isempty(frac_mean_m)
    frac_mean_m = calc_frac_mean_mat(udata.N_fractures,udata.Nf_f, udata.Nf_i);
end

N = size(T);

%-------------------------------------------------------------------------%
%    Initialize                                                           %
%-------------------------------------------------------------------------%
if(~udata.const_viscosity)
    viscosity = 2.414e-5*10.^(247.8./(T+273.15-140));
    viscosity_f = 2.414e-5*10.^(247.8./(Tf+273.15-140));
else
    viscosity = udata.viscosity;
    viscosity_f = udata.viscosity_f;
end

if(~udata.const_density)
    [density_l, density_lf] = calc_density(T,Tf,p,pf);
else
    density_l = udata.density_l;
    density_lf = udata.density_lf;
end

%-------------------------------------------------------------------------%
% Calculate matrix interface permeability & porosity (harmonic average)
%-------------------------------------------------------------------------%
[Kx, Ky] = calc_interface_values(udata.K);
[Kf] = calc_interface_values_fracture(udata,udata.K_f);
[Bf] = calc_interface_values_fracture(udata,udata.b0);


[viscosityx, viscosityy] = calc_interface_values(viscosity);
[viscosityf] = calc_interface_values_fracture(udata,viscosity_f);

[~, densityy] = calc_interface_values(density_l);
[densityf] = calc_interface_values_fracture(udata,density_lf);

%-------------------------------------------------------------------------%
%    Gravity for the matrix                                               %
%-------------------------------------------------------------------------%

Gy    = densityy.*udata.gravity;
g            = Ky./viscosityy.*Gy.*udata.dx(1);  % Gravity term with respect to the interfaces
g(:,1)       = g(:,1).*udata.ibcs(2*N(2)+1:2*N(2)+N(1));
g(:,N(2)+1)  = g(:,N(2)+1).*udata.ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1));

%-------------------------------------------------------------------------%
%    Gravity for the fractures                                            %
%-------------------------------------------------------------------------%
interface_frac_angle = 0.5*frac_mean_m*udata.frac_angle';
Gyf= densityf.*udata.gravity.*sin(interface_frac_angle*pi/180);

gf           = -Kf./viscosityf.*Gyf.*Bf;
Ni=zeros(size(udata.Nf_i));
N_fractures = length(udata.Nf_i);
for i=1:N_fractures
    for j=i:-1:1
        Ni(i) = Ni(i) +udata.Nf_i(j);
    end
end
Nii = 1:1:length(Ni);

gf(Ni-udata.Nf_i+Nii) = 0;
gf(Ni+Nii) = 0;

%-------------------------------------------------------------------------%
%    Transmissivities                                                     %
%-------------------------------------------------------------------------%
Tx           = Kx./viscosityx.*udata.dx(2)./udata.dx(1);                                     % Transmissivities
Ty           = Ky./viscosityy.*udata.dx(1)./udata.dx(2);
Tf           = Kf./viscosityf./udata.dxf.*Bf;

%% Matrix-Fracture transmissivity
%  This assembles the matrix-fracture transmissivity based on the
%  previously computed connectivity index
l = udata.Nf_f;
n = udata.Nf;
X = zeros(length(CI),1);
for i = 1:length(CI)
    indm = CI(i,1);
    indf = CI(i,2);
    lambda_ij = udata.K(indm)/viscosity(indm);
    lambda_k = udata.K_f(indf)/viscosity_f(indf);
    X(i)  =  CI(i,3)*2*lambda_ij*lambda_k/(lambda_ij+lambda_k);   %Harmonic mean for the fracture-matrix transmissivity
    % X(i)  =  CI(i,3)*0.5*(lambda_ij+lambda_k);                  % (Optional) Arithmetic mean for the fracture-matrix transmissivity
    % X(i)  =  CI(i,3)*sqrt(lambda_ij+lambda_k);                  % (Optional) Geometric mean for the fracture-matrix transmissivity
    [ii,jj] = ind2sub(udata.Nf,CI(i,1));
    I(i) = (jj-1)*n(1)+ii;
end
if (isempty(CI))
    Tfm = sparse(zeros(n(1)*n(2),l));
else
    Tfm = sparse(I,CI(:,2),X,n(1)*n(2),l);
end
Tfm = Tfm';

%% Matrix-Fracture heat difffusivity
%  This assembles the matrix-fracture diffusivity based on the
%  previously computed connectivity index
if (udata.lambda_l == 0 || udata.lambda_s == 0)
    DfmT = sparse(zeros(n(1)*n(2),l));
else
    l = udata.Nf_f;
    n = udata.Nf;
    X = zeros(length(CI),1);
    for i = 1:length(CI)
        indm = CI(i,1);
        indf = CI(i,2);
        lambda_ij = (udata.phi(indm).*udata.lambda_l+(1-udata.phi(indm)).*udata.lambda_s);
        lambda_k = (udata.phi_f(indf).*udata.lambda_l+(1-udata.phi_f(indf)).*udata.lambda_s);
        X(i)  =  CI(i,3)*2*lambda_ij*lambda_k/(lambda_ij+lambda_k);                  %Harmonic mean for the fracture-matrix transmissivity
        [ii,jj] = ind2sub(udata.Nf,CI(i,1));
        I(i) = (jj-1)*n(1)+ii;
    end
    if (isempty(CI))
        DfmT = sparse(zeros(n(1)*n(2),l));
    else
        DfmT = sparse(I,CI(:,2),X,n(1)*n(2),l);
    end
end
DfmT = DfmT';

%% Fracture-Fracture transmissivity
%  This assembles the fracture-fracture transmissivity based on the
%  intersection points between the fracture segments
alpha_row = udata.b0(row).*udata.K_f(row)./viscosityf(row)./(0.5*udata.dxf);
alpha_col = udata.b0(col).*udata.K_f(col)./viscosityf(col)./(0.5*udata.dxf);
Tflk = alpha_row.*alpha_col./(alpha_row+alpha_col);
Tff = sparse(row,col,Tflk,l,l);
