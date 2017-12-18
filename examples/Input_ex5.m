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

%% GRID PARAMETERS ---------------------------------------------------------------------------------%
udata.len     = [500 500];                    % physical length of the domain in x and y direction [m]
udata.Nf      = [151 151];  	                % number of cells in x and y direction
udata.dx      = udata.len./udata.Nf;                      % cell length [m]

%% SIMULATION PARAMETER FOR TRANSPORT --------------------------------------------------------------%
udata.timeSim  = 40*365*86400;                         % total simulation time [s]
udata.dt       = 1*365*86400;                           % time step length [s]
udata.tol    = 1.e-3;                       % saturation tolerance on pressure-concentration-heat loop [-]
udata.maxit    = 100;                         % maximum number of pressure concentration-heat loops to converge

%% FRACTURE NETWORK
frac_load_fracman

if (udata.dxf < min(udata.dx))
    error('dxf < dx')
end

%% INITIAL CONDITIONS--------------------------------------------------------------------------------%
udata.T0   = 200*ones(udata.Nf(1),udata.Nf(2));              % Initial matrix temperature [°C]
udata.T0f  = 200*ones(udata.Nf_f,1);                   % Initial fracture temperature [°C]
udata.tmax = 200;                              % maximum temperature [°C] for plotting

udata.p0   = zeros(udata.Nf(1),udata.Nf(2));              % Initial matrix pressure [Pa]
udata.p0f  = zeros(udata.Nf_f,1);                   % Initial fracture pressure [Pa]

%% BC FLUID ----------------------------------------------------------------------------------------%
udata.ibcs = zeros(2*sum(udata.Nf),1);              % type 0:Neumann(N); 1:Dirichlet(D)
udata.Fix  = zeros(2*sum(udata.Nf),1);              % value N [m2/s] (inflow>0); D [Pa]   

udata.ibcs(round(udata.Nf(2)/2):udata.Nf(2)) = 1;
udata.Fix(round(udata.Nf(2)/2):udata.Nf(2)) = 2.5e7;

udata.ibcs(udata.Nf(2)+1:2*udata.Nf(2)) = 1;
udata.Fix(udata.Nf(2)+1:2*udata.Nf(2)) = 0.0;

%% BC TRANSPORT ------------------------------------------------------------------------------------%
udata.flagHeatTransport = 1;

udata.FixT     = 200*ones(2*sum(udata.Nf),1);           % normalized concentration of boundary flow [-]
udata.FixT(round(udata.Nf(2)/2):udata.Nf(2)) = 80;

%% SOURCE TERMS ------------------------------------------------------------------------------------%
Q       = zeros(udata.Nf);                    % source term [m2/s]; inflow positive
QT      = zeros(udata.Nf);                    % normalized concentration for source term [-] 

%% GRAVITY---------------------------------------------------------------------------------%
udata.gravity   = 0;                          % gravity acceleration in y [m/s2]

%% PERMEABILITY ------------------------------------------------------------------------------------%
udata.K       = ones(udata.Nf(1),udata.Nf(2))*1e-15; 	    % permeability field [m2]
udata.K_f     = ones(udata.Nf_f,1)*1e-11; 	        % fracture permeability field [m2]

%% Fracture aperture
udata.b0      = sqrt(12.*udata.K_f).*ones(udata.Nf_f,1);  % fracture aperture field [m]

%% Porosity ----------------------------------------------------------------------------------------%
udata.phi     = ones(udata.Nf(1),udata.Nf(2))*0.05;	    % porosity field
udata.phi_f   = ones(udata.Nf_f,1)*0.5;             % fracture porosity field

%% Fluid density ------------------------------------------------------------------------%
udata.const_density = 1;
if(udata.const_density)
    udata.density_l  = 1000*ones(udata.Nf);       % density of the rock [kg/m3]
    udata.density_lf = 1000*ones(udata.Nf_f,1); % density of the rock [kg/m3]
end

%% Fluid viscosity ------------------------------------------------------------------------%
udata.const_viscosity = 1;
if(udata.const_viscosity)
    udata.viscosity   = 1e-3*ones(udata.Nf);       % density of the rock [kg/m3]
    udata.viscosity_f = 1e-3*ones(udata.Nf_f,1); % density of the rock [kg/m3]
end

%% MECHANIC PROPERTIES
udata.flagIncompressible = 0;                          % If 1 - incompressible fluids are used
udata.compressibility_l = 5e-100;                        % Bulk Modulus of the fluid [Pa]
udata.compressibility_s = 5e-100;                        % Bulk Modulus of the rock  [Pa]

udata.shear_modulus = 27e9;                        % Shear modulus of the rock [Pa]
udata.poisson_ratio = 0.25;                        % Poisson ratio of the rock [-]

udata.friction_coeff = 0.6;                        % Friction coefficient (shear failure) [-]
udata.K_enh = 1e2;				     % Permeability enhancement factor [-]

udata.therm_exp_coeff = 7.9e-6;                    % Thermal expansion coefficient of the rock matrix [-]

udata.flagFracStability = 0;   % If 0 the remaining parameters in this section are not used

udata.sigma_1 = 20e6*ones(udata.Nf_f,1);          % Maximum principal stress [Pa]
udata.sigma_2 = 14e6*ones(udata.Nf_f,1);          % Intermeadiate principal stress [Pa]
udata.sigma_3 = 12e6*ones(udata.Nf_f,1);          % Minimum principal stress [Pa]

udata.stress_trend  = [0 90];  % [TR(S1) TR(S3)]
udata.stress_plunge = 90;      % PL(S1)

udata.frac_az  = 90*ones(size(udata.frac_angle));  % Fracture dip [°]
udata.frac_dip = min(abs(udata.frac_angle),180-abs(udata.frac_angle)); % Fracture azimuth (from N) [°]

udata.use_thermal_stress = 0; % Boolean for thermal stress in calculations

%% ROCK DENSITY ------------------------------------------------------------------------%
udata.density_s = 2650*ones(udata.Nf);       % density of the rock [kg/m3]
udata.density_sf = 2650*ones(udata.Nf_f,1); % density of the rock [kg/m3]

%% THERMAL DIFFUSION ------------------------------------------------------------------------%
udata.lambda_l = 0.5;                           % Thermal conductivity of the fluid [W/(m*K)]
udata.lambda_s = 2.0;                           % Thermal conductivity of the rock [W/(m*K)]
udata.cp_l = 4000;                               % Specific heat capacity of the fluid [J/(kg*K)]
udata.cp_s = 1000;                               % Specific heat capacity of the rock [J/(kg*K)]
udata.ibcD    = zeros(2*sum(udata.Nf),1);           % 1 -> Diffusion on boundary cells
udata.ibcD(round(udata.Nf(2)/2):udata.Nf(2)) = 1;
