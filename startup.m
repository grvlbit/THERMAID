%  Startup of THERMAID code - Runs a simple example to verify the build
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

addpath(genpath('../THERMAID'))

load comsol_bench_h.mat

global k_ratio

% Run with Kf/Km = 1e3
k_ratio = 1e3;
THERMAID('Input_ex1',0)

xf = 0:udata.dxf:5;
xf = xf(1:udata.Nf_i(1));

pm_1e3_EDFM = p(:,floor(udata.Nf(2)/2));
pf_1e3_EDFM = pf(udata.Nf_i(1)+1:udata.Nf_f);

pf_interp_1e3 = interp1(x_pf_1e3,pf_1e3,xf');
pm_interp_1e3 = interp1(x_pm_1e3,pm_1e3,x');
pm_interp_1e3(isnan(pm_interp_1e3)) = 0;

RMSE_pf_1e3 = sqrt(mean((pf_1e3_EDFM-pf_interp_1e3).^2))/(max(pf_interp_1e3)-min(pf_interp_1e3));
RMSE_pm_1e3 = sqrt(mean((pm_1e3_EDFM-pm_interp_1e3).^2))/(max(pm_interp_1e3)-min(pm_interp_1e3));

if (RMSE_pf_1e3 > 0.01 || RMSE_pm_1e3 > 0.01)
    warning('The initialization showed bigger error than expected. Please check your Thermaid installation.')
end

clear all, close all
