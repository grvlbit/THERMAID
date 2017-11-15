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

addpath ../

load comsol_bench_h.mat

global k_ratio

% Run with Kf/Km = 1e3
k_ratio = 1e3;
THERMAID('Input_ex1',0)

pm_1e3_EDFM = p(:,floor(udata.Nf(2)/2));
pf_1e3_EDFM = pf(udata.Nf_i(1)+1:udata.Nf_f);

% Run with Kf/Km = 1e5
k_ratio = 1e5;
THERMAID('Input_ex1',0)

pm_1e5_EDFM = p(:,floor(udata.Nf(2)/2));
pf_1e5_EDFM = pf(udata.Nf_i(1)+1:udata.Nf_f);


xf = 0:udata.dxf:5;
xf = xf(1:udata.Nf_i(1));
figure; 
p = plot(xf+udata.dxf/2,pf_1e3_EDFM,'x',x_pf_1e3,pf_1e3,'k-',xf+udata.dxf/2,pf_1e5_EDFM,'+',x_pf_1e5,pf_1e5,'k--','LineWidth',2.0,'MarkerSize',8);
nummarkers(p,50);
xlabel('x [m]')
ylabel('Fracture Pressure [Pa]')
legend(p,'Thermaid 10^3','Comsol 10^3','Thermaid 10^5','Comsol 10^5')

figure; 
p = plot(x,pm_1e3_EDFM,'x',x_pm_1e3,pm_1e3,'k-',x,pm_1e5_EDFM,'+',x_pm_1e5,pm_1e5,'k--','LineWidth',2.0,'MarkerSize',8);
nummarkers(p,50);
xlabel('x [m]')
ylabel('Matrix Pressure [Pa]')
legend(p,'Thermaid 10^3','Comsol 10^3','Thermaid 10^5','Comsol 10^5')

%% Quantitative analysis

pf_interp_1e3 = interp1(x_pf_1e3,pf_1e3,xf');
pm_interp_1e3 = interp1(x_pm_1e3,pm_1e3,x');
pm_interp_1e3(isnan(pm_interp_1e3)) = 0;

pf_interp_1e5 = interp1(x_pf_1e5,pf_1e5,xf');
pm_interp_1e5 = interp1(x_pm_1e5,pm_1e5,x');
pm_interp_1e5(isnan(pm_interp_1e5)) = 0;

RMSE_pf_1e3 = sqrt(mean((pf_1e3_EDFM-pf_interp_1e3).^2))/(max(pf_interp_1e3)-min(pf_interp_1e3))
RMSE_pm_1e3 = sqrt(mean((pm_1e3_EDFM-pm_interp_1e3).^2))/(max(pm_interp_1e3)-min(pm_interp_1e3))

RMSE_pf_1e5 = sqrt(mean((pf_1e5_EDFM-pf_interp_1e5).^2))/(max(pf_interp_1e5)-min(pf_interp_1e5))
RMSE_pm_1e5 = sqrt(mean((pm_1e5_EDFM-pm_interp_1e5).^2))/(max(pm_interp_1e5)-min(pm_interp_1e5))

MAE_pf_1e3 = mean(abs(pf_1e3_EDFM-pf_interp_1e3))/(max(pf_interp_1e3)-min(pf_interp_1e3))
MAE_pm_1e3 = mean(abs(pm_1e3_EDFM-pm_interp_1e3))/(max(pm_interp_1e3)-min(pm_interp_1e3))

MAE_pf_1e5 = mean(abs(pf_1e5_EDFM-pf_interp_1e5))/(max(pf_interp_1e5)-min(pf_interp_1e5))
MAE_pm_1e5 = mean(abs(pm_1e5_EDFM-pm_interp_1e5))/(max(pm_interp_1e5)-min(pm_interp_1e5))


rmpath ../

