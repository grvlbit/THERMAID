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

load comsol_bench_th2.mat

global k_ratio

% Run with Kf/Km = 1e5
k_ratio = 1e5;
THERMAID('Input_ex2',0)

%% Post analysis
xf = udata.dxf/2:udata.dxf:abs(XY1(1,2)-XY1(end,3));

T_horz = tNewf(udata.Nf_i(1)+1:udata.Nf_f);
T_vert = tNewf(1:udata.Nf_i(1));

figure; 
p = plot(xf,T_vert,'x',x_ref_T_vert,T_ref_vert,'k-',xf,T_horz,'+',x_ref_T_horz,T_ref_horz,'k--','LineWidth',2.0,'MarkerSize',8);
nummarkers(p,50);
xlabel('x [m]')
ylabel('Temperature [°C]')
legend(p,'EDFM vertical','Comsol vertical','EDFM horizontal','Comsol horizontal')

%% Quantitative analysis
x_interp = xf';
T_interp_horz = interp1(x_ref_T_horz,T_ref_horz,x_interp);

% Remove NaN from interpolation
T_interp_horz(end) = T_interp_horz(end-2);
T_interp_horz(end-1) = T_interp_horz(end-2); 

T_interp_vert = interp1(x_ref_T_vert,T_ref_vert,x_interp);
T_interp_horz(end) = T_interp_horz(end-1);
T_interp_vert(end) = T_interp_vert(end-1);

RMSE_T_vert = sqrt(mean((T_vert-T_interp_vert).^2))/(max(T_interp_vert)-min(T_interp_vert));
RMSE_T_horz = sqrt(mean((T_horz-T_interp_horz).^2))/(max(T_interp_horz)-min(T_interp_horz));

MAE_T_vert = mean(abs(T_vert-T_interp_vert))/(max(T_interp_vert)-min(T_interp_vert));
MAE_T_horz = mean(abs(T_horz-T_interp_horz))/(max(T_interp_horz)-min(T_interp_horz));

