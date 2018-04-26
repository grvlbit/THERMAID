function [stress_data,trig_tog,udata] = calc_frac_stability(udata,pf,trig_tog,tNewf)
%  Calculate fracture stability based on slip tendency
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
%  Authors: Gunnar Jansen, University of Neuchatel, 2017
%           Benoît Valley, University of Neuchatel, 2017
%
%  [stress_data,trig_tog,udata] = calc_frac_stability(udata,pf,trig_tog,tNewf)
%
%  Input:
%        udata  [struct]        user data
%        pf     (nf,1)          pressure solution in the fractures
%        trig_tog(nf,1)         keep track of 'failed' H-fracture segments
%        tNewf  (nf,1)          temperature solution in the fractures

%  Output:
%        stress_data(nf,4)      array containing the normal and shear
%                               stress on each fracture segment as well as 
%                               the shear and dilation tendency in columns
%        trig_tog(nf,1)         keep track of 'failed' TH-fracture segments
%        udata  [struct]        user data
%
%  Note: 
%        In order to plot the slip tendency on a mohr diagram use the
%        following commands:
%
%        fig4 = plot_mohr([S1 S2 S3],'Fault stability');
%        hold on
%        scatter(data(:,1)+pf,data(:,2),50,shear_tendency,'filled')


% For the 'standard' 2D case where the principal stresses are perpendicular
% to the sides of the domain the following parameters should be used:

% S1 vertical - S3 horizontal -> plunge = 90 -> trend = [0 0] 
%                             -> frac_az = 0
%                             -> frac_dip = min(abs(udata.frac_angle),180-abs(udata.frac_angle))

% S3 vertical - S1 horizontal -> plunge = 0  -> trend = [0 0]
%                             -> frac_az = 0
%                             -> frac_dip = min(abs(udata.frac_angle),180-abs(udata.frac_angle))

if (size(udata.stress_trend,1) > 1)
    TR1 = udata.stress_trend(:,1);
    TR3 = udata.stress_trend(:,2);
    PL1 = udata.stress_plunge;
else
    TR1 = udata.stress_trend(1)*ones(udata.Nf_f,1);
    TR3 = udata.stress_trend(2)*ones(udata.Nf_f,1);
    PL1 = udata.stress_plunge*ones(udata.Nf_f,1);
end

S1 = udata.sigma_1;
S2 = udata.sigma_2;
S3 = udata.sigma_3;


if (udata.use_thermal_stress)
    % --------------------------Thermal stress-------------------------------------%
    E = 2*udata.shear_modulus*(1+udata.poisson_ratio);

    deltaT = tNewf-udata.T0f.*ones(size(tNewf));
    % Small unexpected temperature changes can be introduced by the QUICK scheme.
    % In this case an expression like deltaT = min(deltaT,0.0); might be useful here.
    St = E/ (1 - 2 * udata.poisson_ratio)*udata.therm_exp_coeff*deltaT;
end

Xg = [udata.frac_az' udata.frac_dip'];     % az and dip of all fracture segments

shear_tendency = zeros(udata.Nf_f,1);
dilation_tendency = zeros(udata.Nf_f,1);

data = zeros(udata.Nf_f,2);

for i=1:length(Xg(:,1))
    
    if (udata.use_thermal_stress)
        [sigma_n, tau]=shearstress_on_fracture([S1(i)-pf(i)+St(i) 0 0; 0 S2(i)-pf(i)+St(i) 0; 0 0 S3(i)-pf(i)+St(i)],TR1(i),PL1(i),TR3(i),Xg(i,1),Xg(i,2));
    else
        [sigma_n, tau]=shearstress_on_fracture([S1(i)-pf(i) 0 0; 0 S2(i)-pf(i) 0; 0 0 S3(i)-pf(i)],TR1(i),PL1(i),TR3(i),Xg(i,1),Xg(i,2));
    end
    
    shear_tendency(i)=tau/sigma_n;
    
    dilation_tendency(i)=(S1(i)-sigma_n)/(S1(i)-S3(i));
    
    data(i,1) = sigma_n; data(i,2) = tau;
    
end

%-----------Evaluate if a point has failed and make sure it has not failed before------------------%
for k = 1:udata.Nf_f;
    if ((shear_tendency(k)>udata.friction_coeff(k)||shear_tendency(k)<-udata.friction_coeff(k)) && trig_tog(k) == 0)
        trig_tog(k) = 1;
        
        % This is the stepwise increase in permeability for fracture segments eligible for slip.
        % Alternatively, the aperture can be modified instead of the permeability.
        udata.K_f(k) = udata.K_enh*udata.K_f(k);
        % udata.b0(k) = udata.K_enh*udata.b0(k);       
    end;
end;

stress_data = [data shear_tendency dilation_tendency];

end




function [n, s]=shearstress_on_fracture(T,tX1,pX1,tX3,az,dip)
% SHEARSTRESS_ON_FRACTURE compute the shear and normal stress on a fracture
% [n s]=shearstress_on_fracture(T,tX1,pX1,tX3,az,dip)
% INPUT:
%   T   : stress tensor
%   tX1 : trend of X1
%   pX1 : plunge of X1
%   tX3 : trend of X3
%   az  : azimuth of the fracture
%   dip : dip of the fracture
% OUTPUT:
%   n : normal stress on the fracture
%   s : shear stress on the fracture

% calculate conversion matrix
Rlmn = DirCosAxes(tX1*pi/180,pX1*pi/180,tX3*pi/180);

% compute stresses in the geographical referentiel  [x to north, y to east, z down]
Txyz=Rlmn'*T*Rlmn;

% calculate conversion matrice for fracture local referentiel
a=az; b=-dip; c=0;
Ruvw=[  cosd(a)*cosd(b)                             sind(a)*cosd(b)                             -sind(b)
    cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c)     sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c)     cosd(b)*sind(c)
    cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c)     sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c)     cosd(b)*cosd(c)];

% stresses in local fracture referentiel  (total stress)
Tuvw=Ruvw*Txyz*Ruvw';
s=sqrt(Tuvw(1,3)^2+Tuvw(2,3)^2);
n=Tuvw(3,3);
end



function dC = DirCosAxes(tX1,pX1,tX3)
%DirCosAxes calculates the direction cosines of a right handed, orthogonal
%X1,X2,X3 cartesian coordinate system of any orientation with respect to
%North-East-Down
%
%   USE: dC = DirCosAxes(tX1,pX1,tX3)
%
%   tX1 = trend of X1
%   pX1 = plunge of X1
%   tX3 = trend of X3
%   dC = 3 x 3 matrix containing the direction cosines of X1 (row 1),
%        X2 (row 2), and X3 (row 3)
%
%   Note: Input angles should be in radians
%
%   DirCosAxes uses function SphToCart
%
%MATLAB script written by Nestor Cardozo for the book Structural
%Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
%this script, please cite this as "Cardozo in Allmendinger et al. (2011)"

%Some constants
east = pi/2.0;
west = 1.5*pi;
tol = 1e-6; % new: tolerance for the difference between tX1 and tX3

%Initialize matrix of direction cosines
dC = zeros(3,3);

%Direction cosines of X1
[dC(1,1),dC(1,2),dC(1,3)] = SphToCart(tX1,pX1,0);

%Calculate plunge of axis 3
%If axis 1 is horizontal
if pX1 == 0.0
    dt = abs(tX1-tX3); % new
    if abs(dt - east) < tol || abs(dt - west) < tol % new
        pX3 = 0.0;
    else
        pX3 = east;
    end
    %Else
else
    %From Equation 2.14 and with theta equal to 90 degrees
    pX3 = atan(-(dC(1,1)*cos(tX3)+dC(1,2)*sin(tX3))/dC(1,3));
end

%Direction cosines of X3
[dC(3,1),dC(3,2),dC(3,3)] = SphToCart(tX3,pX3,0);

%Compute direction cosines of X2 by the cross product of X3 and X1
dC(2,1) = dC(3,2)*dC(1,3) - dC(3,3)*dC(1,2);
dC(2,2) = dC(3,3)*dC(1,1) - dC(3,1)*dC(1,3);
dC(2,3) = dC(3,1)*dC(1,2) - dC(3,2)*dC(1,1);
% Convert X2 to a unit vector
r = sqrt(dC(2,1)*dC(2,1)+dC(2,2)*dC(2,2)+dC(2,3)*dC(2,3));
for i = 1:3
    dC(2,i) = dC(2,i)/r;
end

end




function [cn,ce,cd] = SphToCart(trd,plg,k)
%SphToCart converts from spherical to cartesian coordinates
%
%[cn,ce,cd] = SphToCart(trd,plg,k) returns the north (cn),
% east (ce), and down (cd) direction cosines of a line.
%
% k is an integer to tell whether the trend and plunge of a line
% (k = 0) or strike and dip of a plane in right hand rule
% (k = 1) are being sent in the trd and plg slots. In this
% last case, the direction cosines of the pole to the plane
% are returned
%
% NOTE: Angles should be entered in radians

%If line (see Table 2.1)
if k == 0
    cd = sin(plg);
    ce = cos(plg) * sin(trd);
    cn = cos(plg) * cos(trd);
    %Else pole to plane (see Table 2.1)
elseif k == 1
    cd = cos(plg);
    ce = -sin(plg) * cos(trd);
    cn = sin(plg) * sin(trd);
end
end
