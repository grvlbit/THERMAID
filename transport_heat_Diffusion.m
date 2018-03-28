function [Di,ri] = transport_heat_Diffusion(udata,phix,phiy,phif,Bf,Dfm)
%  Creates the diffusion matrix  
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
%  transport_heat_Diffusion(udata,phix,phiy,phif,Kf,Dfm)
%
%  Input: 
%        udata  [struct]        user data
%        phix   (nx+1,ny)       matrix porosity at the interface in
%                               x-direction
%        phiy   (nx,ny+1)       matrix porosity at the interface in
%                               y-direction
%        phif   (nf+Nf,1)       fracture porosity at the interface
%        Kf     (nf+Nf,1)       fracture permeability at the interface
%        Dfm    (nf,nx*ny)      fracture-matrix diffusivity
%
%  Output: 
%        Di     (nx*ny+nf,nx*ny+nf) diffusion matrix of the transport system
%        ri     (nx*ny+nf,1)    diffusion rhs of the transport system

Nf = udata.Nf;
Nf_f = udata.Nf_f;
Nf_i = udata.Nf_i;
dx = udata.dx;

Dx = dx(2).*(phix.*udata.lambda_l+(1-phix).*udata.lambda_s)./dx(1);
Dy = dx(1).*(phiy.*udata.lambda_l+(1-phiy).*udata.lambda_s)./dx(2);
Df = Bf.*(phif.*udata.lambda_l+(1-phif).*udata.lambda_s)./udata.dxf;

Tdeast  = sparse(Nf(1),Nf(2));
Tdwest  = sparse(Nf(1),Nf(2));
Tdsouth = sparse(Nf(1),Nf(2));
Tdnorth = sparse(Nf(1),Nf(2));

N_fractures = length(Nf_i); 
[frac_left_flux_mat, frac_right_flux_mat] = calc_frac_flux_mat(N_fractures,Nf_f, Nf_i);

%-------------------------------------------------------------%
%        preparation & setting of matrix Di                   %
%-------------------------------------------------------------%
Tdeast (2:Nf(1),:)   = Dx(2:Nf(1),:); Tdeast(1,:)      = 0;
Tdnorth(:,2:Nf(2))   = Dy(:,2:Nf(2)); Tdnorth(:,1)     = 0;
Tdwest (1:Nf(1)-1,:) = Dx(2:Nf(1),:); Tdwest(Nf(1),:)  = 0;
Tdsouth(:,1:Nf(2)-1) = Dy(:,2:Nf(2)); Tdsouth(:,Nf(2)) = 0;

Ds      = [Tdsouth(:) Tdwest(:) zeros(prod(Nf),1) Tdeast(:) Tdnorth(:)];
Ds(:,3) = -sum(Ds,2);
Di      = spdiags(-Ds,[-Nf(1),-1,0,1,Nf(1)],Nf(1)*Nf(2),Nf(1)*Nf(2));

Dsf      = [frac_right_flux_mat*Df zeros(Nf_f,1) frac_left_flux_mat*Df];
Dsf(:,2) = -sum(Dsf,2);
Di_f       = spdiags(-Dsf,[-1,0,1],Nf_f,Nf_f);

%-------------------------------------------------------------%
%        boundary conditions (matrix only)                    %
%-------------------------------------------------------------%

ri   = sparse(prod(Nf),1);
rif  = zeros(Nf_f,1);

i1   = 1:Nf(2);                      i2  = Nf(2) + (1:Nf(2)); 
i3   = 2*Nf(2) + (1:Nf(1));          i4  = 2*Nf(2) + Nf(1) + (1:Nf(1));

ic1  = 1:Nf(1):prod(Nf);             ic2 = Nf(1):Nf(1):prod(Nf); 
ic3  = 1:Nf(1);                      ic4 = (Nf(2)-1)*Nf(1)+1:prod(Nf);

t1(1:Nf(2),1)  = Dx(1,1:Nf(2))'.*udata.ibcD(i1);   
t2(1:Nf(2),1)  = Dx(Nf(1)+1,1:Nf(2))'.*udata.ibcD(i2);
t3(1:Nf(1),1)  = Dy(1:Nf(1),1) .*udata.ibcD(i3);   
t4(1:Nf(1),1)  = Dy(1:Nf(1),Nf(2)+1) .*udata.ibcD(i4);   

iD  = [ic1';ic2';ic3';ic4'];
tD  = [t1;t2;t3;t4];
Di  = Di + sparse(iD,iD,tD,prod(Nf),prod(Nf));                     

ri(ic1,1) = ri(ic1,1) + t1.*udata.FixT(i1,1);
ri(ic2,1) = ri(ic2,1) + t2.*udata.FixT(i2,1);
ri(ic3,1) = ri(ic3,1) + t3.*udata.FixT(i3,1);
ri(ic4,1) = ri(ic4,1) + t4.*udata.FixT(i4,1);


%------------------------%
%  internal pressure BC  %
%------------------------%
% if ~(isempty(udata.ibcp))
%     ind = sub2ind(Nf,udata.ibcp(:,1),udata.ibcp(:,2)); 
% 
%     Di(ind) = Di(ind) + Dx(ind);
%     ri(ind) = ri(ind) + Dx(ind).*udata.ibcp(:,5);
% end


if ~(isempty(udata.ibcp))
    for i = 1:length(udata.ibcp(:,1))
        ind = udata.ibcp(i,3);
        if (ind > 0 && udata.ibcp(i,5) > 0)
            Di_f(ind,:) = 0;
            Di_f(ind,ind) = 1;
            rif(ind) = udata.ibcp(i,5);
        end 
    end
end

%-------------------------------------------------------------%
%        merging matrix and fracture matrix and rhs           %
%-------------------------------------------------------------%

Dmf = Dfm';
Ds = sum(Dmf,2);
DsT = sum(Dfm,2);
[m,n]=size(Di);
B = spdiags(Ds,0,m,n);                                                     % Diagonal contribution of Amf to the main diagonal of A 
[m,n]=size(Di_f);
Bf = spdiags(DsT,0,m,n);   

%% Matrix-Fracture solution
Di = Di + B;
Di_f = Di_f + Bf;
Di = [Di -Dmf; -Dfm Di_f];  

%% Fracture only solution
% Di_f = spdiags(zeros(Nf_f,1),0,Nf_f,Nf_f);
% Di = blkdiag(Di, Di_f);
 
ri = vertcat(ri,rif);
