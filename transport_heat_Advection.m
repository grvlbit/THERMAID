function [Up,rhs] = transport_heat_Advection(udata,vx,vy,vf,Vfm,Vmf,Vff,Q,QT,T,Tf,cprho,cprho_f)
%  Calculates the stabilized heat advection term matrix and rhs
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
%  Acknowledgement:  Thanks are due to Hadi Hajibeygi for contributing to
%                    the very early development of the code.
%
%  transport_heat_Advection(udata,vx,vy,vf,Vfm,Vmf,Vff,Q,QT,T,Tf,cprho,cprho_f)
%
%  Input:
%        udata  [struct]        user data
%        vx     (nx+1,ny)       matrix velocity in the x direction
%        vy     (nx,ny+1)       matrix velocity in the y direction
%        vf     (nf,1)          fracture velocity
%        Vfm    (nf,nx*ny)      velocity matrix between fracture and matrix
%        Vmf    (nx*ny,nf)      velocity matrix between matrix and fracture
%        Vff    (nf*nf,nf*nf)   fracture-fracture intersection velocity
%        Q      (nx,ny)         source term of the pressure problem
%        QT     (nx,ny)         source term of the temperature problem
%        T      (nx,ny)         old transport solution on the matrix
%        Tf     (nf,1)          old transport solution in the fractures
%        cprho  (nx,ny)         heat capacity * density in the matrix
%        cprho_f(nf,1)          heat capacity * density in the fractures
%
%  Output:
%        Up    (nx*ny+nf,nx*ny+nf) advection matrix of the transport system
%        rhs   (nx*ny+nf,1)      advection rhs of the transport system

n    = udata.Nf;
Nf   = udata.Nf;
nf   = udata.Nf_f;

use_quick_scheme = 0;                                                                % Set the numerical scheme for the advection
% 0 - Upwind / 1 - QUICK

ix      = (1+sign(vx))/2;                                                  % Velocity indicator for upwinding
iy      = (1+sign(vy))/2;                                                  % for aritmetic mean set ix = iy = 1/2

%-------------------------------------------------------------------------%
% Calculate interface values (harmonic average)
%-------------------------------------------------------------------------%
[cprhox, cprhoy] = calc_interface_values(cprho);
[cprhof] = calc_interface_values_fracture(udata,cprho_f);

cprhovx = cprhox.*vx;
cprhovy = cprhoy.*vy;


if (use_quick_scheme == 1) 

   warning('The second order accurate flux limited QUICK scheme is currently unavailable due to licensing concerns. An alternative is in preparation. Currently the standard upwind method is used as fallback.')
    %-------------------------------------------------------------------------%
    %    Upwind Matrix                                                        %
    %-------------------------------------------------------------------------%
    
    upw(1:Nf(1),:) = -ix(1:Nf(1),:)   .*cprhovx(1:Nf(1),:);
    ups(:,1:Nf(2)) = -iy(:,1:Nf(2))   .*cprhovy(:,1:Nf(2));
    upe(1:Nf(1),:) =  (1-ix(2:Nf(1)+1,:)) .*cprhovx(2:Nf(1)+1,:);
    upn(:,1:Nf(2)) =  (1-iy(:,2:Nf(2)+1)) .*cprhovy(:,2:Nf(2)+1);
    
    Txeast(2:n(1),:)   = upe(1:n(1)-1,:);    Txeast(1,:)     = 0;
    Tynorth(:,2:n(2))  = upn(:,1:n(2)-1);    Tynorth(:,1)    = 0;
    Txwest(1:n(1)-1,:) = upw(2:n(1),:);      Txwest(n(1),:)  = 0;
    Tysouth(:,1:n(2)-1)= ups(:,2:n(2));      Tysouth(:,n(2)) = 0;
    
    upd(1:Nf(1),1:Nf(2)) = ix(2:Nf(1)+1,1:Nf(2))   .*cprhovx(2:Nf(1)+1,1:Nf(2))...
        +iy(1:Nf(1),2:Nf(2)+1)   .*cprhovy(1:Nf(1),2:Nf(2)+1)...
        -(1-ix(1:Nf(1),1:Nf(2))) .*cprhovx(1:Nf(1),1:Nf(2))...
        -(1-iy(1:Nf(1),1:Nf(2))) .*cprhovy(1:Nf(1),1:Nf(2));
    Ds      = [Tysouth(:) Txwest(:) upd(:) Txeast(:) Tynorth(:)];
    
    Up   = spdiags(Ds,[-n(1) -1 0 1 n(1)],n(1)*n(2),n(1)*n(2));    


    
else % This is the first order upwind scheme
    %-------------------------------------------------------------------------%
    %    Upwind Matrix                                                        %
    %-------------------------------------------------------------------------%
    
    upw(1:Nf(1),:) = -ix(1:Nf(1),:)   .*cprhovx(1:Nf(1),:);
    ups(:,1:Nf(2)) = -iy(:,1:Nf(2))   .*cprhovy(:,1:Nf(2));
    upe(1:Nf(1),:) =  (1-ix(2:Nf(1)+1,:)) .*cprhovx(2:Nf(1)+1,:);
    upn(:,1:Nf(2)) =  (1-iy(:,2:Nf(2)+1)) .*cprhovy(:,2:Nf(2)+1);
    
    Txeast(2:n(1),:)   = upe(1:n(1)-1,:);    Txeast(1,:)     = 0;
    Tynorth(:,2:n(2))  = upn(:,1:n(2)-1);    Tynorth(:,1)    = 0;
    Txwest(1:n(1)-1,:) = upw(2:n(1),:);      Txwest(n(1),:)  = 0;
    Tysouth(:,1:n(2)-1)= ups(:,2:n(2));      Tysouth(:,n(2)) = 0;
    
    upd(1:Nf(1),1:Nf(2)) = ix(2:Nf(1)+1,1:Nf(2))   .*cprhovx(2:Nf(1)+1,1:Nf(2))...
        +iy(1:Nf(1),2:Nf(2)+1)   .*cprhovy(1:Nf(1),2:Nf(2)+1)...
        -(1-ix(1:Nf(1),1:Nf(2))) .*cprhovx(1:Nf(1),1:Nf(2))...
        -(1-iy(1:Nf(1),1:Nf(2))) .*cprhovy(1:Nf(1),1:Nf(2));
    Ds      = [Tysouth(:) Txwest(:) upd(:) Txeast(:) Tynorth(:)];
    
    Up   = spdiags(Ds,[-n(1) -1 0 1 n(1)],n(1)*n(2),n(1)*n(2));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     CONSTRUCT  RHS                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhs = sparse(zeros(n(1)*n(2),1));

%------------------------------------------------------%
%                   boundary conditions                %
%------------------------------------------------------%

i1   = 1:n(2);                      i2  = n(2) + (1:n(2));
i3   = 2*n(2) + (1:n(1));           i4  = 2*n(2) + n(1) + (1:n(1));
ic1  = 1:n(1):prod(n);              ic2 = n(1):n(1):prod(n);
ic3  = 1:n(1);                      ic4 = (n(2)-1)*n(1)+1:prod(n);

%---------------------%
% Dirichlet Transport %
%---------------------%

iD1  = ic1';              iD2 = ic2';
iD3  = ic3';              iD4 = ic4';
t1         = upw(1,:)';
t2         = upe(n(1),:)';
t3         = ups(:,1);
t4         = upn(:,n(2));

rhs(iD1,1) = rhs(iD1,1) - t1.*udata.FixT(i1,1);
rhs(iD2,1) = rhs(iD2,1) - t2.*udata.FixT(i2,1);
rhs(iD3,1) = rhs(iD3,1) - t3.*udata.FixT(i3,1);
rhs(iD4,1) = rhs(iD4,1) - t4.*udata.FixT(i4,1);


%clear t1 t2 t3 t4

%------------------------------------------------------%
%    source terms                                      %
%------------------------------------------------------%
iQ  = sparse((1+sign(Q(:)))/2);                                              % Indicator for in or outflow

out = spdiags(cprho(:).*(1-iQ).*Q(:),0,prod(n),prod(n));
in  = cprho(:).*iQ.*Q(:).*QT(:);

rhs = sparse(rhs + in);
Up  = sparse(Up  - out);


%%-------------------------------------------------------------------------%
%%    Fracture                                                             %
%%-------------------------------------------------------------------------%
ifv      = (1+sign(vf))/2;                                             % Velocity indicator for upwinding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     CONSTRUCT  Upwind Matrix                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Upfeast = zeros(nf+udata.N_fractures,1);
Upfwest = zeros(nf+udata.N_fractures,1);
updf   = zeros(nf+udata.N_fractures,1);


A = []; B =[];

ia = 1;
ib = udata.Nf_i(1);
for i = 2:udata.N_fractures
    A = [A ia];
    B = [B ib];
    
    ia = ib + 2;
    ib = ia+ udata.Nf_i(i) -1;
    
end
A = [A ia];
B = [B ib];

for i = 1:udata.N_fractures
    ia = A(i);
    ib = B(i);
    
    upwf = -ifv(ia:ib)       .*cprhof(ia:ib).*vf(ia:ib)  ;
    upef =  (1-ifv(ia+1:ib+1)).*cprhof(ia+1:ib+1).*vf(ia+1:ib+1);
    
    updf(ia:ib)      = ifv(ia:ib)   .*cprhof(ia:ib).*vf(ia:ib) ...
        -(1-ifv(ia+1:ib+1)) .*cprhof(ia+1:ib+1).*vf(ia+1:ib+1);
    
    Upfeast(ia+1:ib) = upef(1:end-1);
    Upfwest(ia:ib-1) = upwf(2:end);
end

Dsf  = [Upfwest updf Upfeast];
Upf  = spdiags(Dsf,[-1,0,1],nf,nf);
rhsf = sparse(zeros(nf,1));


%------------------------%
%  internal pressure BC  %
%------------------------%
% if ~(isempty(udata.ibcp))
%     ind = sub2ind(n,udata.ibcp(:,1),udata.ibcp(:,2));
%
%     Up(ind,:) = 0;
%     Up(ind,ind) = 1;
%     rhs(ind) = udata.ibcp(:,5);
% end

if ~(isempty(udata.ibcp))
    for i = 1:length(udata.ibcp(:,1))
        ind = udata.ibcp(i,3);
        if (ind > 0)
            Upf(ind,:) = 0;
            Upf(ind,ind) = 1;
            rhsf(ind) = udata.ibcp(i,5);
        end
    end
end

%-------------------------------------------------------------------------%
%    Matrix-Fracture                                                      %
%-------------------------------------------------------------------------%
cprhoVmf = bsxfun(@times, cprho_f', Vmf);
Upmf = -min(cprhoVmf,0);
Ds =  sum(max(cprhoVmf,0),2);

[m,n]=size(Up);
B = spdiags(Ds,0,m,n);                                                     % Diagonal contribution of Amf to the main diagonal of A

%-------------------------------------------------------------------------%
%    Fracture-Matrix                                                      %
%-------------------------------------------------------------------------%
cprhoVfm = bsxfun(@times, cprho_f, Vfm);
Upfm = -min(cprhoVfm,0);
DsT =  sum(max(-cprhoVfm,0),2);

[m,n]=size(Upf);

%-------------------------------------------------------------------------%
%    Fracture-Fracture                                                    %
%-------------------------------------------------------------------------%
cprhoVff = bsxfun(@times, cprho_f, Vff);
DsTT = -sum(max(cprhoVff,0),2);
Tff = max(cprhoVff,0);
DsT = DsT + DsTT;
Bf = spdiags(DsT,0,m,n);                                                   % Diagonal contribution of Afm and Aff to the main diagonal of A

%-------------------------------------------------------------%
%        merging matrix and fracture matrix and rhs           %
%-------------------------------------------------------------%

% Fracture-only solution
% Up = blkdiag(Up,Upf);

% Fracture-matrix solution
Up = Up + B;
Upf = Upf + Bf + Tff;
Up = [Up -Upmf; -Upfm Upf];

rhs = vertcat(rhs,rhsf);                                                   % Concenate the RHS vectors of matrix and fracture
