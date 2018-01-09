function THERMAID(inputFile,showPlot,attach_pre,attach_post)
% Thermaid
%    Numerical code to solve flow, heat and tracer transport
%    in fractured porous media
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
%  Acknowledgement:  thanks to Manav Tyagi, Brad Mallison and Hadi Hajibeygi
%                    for contributing to the early development of the code. 
%
%  Thermaid: Numerical code to solve flow and tracer transport in
%             fractured porous media using the embdedded discrete fracture model
%
%  THERMAID(inputFile,showPlot,attach_pre,attach_post)
%
%  inputFile    ->  must be a string (optional, default 'inputFile')
%  showPlot     ->  must be a boolean(optional, default '1')
%  attach_pre   ->  must be a string (optional, default 'attach_pre_timestep')
%  attach_post  ->  must be a string (optional, default 'attach_post_timestep')
%
%  THERMAID.m is the main program "driving" the simulation
%  It contains the time loop and calls the neccessary functions.
%
%  To see the structure of the input file, open InputFile.m
% 
%%-------------------------------------------------------------------------%
%
%  Conventions used in the code:
%
%  1) fluxes:  outgoing fluxes are negative; ingoing fluxes are positive
%  2) gravity: positive g -> gravity is directed downwards
%  3) scalar fields: p(x,y), s(x,y), Kx(x,y), Ky(x,y), Q(x,y), QT(x,y), ...
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                             --------------------------
%                             |(1,3)  (2,3) (3,3) (4,3)|
%                          Y  |(1,2)  (2,2) (3,2) (4,2)|
%                             |(1,1)  (2,1) (3,1) (4,1)|
%                             --------------------------
%                                          X
%
%
%  4) boundary conditions (b.c.):
%
%                   ibcs(i) and ibcT(i) define the b.c. type:
%                   1 -> Dirichlet; 0 -> Neumann 
%
%                   Fix(i) and FixT(i) contain the assigned value:
%                   if ibcs(i) = 1 -> Fix(i) is the pressure [Pa]
%                   if ibcs(i) = 0 -> Fix(i) is the flux (per unit 
%                                         of transversal length) [m2/s]
%                   if ibcT(i) = 1 -> FixT(i) is the concentration [kg/m3]
%
%                   the boundary conditions are assigned by a vector that
%                   represents all the cells on the perimeter.
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%
%                              i= 11  i= 12 i= 13 i= 14
%                             --------------------------
%                       i=  3 |(1,3)  (2,3) (3,3) (4,3)| i= 6
%                       i=  2 |(1,2)  (2,2) (3,2) (4,2)| i= 5
%                       i=  1 |(1,1)  (2,1) (3,1) (4,1)| i= 4
%                             --------------------------  
%                              i= 7   i= 8  i=  9 i= 10
%
%  5) vector fields: vx(x,y), vy(x,y), Kx(x,y), Ky(x,y), ...
%                                     [Kx,Ky are diagonal matrix, they can 
%                                              be represented as vector...]
%
%                   vectors are defined at cell interfaces, e.g.,
%                   
%                   vx(i,j) is the velocity between (i,j) and (i+1,j)
%                   vy(i,j) is the velocity between (i,j) and (i,j+1)
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                   vx:
%                             --------------------------
%                           (1,3)  (2,3) (3,3) (4,3) (5,4)
%                        Y  (1,2)  (2,2) (3,2) (4,2) (5,2)
%                           (1,1)  (2,1) (3,1) (4,1) (5,1)
%                             --------------------------
%                                          X
%                   vy:
%                             |(1,4)  (2,4) (3,4) (4,4)|
%                             |(1,3)  (2,3) (3,3) (4,3)|
%                          Y  |(1,2)  (2,2) (3,2) (4,2)|
%                             |(1,1)  (2,1) (3,1) (4,1)|
%                                          X
%
%  6) linear algebra:
%                   for the solution of the linear system, the scalar
%                   fields are transformed into vectors following the
%                   convention:
%
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                             --------------------------
%                             |  9     10    11    12  |
%                          Y  |  5      6     7     8  |
%                             |  1      2     3     4  |
%                             --------------------------
%                                          X
%
%                   thus the vector is
%                                       i= 1 |(1,1)|
%                                       i= 2 |(2,1)|
%                                       i= 3 |(3,1)|
%                                       i= 4 |(4,1)|
%                                       i= 5 |(1,2)|
%                                       i= 6 |(2,2)|
%                                       i= 7 |(3,2)|
%                                       i= 8 |(4,2)|
%                                       i= 9 |(1,3)|
%                                       i=10 |(2,3)|
%                                       i=11 |(3,3)|
%                                       i=12 |(4,3)|
%
%
%-------------------------------------------------------------------------%

% Add path (at beginning of script)
addpath(genpath(pwd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULTS AND GLOBAL PARAMETERS                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global innerIter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         INITIALIZATION                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
udata.ibcp = []; % Empty initialization of internal pressure BC.

if (nargin<1)
    warning('No input file given. Using default: InputFile.')
    inputFile  = 'InputFile'; 
end
if (nargin<2), 
    showPlot = 1; end
if (nargin<3), 
    attach_pre = 'attach_pre_timestep'; end
if (nargin<4),
    attach_post = 'attach_post_timestep'; end

run(inputFile)                                                             % Load Input Data 

x = linspace(0,udata.len(1),2*udata.Nf(1)+1); x = x(2:2:end);              % x grid vector 
y = linspace(0,udata.len(2),2*udata.Nf(2)+1); y = y(2:2:end);              % y grid vector

tNew  = udata.T0  .* ones(udata.Nf);                                       % initialize matrix   temperature solution vector
tNewf = udata.T0f .* ones(udata.Nf_f,1);                                   % initialize fracture temperature solution vector
tN    = [tNew(:); tNewf(:)];

p  = udata.p0  .* ones(udata.Nf);                                          % initialize matrix   pressure saturation solution vector
pf = udata.p0f .* ones(udata.Nf_f,1);                                      % initialize fracture pressure saturation solution vector
pN    = [p(:); pf(:)];

trig_tog = zeros(udata.Nf_f,1);                                            % initialize trigger array for fracture segments


%% Find the intersection of each segment with the grid and compute the
%% Connectivity Index based on the segments lengths and the average distance
%% to the fracture
CI = intersectionsGrid(udata,x,y,XY1);
[row,col] = intersectionsSegments(XY1,XY1);

%% In case of pressure controlled well BC find the fracture segments that
%% lie in the given matrix cell where the BC is applied to
if ~(isempty(udata.ibcp))
    ind_ij = sub2ind(udata.Nf,udata.ibcp(:,1),udata.ibcp(:,2));
    for i = 1:length(ind_ij)
        ind = find(CI(:,1) == ind_ij(i),1);
        if ~(isempty(ind))
            udata.ibcp(i,3) = CI(ind,2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         SIMULATION                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = udata.dt;                                                           % Set the time for the first timestep                                             
dt0 = udata.dt;
i = 0;                                                                     % Count all timesteps
while (time <= udata.timeSim)                                              % Time loop
    
    if showPlot, fprintf(['Simulation time = ' time2str(time) ' / dt = ' time2str(udata.dt) '\n']);end

    i      = i+1;
    
    tOld      = tNew;                                                      % Update saturation of previous iteration (time loop)
    tOldf     = tNewf;
    
    pOld      = p;                                                         % Update saturation of previous iteration (time loop)
    pOldf     = pf;
    
    epsP       = inf;
    epsT       = inf;
    
    run(attach_pre)
        
    innerIter = 1;
    while ((epsP >= udata.tol || epsT >= udata.tol) && innerIter <= udata.maxit) % Inner loop due to saturation dependence of gravity and viscosity
        tIt = tN;                                                          % Update saturation of previous iteration (inner loop)
        pIt = pN;
        
        [Tx,Ty,Tf,Tfm,Tff,DfmT,g,gf,density_l,density_lf]  = initialize(udata,tNew,tNewf,p,pf,CI,row,col);
                                                                           % Initialize (update) the transmissivity based on new temperature
                                                                           % If the viscosity and density are not temperature dependent, this
                                                                           % function could be moved outside of the timeloop
 
                                                                           % If the fractures evolve with time (i.e. fracture creation)
                                                                           % the initializeDFN function also needs to go here.

                                                                           
        [pN]= pressureSystem(udata,pOld,pOldf,Tx,Ty,Tf,g,gf,Q,Tfm,Tff);    % Construct pressure matrix and rhs vectors

        
        p = reshape(pN(1:prod(udata.Nf)),udata.Nf(1),udata.Nf(2));         % Assign solution to matrix solution array 
        pf = pN(prod(udata.Nf)+1:length(pN));                              % Assign solution to fracture solution vector

        [vx,vy,vf,Vfm,Vmf,Vff] = calcVelocity(udata,p,pf,Tx,Ty,Tf,Tfm,Tff,g,gf);  % Calculate Darcy velocities     
        
        if(udata.flagHeatTransport)
            [tN]  = transport_heat_System(udata,tOld,tOldf,vx,vy,vf,Vfm,Vmf,Vff,DfmT,Q,QT,density_l,density_lf); % Solve mass transport equation 

            tNew = tN(1:prod(udata.Nf));                     
            tNew= reshape(tNew,udata.Nf(1),udata.Nf(2));                   % Assign transport solution to matrix solution array 
            tNewf = tN(prod(udata.Nf)+1:length(tN));                       % Assign transport solution to fracture solution array 
        end
        
        if(udata.flagFracStability)
            [stress_data,trig_tog,udata] = calc_frac_stability(udata,pf,trig_tog,tNewf);    
        end
        
        epsP = norm((abs(pN(:) - pIt(:))),inf); 
        epsT = norm((abs(tN(:) - tIt(:))),inf); 

if (innerIter < 2) epsP0 = epsP; end
        
        if showPlot
            fprintf('\t Residual at %d. loop: %d and %d\n', innerIter,epsP, epsT);       
        end
        innerIter = innerIter+1;


        if (innerIter == udata.maxit), error('outer finescale loop did not converge'), end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             Output             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run(attach_post)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     Time Step Control                                             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (innerIter > 25)
        udata.dt = 0.6*udata.dt;
%     elseif (innerIter < 4) % REMOVED THIS FOR RENO SIMULATIONS
%         udata.dt = min(1.2*udata.dt,dt0);
    end
        
    if (time+udata.dt > udata.timeSim) 
    udata.dt = udata.timeSim - time;
    if udata.dt == 0; udata.dt = nan; end
    end

    time=time+udata.dt;
end

% Since THERMAID.m is a function, the variables used are not written into the
% workspace after the simulation has ended. This can be achieved by the
% following if relevant.
ListOfVariables = who;
for k = 1:length(ListOfVariables)
   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
end

try
    close(vidObj);
catch
     %warning('No video object to close')
end
