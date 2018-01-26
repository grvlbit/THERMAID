function test_suite=test_pressureSystem
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end


function test_validate_gaussian
%  Validates the time-dependent pressure solver against the analytical
%  solution for diffusion of a gaussian
evalin('base','clear  calcVelocity initialize pressureSystem calc_interface_values_fracture');

udata.len     = [10 10];                            % physical length of the domain in x and y direction [m]
udata.Nf      = [100 100];  	                    % number of cells in x and y direction
udata.dx      = udata.len./udata.Nf;                % cell length [m]

frac_cross;

udata.ibcs = zeros(2*sum(udata.Nf),1);              % type 0:Neumann(N); 1:Dirichlet(D)
udata.Fix  = zeros(2*sum(udata.Nf),1);              % value N [m2/s] (inflow>0); D [Pa]   

Q       = zeros(udata.Nf);                          % source term [m2/s]; inflow positive
udata.gravity   = 0;                                % gravity acceleration in y [m/s2]

udata.K       = ones(udata.Nf(1),udata.Nf(2))*1; 	% permeability field [m2]
udata.K_f     = ones(udata.Nf_f,1)*1e-18; 	        % fracture permeability field [m2]

udata.b0      = 1*ones(udata.Nf_f,1);               % fracture aperture field [m]

udata.phi     = ones(udata.Nf(1),udata.Nf(2));	    % porosity field
udata.phi_f   = 0*ones(udata.Nf_f,1);               % fracture porosity field

udata.const_density = 1;
if(udata.const_density)
    udata.density_l  = 1*ones(udata.Nf);            % density of the rock [kg/m3]
    udata.density_lf = 1*ones(udata.Nf_f,1);        % density of the rock [kg/m3]
end

udata.const_viscosity = 1;
if(udata.const_viscosity)
    udata.viscosity   = 1*ones(udata.Nf);           % density of the rock [kg/m3]
    udata.viscosity_f = 1*ones(udata.Nf_f,1);       % density of the rock [kg/m3]
end

udata.flagIncompressible = 0;                       % If 1 - incompressible fluids are used
udata.compressibility_l = 1;                        % Bulk Modulus of the fluid [Pa]
udata.compressibility_s = 0;                        % Bulk Modulus of the rock  [Pa]

udata.density_s = 1*ones(udata.Nf);                 % density of the rock [kg/m3]
udata.density_sf = 1*ones(udata.Nf_f,1);            % density of the rock [kg/m3]

udata.lambda_l = 0;                                 % Thermal conductivity of the fluid [W/(m*K)]
udata.lambda_s = 0;                                 % Thermal conductivity of the rock [W/(m*K)]
udata.cp_l = 0;                                     % Specific heat capacity of the fluid [J/(kg*K)]
udata.cp_s = 0;                                     % Specific heat capacity of the rock [J/(kg*K)]

udata.ibcp = []; % Empty initialization of internal pressure BC.

tNew  = ones(udata.Nf);                             % initialize matrix   temperature solution vector
tNewf = ones(udata.Nf_f,1);                         % initialize fracture temperature solution vector

x = linspace(0,udata.len(1),2*udata.Nf(1)+1); x = x(2:2:end);                          % x grid vector (used in plots)
y = linspace(0,udata.len(2),2*udata.Nf(2)+1); y = y(2:2:end);                          % y grid vector (used in plots)

CI = intersectionsGrid(udata,x,y,XY1);
[row,col] = intersectionsSegments(XY1,XY1);

% Numerical input
W       =   udata.len(1);
H       =   udata.len(2);
nx      =   udata.Nf(1);
ny      =   udata.Nf(2);
sigma   =   1.0;
A       =   1.5;

% Initial gaussian wave input
[x,y] =   meshgrid([-W/2:W/(nx-1):W/2],[-H/2:H/(ny-1):H/2]);
x = x';
y = y';

p0      =   A*exp(-(x.^2 + y.^2)./(sigma^2));
p0f     =   zeros(size(tNewf));

p_old   =   p0;   
p_oldf  =   p0f;
p       =   p0;   
pf      =   p0f;
pN      =   [p_old(:); p_oldf(:)];

[Tx,Ty,Tf,Tfm,Tff,~,g,gf,~,~]  = initialize(udata,tNew,tNewf,p,pf,CI,row,col);

tol = 1e-5;
maxit = 100;

count = 1;
dts = [];
for dt=0.02:0.02:0.1
    udata.dt = dt;
    
    p =   p0;
    pf=   p0f;
    for time    =   0:dt:1
        p_old      =   p;
        p_oldf     =   pf;
        epsP       = inf;
        innerIter = 1;
        while (epsP >= tol && innerIter <= maxit)                            % Inner loop due to saturation dependence of gravity and viscosity
            pIt = pN;

            P_anal     =   A./(1+4*time/sigma^2).*exp(-(x.^2 + y.^2)./(4*time+sigma^2));      
            [pN] = pressureSystem(udata,p_old,p_oldf,Tx,Ty,Tf,g,gf,Q, Tfm,Tff);
            
            p = reshape(pN(1:prod(udata.Nf)),udata.Nf(1),udata.Nf(2));                            % Assign solution to matrix solution array 
            pf = pN(prod(udata.Nf)+1:length(pN));                                      % Assign solution to fracture solution vector
        
            epsP = norm((abs(pN(:) - pIt(:))),inf); 

            innerIter = innerIter+1;
        end

%         figure(111)
%         subplot(221)
%         pcolor(x,y,p), shading interp, colorbar,title('P calculated'); title(time)
%         subplot(222)
%         pcolor(x,y,P_anal), shading interp, colorbar,title('P analytical');
%         subplot(223)
%         pcolor(x,y,p-P_anal), shading interp, colorbar,title('Error');
%         drawnow
    end
    P_anal     =   A./(1+4*time/sigma^2).*exp(-(x.^2 + y.^2)./(4*time+sigma^2));      
    perror(count) = sqrt(mean((p-P_anal).^2))/(max(P_anal)-min(P_anal));

    count = count +1;
  
    dts = [dts dt];
    
%     subplot(224)
%     plot(dts,100*perror/A,'x-')
%     xlabel('time step [s]')
%     ylabel('Error [%]')
end

assertTrue(perror(1)<=0.0087);
assertTrue(perror(2)<=0.0115);
assertTrue(perror(3)<=0.0145);
assertTrue(perror(4)<=0.0174);
assertTrue(perror(5)<=0.0201);
end
