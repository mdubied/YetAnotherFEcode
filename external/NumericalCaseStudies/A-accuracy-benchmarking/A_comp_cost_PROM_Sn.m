% ------------------------------------------------------------------------ 
% Accuracy benchmarking: computational cost of the PROM model (Sn).
% Used element type: TRI3.
% 
% Last modified: 16/04/2023, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

whichModel = 'ABAQUS';
elementType = 'TRI3';

FORMULATION = 'N1'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 0;

%% PREPARE MODEL                                                    
tic
% DATA ____________________________________________________________________
E       = 2600000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 
thickness = .1;         % [m] out-of-plane thickness

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain

% element
switch elementType
    case 'TRI3'
        myElementConstructor = @()Tri3Element(thickness, myMaterial);
end

% PROM parameters
xi1 = 0.2;

% MESH ____________________________________________________________________

% nominal mesh
switch upper( whichModel )
    case 'ABAQUS'
        filename = 'naca0012TRI3_90Elements';
        [nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
end

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil


%%
% boundary conditions of nominal mesh
nel = size(elements,1);
nset = {};
for el=1:nel   
    for n=1:size(elements,2)
        if  nodes(elements(el,n),1)<Lx*0.15 && ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
    end   
end
for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],1:2,0)
end

% shape variations (only 1 in this example)
% (1) thinner airfoil 
nodes_projected = [nodes(:,1), nodes(:,2)*0];   % projection on x-axis
yDif = nodes_projected(:,2) - nodes(:,2);       % y-difference projection vs nominal
thinAirfoil = zeros(numel(nodes),1);            % create a single long vectors [x1 y1 x2 y2 ...]^T
thinAirfoil(2:2:end) = yDif;                    % fill up all y-positions

% shape variations basis
U = thinAirfoil;    % shape variations basis
xi = xi1;           % shape variations parameter
m = length(xi);     % number of shape variations parameters

% shape-varied mesh 
d = U*xi;                       % displacement field introduced by shape variations
dd = [d(1:2:end) d(2:2:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)


% ASSEMBLY ________________________________________________________________

% nominal
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;



%% DAMPING ________________________________________________________________
alfa = 0.912;
beta = 0.002;

% nominal
Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
NominalAssembly.DATA.C = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);


toc
tic

% %% EIGENMODES - VIBRATION MODES (VMs) _____________________________________
n_VMs = 4;  % number of vibration modes to include in the ROM

% NOMINAL _________________________________________________________________
% eigenvalue problem
Kc = NominalAssembly.constrain_matrix(Kn);
Mc = NominalAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssembly.unconstrain_vector(VMn);




%% MODAL DERIVATIVES (MDs) ________________________________________________                   

% nominal
[MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);

% shape variation sensitivities
[DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    FORMULATION);

%% ROM TENSORS (INTERNAL FORCES) __________________________________________                       
% define reduced order basis

V  = [VMn MDn DS]; 	% reduced order basis (PROM) 

% orthonormalize reduction basis

V  = orth(V);	% PROM
toc

% PROM: parametric formulation for shape variations
tensors_PROM = reduced_tensors_DpROM(NominalAssembly, elements, ...
    V, U, FORMULATION, VOLUME, USEJULIA); %compute tensors

% evaluate the shape-varied tensors at xi
[Q2, Q3, Q4, Q3t, Q4t, M] = DefectedTensors(tensors_PROM, xi);

%% REDUCED ASSEMBLIES _____________________________________________________


% % PROM ____________________________________________________________________
PROM_Assembly = ReducedAssembly(MeshNominal, V);
PROM_Assembly.DATA.M = PROM_Assembly.mass_matrix();     % reduced mass matrix 
PROM_Assembly.DATA.C = V.'*Dn*V;                        % reduced damping matrix 
PROM_Assembly.DATA.K = V.'*Kn*V;                        % reduced stiffness matrix 

%% ROM TENSORS - HYDRODYNAMIC FORCES ______________________________________

[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);
vwater = [0.5;0.5];   % water velocity vector
rho = 997*0.01;
c = 0.2;

% PROM w/o 4th order tensors
FOURTHORDER = 0;
tensors_hydro_PROM = reduced_tensors_hydro_PROM(NominalAssembly, elements, V, U, FOURTHORDER, skinElements, skinElementFaces, vwater, rho, c);


%% TIME INTEGRATION _______________________________________________________

% Parameters' initialization for all models _______________________________
% time step for integration
h1 = 0.01;
tmax = 1.0; 




%% PROM - nominal (defect=0) ______________________________________________
% initial condition: equilibrium
q0 = zeros(size(V,2),1);
qd0 = zeros(size(V,2),1);
qdd0 = zeros(size(V,2),1);

fprintf('solver\n')
tic
% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_PROM.Tr1) + ...
    double(tensors_hydro_PROM.Tru2*q) + double(tensors_hydro_PROM.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_PROM = ImplicitNewmark('timestep',h1,'alpha',0.005,'MaxNRit',60,'RelTol',1e-8);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,PROM_Assembly,F_ext,tensors_hydro_PROM);

% nonlinear Time Integration
TI_NL_PROM.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution
toc
