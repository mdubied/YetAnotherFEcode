% ------------------------------------------------------------------------ 
% Accuracy benchmarking: computational cost of the FOM-t model.
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
%% ROM TENSORS - HYDRODYNAMIC FORCES ______________________________________
tic
[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);
vwater = [0.5;0.5];   % water velocity vector
rho = 997*0.01;
c = 0.2;

%% FOM TENSORS - HYDRODYNAMIC FORCES (optional) ___________________________

FOM = 1;    % FOM=1 for assembling the hydrodynamic tensors at the assembly level (FOM)

if FOM == 1
    tensors_hydro_FOM = unreduced_tensors_hydro_FOM(NominalAssembly, elements, skinElements, skinElementFaces, vwater, rho, c);
end

toc
%% TIME INTEGRATION _______________________________________________________

% Parameters' initialization for all models _______________________________
% time step for integration
h1 = 0.01;
tmax = 1.0; 

%% FOM-II (nominal) _______________________________________________________
% initial condition: equilibrium
fprintf('solver \n')
nUncDOFs = size(MeshNominal.EBC.unconstrainedDOFs,2);
q0 = zeros(nUncDOFs,1);
qd0 = zeros(nUncDOFs,1);
qdd0 = zeros(nUncDOFs,1);
tic
% hydrodynamic forces
T1 = NominalAssembly.constrain_vector(double(tensors_hydro_FOM.T1));
Tu2 = NominalAssembly.constrain_matrix(double(tensors_hydro_FOM.Tu2));
Tudot2 = NominalAssembly.constrain_matrix(double(tensors_hydro_FOM.Tudot2));
Tuu3 = tensor(NominalAssembly.constrain_tensor(double(tensors_hydro_FOM.Tuu3)));
Tuudot3 = tensor(NominalAssembly.constrain_tensor(double(tensors_hydro_FOM.Tuudot3)));
Tudotudot3 = tensor(NominalAssembly.constrain_tensor(double(tensors_hydro_FOM.Tudotudot3)));
F_ext = @(t,q,qd) (T1 + Tu2*q + Tudot2*qd + double(ttv(ttv(Tuu3,q,3), q,2)) + ...
                    double(ttv(ttv(Tuudot3,qd,3),q,2))+ double(ttv(ttv(Tudotudot3,qd,3), qd,2)));

% instantiate object for nonlinear time integration
TI_NL_FOMn = ImplicitNewmark('timestep',h1,'alpha',0.005,'MaxNRit',200,'RelTol',1e-6);

% modal nonlinear Residual evaluation function handle
Residual_NL = @(q,qd,qdd,t)residual_nonlinear_hydro(q,qd,qdd,t,NominalAssembly,F_ext,Tu2,Tudot2,Tuu3,Tuudot3,Tudotudot3);

% nonlinear Time Integration
TI_NL_FOMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
TI_NL_FOMn.Solution.u = zeros(NominalAssembly.Mesh.nDOFs,size(TI_NL_FOMn.Solution.q,2));
for t=1:size(TI_NL_FOMn.Solution.q,2)
    TI_NL_FOMn.Solution.u(:,t) = NominalAssembly.unconstrain_vector(TI_NL_FOMn.Solution.q(:,t));
end
toc




