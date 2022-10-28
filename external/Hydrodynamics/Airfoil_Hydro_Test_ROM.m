% ------------------------------------------------------------------------ 
% Script to test the implementation of hydrodynamic forces on 2D structures
% with TRI3 elements, using a ROM formulation
% 
% Last modified: 28/10/2022, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

%whichModel = 'CUSTOM'; % or "ABAQUS"
whichModel = 'ABAQUS';
%elementType = 'QUAD4';
elementType = 'TRI3';
% elementType = 'QUAD8'; % only QUAD4 is implemented for now

FORMULATION = 'N1'; % N1/N1t/N0
VOLUME = 0;         % integration over defected (1) or nominal volume (0)

USEJULIA = 0;

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 263824;     % Young's modulus [Pa]
rho     = 1070;     % density [kg/m^3]
nu      = 0.499;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
switch elementType
    case 'QUAD4'
        myElementConstructor = @()Quad4Element(thickness, myMaterial);
    case 'QUAD8'
        myElementConstructor = @()Quad8Element(thickness, myMaterial);
    case 'TRI3'
        myElementConstructor = @()Tri3Element(thickness, myMaterial);
end

% MESH_____________________________________________________________________
Lx = 3;
Ly = .2;
nx = 5;
ny = 6;
switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,elementType);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'naca0012TRI';
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
end

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

% boundary conditions of nominal mesh: front and back nodes fixed
frontNode = find_node_2D(0,0,nodes);
Lx = 0.15;
backNode = find_node_2D(Lx,0,nodes);
nset = {frontNode,backNode};
MeshNominal.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)


% ASSEMBLY ________________________________________________________________
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);

% store matrices
NominalAssembly.DATA.K = Kn;
NominalAssembly.DATA.M = Mn;
%f = QuadAssembly.vector('drag_force', u, ud);


%% PLOT MESH WITH NODES AND ELEMENTS
elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 1);


%% EIGENMODES - VIBRATION MODES  

n_VMs = 2;

% EIGENVALUE PROBLEM_______________________________________________________
% Vibration Modes (VM): nominal
Kc = NominalAssembly.constrain_matrix(Kn);
Mc = NominalAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssembly.unconstrain_vector(VMn);

% PLOT ____________________________________________________________________
mod = 1;
elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMn(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0n(mod),3) ' Hz']);

% DAMPING _________________________________________________________________
alfa = 3.1;
beta = 6.3*1e-6;
Dn = alfa*Mn + beta*Kn; % Rayleigh damping
NominalAssembly.DATA.D = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);

%% MODAL DERIVATIVES                      

% nominal
[MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);
% % defected
% MDd = modal_derivatives(DefectedAssembly, elements, VMd);
% 
% % defect sensitivities
% [DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
%     FORMULATION);

%% ROM TENSORS (INTERNAL FORCES)                         
% define reduced order basis
Vn = [VMn MDn];     % reduced order basis (ROM-n)
% V  = [VMn MDn DS]; 	% reduced order basis (DpROM)
% Vd = [VMd MDd];   	% reduced order basis (ROM-d)

% orthonormalize reduction basis
Vn = orth(Vn);	% ROM-n
% V  = orth(V);	% DpROM
% Vd = orth(Vd);	% ROM-d

% standard reduced order model (no defects in the mesh)
tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, Vn, USEJULIA);

% % standard reduced order model (defects in the mesh)
% tensors_ROMd = reduced_tensors_ROM(DefectedAssembly, elements, Vd, USEJULIA);
% tensors_ROMd.xi = xi; % save for which xi ROMd is computed
% 
% % parametric formulation for defects
% tensors_DpROM = reduced_tensors_DpROM(NominalAssembly, elements, ...
%     V, U, FORMULATION, VOLUME, USEJULIA); %compute tensors
% 
% % evaluate the defected tensors at xi
% [Q2, ~, ~, ~, ~, Mxi] = DefectedTensors(tensors_DpROM, xi);

%% Reduced Assembly
ROMn_Assembly = ReducedAssembly(MeshNominal, Vn);

ROMn_Assembly.DATA.M = ROMn_Assembly.mass_matrix();         % reduced mass matrix (ROM-n)
ROMn_Assembly.DATA.C = Vn.'*Dn*Vn;      % reduced damping matrix (ROM-n), using C as needed by the residual function of Newmark integration
ROMn_Assembly.DATA.K = Vn.'*Kn*Vn;    % reduced stiffness matrix (ROM-n)

%% ROM TENSORS - HYDRODYNAMIC FORCES
[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);
vwater = [1;0.1];
rho = 1;
tensors_hydro_ROM = reduced_tensors_hydro_ROM(NominalAssembly, elements, Vn, skinElements, skinElementFaces, vwater, rho);

%% TIME INTEGRATION
F_ext = @(t,q,qd) (tensors_hydro_ROM.Tr1 + tensors_hydro_ROM.Tr2u*q + tensors_hydro_ROM.Tr2udot*qd); % q, qd are reduced order DOFs

% time step for integration
h = 0.05;

% Initial condition: equilibrium
q0 = zeros(size(Vn,2),1);
qd0 = zeros(size(Vn,2),1);
qdd0 = zeros(size(Vn,2),1);

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% Modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMn_Assembly,F_ext);

% Nonlinear Time Integration
tmax = 4.0; 
TI_NL.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL.Solution.u = Vn * TI_NL.Solution.q; % get full order solution

%TI_NL_sol = reshape(TI_NL.Solution.u(:,36), 2, []).';

%% Visualize
%PlotMesh(nodes, elementPlot, 0)
%PlotFieldonDeformedMesh(nodes, elementPlot, TI_NL_sol, 'factor', 100)
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL.Solution.u,'factor',100,'index',1:2,'filename','result_video')



%% BELOW: TESTING OF VISUAL APPEARANCE OF THE FORCES
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET OUTER SURFACE
[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);


%% APPLY FORCE TO SKIN ELEMENTS
%f = BeamAssembly.vector('drag_force', u, ud);
vwater = [1;0];
rho = 1;
%f = AirfoilAssembly.skin_force('force_length_prop_skin_normal', 'weights',
%skinElements, skinElementFaces); % change in assembly -> skin_force needed
%for this function
f = NominalAssembly.skin_force('drag_force', 'weights', skinElements, skinElementFaces, vwater, rho);


%% PLOT FORCES ON MESH
PlotMeshandForce(nodes, elementPlot, 1,30*f);


