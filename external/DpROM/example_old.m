% EXAMPLE: beam meshed with 2D element
clear; 
close all; 
clc


%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .2;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 2;
Ly = .05;
nx = 20;
ny = 2;
[nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny);
    

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
myMesh.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)


% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(myMesh);
M = BeamAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;


%% Eigenmodes                                                       

% Eigenvalue problem_______________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[Phi,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
Phi = Phi(:,ind);
for ii = 1:n_VMs
    Phi(:,ii) = Phi(:,ii)/max(sqrt(sum(Phi(:,ii).^2,2)));
end
Phi = BeamAssembly.unconstrain_vector(Phi);

% PLOT
mod = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(Phi(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', 2*max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])


%% (Dp)ROM                                                          

V = Phi;        % reduced order basis (only VMs for now)

Mr = V'*M*V;    % reduced mass matrix

% standard reduced order model (nominal, no defects)
tensors_ROM = ROM_reduced_tensors(BeamAssembly, elements, V);

% parametric formulation for defects
FORMULATION = 'N1'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)
U = V(:,1:2);       % defect basis
tensors_DpROM = DpROM_reduced_tensors(FORMULATION, VOLUME, ...
    BeamAssembly, elements, V, U);

% evaluate the defected tensors at xi
xi = rand(size(U,2),1)*0;
[Q2, Q3, Q4] = DefectedTensors(tensors_DpROM, xi);

% check eigenfrequencies
f0_DpROM = sort(sqrt(eig(Mr\Q2))/2/pi);

% add the reduced tensors to the DATA field in Assembly (use this syntax
% for the linear, quadratic and cubic stiffness tensors):
BeamAssembly.DATA.Kr2 = Q2;
BeamAssembly.DATA.Kr3 = Q3;
BeamAssembly.DATA.Kr4 = Q4;