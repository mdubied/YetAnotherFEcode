% EXAMPLE: beam meshed with 2D element
clear; 
close all; 
clc
format short g

FORMULATION = 'N1'; % N1/N1t/N0


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
nx = 40;
ny = 2;
[nodes, elements, nset] = mesh_2Drectangle(Lx, Ly, nx, ny);

% nominal mesh
MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);
MeshNominal.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)

% defected mesh
    % arch defect
    xi = 2;                                     % defect amplitude
    yd = Ly * sin(pi/Lx * nodes(:,1));          % y-displacement 
    nodes_defected = nodes + [yd*0 yd]*xi;   	% new nodes
    arc_defect = zeros(numel(nodes),1);         
    arc_defect(2:2:end) = yd;                   % vectorized defect-field
    U = arc_defect;                             % defect basis
MeshDefected = Mesh(nodes_defected);
MeshDefected.create_elements_table(elements,myElementConstructor);
MeshDefected.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
fprintf(' Arc defect: %.1f * beam thickness \n\n', xi)

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

% defected
DefectedAssembly = Assembly(MeshDefected);
Md = DefectedAssembly.mass_matrix();
[Kd,~] = DefectedAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    DefectedAssembly.DATA.K = Kd;
    DefectedAssembly.DATA.M = Md;


%% Eigenmodes                                                       

% Eigenvalue problem_______________________________________________________
n_VMs = 3; % first n_VMs modes with lowest frequency calculated

% Vibration Modes (VM): nominal
Knc = DefectedAssembly.constrain_matrix(Kn);
Mnc = DefectedAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(Knc, Mnc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssembly.unconstrain_vector(VMn);

% Vibration Modes (VM): defected
Kdc = DefectedAssembly.constrain_matrix(Kd);
Mdc = DefectedAssembly.constrain_matrix(Md);
[VMd,om] = eigs(Kdc, Mdc, n_VMs, 'SM');
[f0d,ind] = sort(sqrt(diag(om))/2/pi);
VMd = VMd(:,ind);
for ii = 1:n_VMs
    VMd(:,ii) = VMd(:,ii)/max(sqrt(sum(VMd(:,ii).^2,2)));
end
VMd = NominalAssembly.unconstrain_vector(VMd);

% PLOT (defected)
mod = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure
PlotMesh(nodes_defected, elementPlot, 0);
v1 = reshape(VMd(:,mod), 2, []).';
S = 2*max(nodes_defected(:,2));
PlotFieldonDeformedMesh(nodes_defected, elementPlot, v1, 'factor', S);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0d(mod),3) ' Hz'])
axis on; grid on; box on


%% Modal Derivatives & Defect Sensitivities                         

[MDn, MDnames] = modal_derivatives(NominalAssembly, elements, VMn); % nominal
MDd = modal_derivatives(DefectedAssembly, elements, VMd);           % defected

[DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    FORMULATION); % for DpROM

% PLOT a MD
nwho = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure
v1 = reshape(MDn(:, nwho), 2, []).';
S = Ly/2;
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S, 'component', 'U');
title(['\theta_{' num2str(MDnames(1,nwho)) num2str(MDnames(2,nwho)) '}'])
axis on; grid on; box on


%% (Dp)ROM                                                          

Vn = [VMn MDn DS]; 	% reduced order basis (DpROM)
Vd = [VMd MDd];   	% reduced order basis (ROM-d)

Mnr = Vn'*Mn*Vn; 	% reduced mass matrix (DpROM)
Mdr = Vd'*Md*Vd; 	% reduced mass matrix (ROM-d)

% standard reduced order model (defects in the mesh)
tensors_ROM = reduced_tensors_ROM(DefectedAssembly, elements, Vd);

% parametric formulation for defects
VOLUME = 1;         % integration over defected (1) or nominal volume (0)
tensors_DpROM = reduced_tensors_DpROM(NominalAssembly, elements, Vn, U, ...
    FORMULATION, VOLUME);

% evaluate the defected tensors at xi
[Q2, Q3, Q4, Q3t, Q4t] = DefectedTensors(tensors_DpROM, xi);

% check eigenfrequencies
f0_ROMd = sort(sqrt(eig(Mdr\tensors_ROM.Q2))/2/pi);
f0_DpROM = sort(sqrt(eig(Mnr\Q2))/2/pi);
id = 1 : n_VMs;
f0_ROMd = f0_ROMd(id);
f0_DpROM = f0_DpROM(id);
disp(table(f0n, f0d, f0_ROMd, f0_DpROM))

% add the reduced tensors to the DATA field in Assembly (use this syntax
% for the linear, quadratic and cubic stiffness tensors):
NominalAssembly.DATA.Kr2 = Q2;
NominalAssembly.DATA.Kr3 = Q3;
NominalAssembly.DATA.Kr4 = Q4;









