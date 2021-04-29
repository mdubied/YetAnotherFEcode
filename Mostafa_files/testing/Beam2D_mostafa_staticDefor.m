% EXAMPLE: beam meshed with 3D element
clear; 
close all; 
clc

whichModel = 'CUSTOM'; % or "ABAQUS"


%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 3;
Ly = .3;
nx = 30;
ny = 3;
switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'Job-BeamQuad';
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
end

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
switch upper( whichModel )
    case 'CUSTOM'
        myMesh.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
    case 'ABAQUS'
        myMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
end

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(myMesh);
M = BeamAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;
C= 0*M+0*K;

%% EXAMPLE 1                                                        

% Eigenvalue problem_______________________________________________________
n_VMs = 2; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = BeamAssembly.unconstrain_vector(V0);

% PLOT
mod = 2;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 1);
v1 = reshape(V0(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', Ly*1.1);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])


%% EXAMPLE 2                                                        

% Define external force:
% % Body force
% Pressure = 1e6;
% F = Pressure*BeamAssembly.uniform_body_force();

% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(2)) = 2e6;

u_lin = BeamAssembly.solve_system(K, F);
ULIN = reshape(u_lin,2,[]).';	% Linear response
u = static_equilibrium(BeamAssembly, F, 'display', 'iter-detailed');
UNL = reshape(u,2,[]).';        % Nonlinear response

fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ... 
    '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 5;
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,UNL,'factor',scale,'color','k');
colormap jet
title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])

figure('units','normalized','position',[.2 .1 .6 .8])
scale = 5;
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,ULIN,'factor',scale,'color','k');
colormap jet
title(['LINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])



%% Dynamic response using Implicit Newmark
% forcing frequency of the average of first two natural frequencies
omega_ext = 2*2*pi*f0(2); 
T =  2*pi/omega_ext; % time period of forcing

% load amplification factor
amplification_factor = 1;

% forcing function
F_ext = @(t) amplification_factor * F * sin(omega_ext * t);

% Initial condition: equilibrium
u0 = zeros(BeamAssembly.Mesh.nDOFs, 1);
v0 = zeros(BeamAssembly.Mesh.nDOFs, 1);
a0 = zeros(BeamAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0)) 

q0 = BeamAssembly.constrain_vector(u0);
qd0 = BeamAssembly.constrain_vector(v0);
qdd0 = BeamAssembly.constrain_vector(a0);

% time step for integration
h = T/20;

% Precompute data for Assembly object
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.C = C; %rayleigh
%%
% Instantiate object for linear time integration
TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,BeamAssembly,F_ext);

% Linearized Time Integration
tmax = 1000*T; 
%tmax=0.002;
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin.Solution.u = BeamAssembly.unconstrain_vector(TI_lin.Solution.q);

% Animate solution on Mesh (very slow)
%AnimateFieldonDeformedMesh(myMesh.nodes,myMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:2,'filename','lineardisp')
%%
% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamAssembly,F_ext);

% Nonlinear Time Integration
 tmax = 1000*T; 
%tmax=0.002;

TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = BeamAssembly.unconstrain_vector(TI_NL.Solution.q);

%  save('TI_NL.mat','TI_NL');
%% Generalized alpha scheme
% linear
TI_lin_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7, 'linear',true);
TI_lin_alpha.Integrate(q0,qd0,qdd0,tmax,residual_lin);
TI_lin_alpha.Solution.u = BeamAssembly.unconstrain_vector(TI_lin_alpha.Solution.q);

% nonlinear
TI_NL_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
TI_NL_alpha.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL_alpha.Solution.u = BeamAssembly.unconstrain_vector(TI_NL_alpha.Solution.q);
%% Comparison
% Linear
%dof = 66; % random degree of freedom at which time response is compared
% 
% nfsense = find_node(9.25e-07,91.5e-6,[],nodes); % node where to see the results
% dof = get_index(nfsense, myMesh.nDOFPerNode );
% dof=dof(1);
% nfsense = find_node(148.85e-6,182.07e-6,[],nodes); % node where to see the results
% dof = get_index(nfsense, myMesh.nDOFPerNode );
% dof=dof(2);
% nfsense = find_node(16.775e-6,22.85e-6,[],nodes); % node where to see the results
% dof = get_index(nfsense, myMesh.nDOFPerNode );
% dof=dof(1);
nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);


figure;
plot(TI_lin.Solution.time, TI_lin.Solution.u(dof,:),'DisplayName', 'Full linear (Newmark)')
hold on
plot(TI_lin_alpha.Solution.time, TI_lin_alpha.Solution.u(dof,:),'DisplayName', 'Full linear (Generalized-$$\alpha$$)')
%plot(TI_lin_red.Solution.time, TI_lin_alpha.Solution.u(dof,:),'DisplayName', 'Reduced linear (Newmark)')
xlabel('q'); ylabel('time'); grid on; axis tight; legend('show')

% Nonlinear
% figure;
plot(TI_NL.Solution.time, TI_NL.Solution.u(dof,:),'DisplayName', 'Full nonlinear (Newmark)')
hold on
plot(TI_NL_alpha.Solution.time, TI_NL_alpha.Solution.u(dof,:),'DisplayName', 'Full nonlinear (Generalized-$$\alpha$$)')
%plot(TI_NL_alpha_red.Solution.time, TI_NL_alpha_red.Solution.u(dof,:),'DisplayName', 'Reduced nonlinear (Generalized-$$\alpha$$)')
xlabel('q'); ylabel('time'); grid on; axis tight; legend('show')

