% EXAMPLE: beam meshed with 3D element
clear;
close all;
clc

whichModel = 'ANSYS'; % or "ABAQUS" or %CUSTOM


%% PREPARE MODEL

% DATA ____________________________________________________________________
E       = 148e9;     % Young's modulus [Pa]
rho     = 2329;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio
thickness = 1e-5;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________

switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'Job-BeamQuad';
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
    case 'ANSYS'
        % Alternatively, one can write an input file in ANSYS and read it as:
        
        meshfile='MeshFreqDividerFull_QUAD8_metric.dat';
        meshtype='QUAD8';
        model=Ansys2Matlab_mesh(meshfile,meshtype);
        nodes=model.nodes(:,2:3);          % 2:3 because 2D
        elements=model.elements(:,2:end);    
end

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
switch upper( whichModel )
    case 'CUSTOM'
        myMesh.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
    case 'ABAQUS'
        myMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
    case 'ANSYS'
        fixed_nodes = model.fixedNodes(:)';
        myMesh.set_essential_boundary_condition((fixed_nodes),1:2,0);
end

% plot
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
hold on
plot(nodes(fixed_nodes,1),nodes(fixed_nodes,2),'r*')
%plot3(nodes(fixed_nodes,1),nodes(fixed_nodes,2),nodes(fixed_nodes,3),'r*')
hold off
%%
% ASSEMBLY ________________________________________________________________
FreqDivAssembly = Assembly(myMesh);
M = FreqDivAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = FreqDivAssembly.tangent_stiffness_and_force(u0);

 C = 761*M + 1.96e-10*K;


% store matrices
FreqDivAssembly.DATA.K = K;
FreqDivAssembly.DATA.M = M


%% EXAMPLE 1

% Eigenvalue problem_______________________________________________________
n_VMs = 8; % first n_VMs modes with lowest frequency calculated
Kc = FreqDivAssembly.constrain_matrix(K);
Mc = FreqDivAssembly.constrain_matrix(M);
[V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = FreqDivAssembly.unconstrain_vector(V0);

% PLOT
for mod=1:5
% mod = 4;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(V0(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', 0.00001);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])

end
%% EXAMPLE 2

% Define external force:
% % Body force
% Pressure = 1e6;
% F = Pressure*FreqDivAssembly.uniform_body_force();

% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(9.25e-07,0.000183,[],nodes); % node where to put the force
%nf=33;
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(2)) = 3e-4;

u_lin = FreqDivAssembly.solve_system(K, F);
ULIN = reshape(u_lin,2,[]).';	% Linear response
u = static_equilibrium(FreqDivAssembly, F, 'display', 'iter-detailed');
UNL = reshape(u,2,[]).';        % Nonlinear response

fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ...
    '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 7.9e3;
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,UNL,'factor',scale,'color','k');
colormap jet
title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])

%% Dynamic response using Implicit Newmark
% forcing frequency of the average of first two natural frequencies
omega_ext = 2*pi*2*f0(4); 
T =  2*pi/omega_ext; % time period of forcing

% load amplification factor
amplification_factor = 1;

% forcing function
F_ext = @(t) amplification_factor * F * sin(omega_ext * t);

% Initial condition: equilibrium
u0 = zeros(FreqDivAssembly.Mesh.nDOFs, 1);
v0 = zeros(FreqDivAssembly.Mesh.nDOFs, 1);
a0 = zeros(FreqDivAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0)) 

q0 = FreqDivAssembly.constrain_vector(u0);
qd0 = FreqDivAssembly.constrain_vector(v0);
qdd0 = FreqDivAssembly.constrain_vector(a0);

% time step for integration
h = T/20;

% Precompute data for Assembly object
FreqDivAssembly.DATA.M = M;
FreqDivAssembly.DATA.K = K;
FreqDivAssembly.DATA.C = C; %rayleigh
%%
% Instantiate object for linear time integration
TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,FreqDivAssembly,F_ext);

% Linearized Time Integration
tmax = 10*T; 
%tmax=0.002;
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin.Solution.u = FreqDivAssembly.unconstrain_vector(TI_lin.Solution.q);

% Animate solution on Mesh (very slow)
%AnimateFieldonDeformedMesh(myMesh.nodes,myMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:2,'filename','lineardisp')
%%
% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,FreqDivAssembly,F_ext);

% Nonlinear Time Integration
%  tmax = 10*T; 
tmax=0.0016;

TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = FreqDivAssembly.unconstrain_vector(TI_NL.Solution.q);

NL_SOLUTION=TI_NL.Solution;
save('TI_NL','-struct','NL_SOLUTION');
%% Generalized alpha scheme
% linear
TI_lin_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7, 'linear',true);
TI_lin_alpha.Integrate(q0,qd0,qdd0,tmax,residual_lin);
TI_lin_alpha.Solution.u = FreqDivAssembly.unconstrain_vector(TI_lin_alpha.Solution.q);

% nonlinear
TI_NL_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
TI_NL_alpha.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL_alpha.Solution.u = FreqDivAssembly.unconstrain_vector(TI_NL_alpha.Solution.q);
%% Comparison
% Linear
%dof = 66; % random degree of freedom at which time response is compared
% 
nfsense = find_node(9.25e-07,91.5e-6,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(1);
% nfsense = find_node(148.85e-6,182.07e-6,[],nodes); % node where to see the results
% dof = get_index(nfsense, myMesh.nDOFPerNode );
% dof=dof(2);
% nfsense = find_node(16.775e-6,22.85e-6,[],nodes); % node where to see the results
% dof = get_index(nfsense, myMesh.nDOFPerNode );
% dof=dof(1);


figure;
plot(TI_lin.Solution.time, TI_lin.Solution.u(dof,:),'DisplayName', 'Full linear (Newmark)')
hold on
plot(TI_lin_alpha.Solution.time, TI_lin_alpha.Solution.u(dof,:),'DisplayName', 'Full linear (Generalized-$$\alpha$$)')
%plot(TI_lin_red.Solution.time, TI_lin_alpha.Solution.u(dof,:),'DisplayName', 'Reduced linear (Newmark)')
xlabel('q'); ylabel('time'); grid on; axis tight; legend('show')

% Nonlinear
figure;
plot(TI_NL.Solution.time, TI_NL.Solution.u(dof,:),'DisplayName', 'Full nonlinear (Newmark)')
hold on
plot(TI_NL_alpha.Solution.time, TI_NL_alpha.Solution.u(dof,:),'DisplayName', 'Full nonlinear (Generalized-$$\alpha$$)')
%plot(TI_NL_alpha_red.Solution.time, TI_NL_alpha_red.Solution.u(dof,:),'DisplayName', 'Reduced nonlinear (Generalized-$$\alpha$$)')
xlabel('q'); ylabel('time'); grid on; axis tight; legend('show')

