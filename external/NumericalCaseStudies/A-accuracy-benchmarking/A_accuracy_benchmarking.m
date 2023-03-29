% ------------------------------------------------------------------------ 
% Accuracy benchmarking.
% Script to benchmark the accuracy of the ROM and PROM formulations in
% comparison to the FOM solutions for 2D structures under hydrodynamic
% forces. Used element type: TRI3.
% 
% Last modified: 28/03/2023, Mathieu Dubied, ETH Zurich
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

% DATA ____________________________________________________________________
E       = 2600000;%263824;   % Young's modulus [Pa]
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
        filename = 'naca0012TRI_medium_mesh';%'naca0012TRI3_268Elements';
        [nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
end

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil

%% plot nominal mesh
elementPlot = elements(:,1:3); 
figure('units','normalized','position',[.2 .1 .6 .4],'name','Nominal mesh with element and node indexes')
PlotMesh(nodes, elementPlot, 0);
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

% plot the two meshes
figure('units','normalized','position',[.2 .3 .6 .4],'name','Shape-varied mesh');
elementPlot = elements(:,1:3); hold on 
PlotMesh(nodes_sv, elementPlot, 0); 
PlotMesh(nodes,elementPlot,0);
v1 = reshape(U*xi, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow

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

% shape-varied
svAssembly = Assembly(svMesh);
Msv = svAssembly.mass_matrix();
[Ksv,~] = svAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    svAssembly.DATA.K = Ksv;
    svAssembly.DATA.M = Msv;

%% DAMPING ________________________________________________________________
alfa = 0.912;
beta = 0.002;

% nominal
Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
NominalAssembly.DATA.C = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);

% shape-varied
Dsv = alfa*Msv + beta*Ksv; % Rayleigh damping 
svAssembly.DATA.D = Dsv;
svAssembly.DATA.C = Dsv;
Dsvc = svAssembly.constrain_matrix(Dsv);

%% EIGENMODES - VIBRATION MODES (VMs) _____________________________________

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

% plot
mod = 1;
elementPlot = elements(:,1:3); 
figure('units','normalized','position',[.2 .1 .6 .4],'name','Vibration mode for nominal mesh')
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMn(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0n(mod),3) ' Hz']);

% SHAPE-VARIED ____________________________________________________________
% eigentvalue problem
Ksvc = svAssembly.constrain_matrix(Ksv);
Msvc = svAssembly.constrain_matrix(Msv);
[VMsv,om] = eigs(Ksvc, Msvc, n_VMs, 'SM');
[f0d,ind] = sort(sqrt(diag(om))/2/pi);
VMsv = VMsv(:,ind);
for ii = 1:n_VMs
    VMsv(:,ii) = VMsv(:,ii)/max(sqrt(sum(VMsv(:,ii).^2,2)));
end
VMsv = svAssembly.unconstrain_vector(VMsv);

% plot
mod = 1;
elementPlot = elements(:,1:3); 
figure('units','normalized','position',[.2 .1 .6 .4],'name','Vibration mode for shape-varied mesh')
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMsv(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0d(mod),3) ' Hz'])

%% MODAL DERIVATIVES (MDs) ________________________________________________                   

% nominal
[MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);
% shape-varied
MDsv = modal_derivatives(svAssembly, elements, VMsv);
 
% shape variation sensitivities
[DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    FORMULATION);

%% ROM TENSORS (INTERNAL FORCES) __________________________________________                       
% define reduced order basis
Vn = [VMn MDn];     % reduced order basis (ROM-n)
Vsv = [VMsv MDsv];   	% reduced order basis (ROM-sv)
V  = [VMn MDn DS]; 	% reduced order basis (PROM) 

% orthonormalize reduction basis
Vn = orth(Vn);	% ROM-n
Vsv = orth(Vsv);	% ROM-sv
V  = orth(V);	% PROM

% ROM-n: reduced order model for nominal mesh
tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, Vn, USEJULIA);

% ROM-sv: reduced order model for shape-varied mesh
tensors_ROMsv = reduced_tensors_ROM(svAssembly, elements, Vsv, USEJULIA);
tensors_ROMsv.xi = xi; % save for which xi ROM-sv is computed

% PROM: parametric formulation for shape variations
tensors_PROM = reduced_tensors_DpROM(NominalAssembly, elements, ...
    V, U, FORMULATION, VOLUME, USEJULIA); %compute tensors

% evaluate the shape-varied tensors at xi
[Q2, Q3, Q4, Q3t, Q4t, M] = DefectedTensors(tensors_PROM, xi);

%% REDUCED ASSEMBLIES _____________________________________________________

% ROM-n ___________________________________________________________________
ROMn_Assembly = ReducedAssembly(MeshNominal, Vn);
ROMn_Assembly.DATA.M = ROMn_Assembly.mass_matrix();     % reduced mass matrix 
ROMn_Assembly.DATA.C = Vn.'*Dn*Vn;                      % reduced damping matrix 
ROMn_Assembly.DATA.K = Vn.'*Kn*Vn;                      % reduced stiffness matrix 

% ROM-sv __________________________________________________________________
ROMsv_Assembly = ReducedAssembly(svMesh, Vsv);
ROMsv_Assembly.DATA.M = ROMsv_Assembly.mass_matrix();   % reduced mass matrix 
ROMsv_Assembly.DATA.C = Vsv.'*Dsv*Vsv;                  % reduced damping matrix
ROMsv_Assembly.DATA.K = Vsv.'*Ksv*Vsv;                  % reduced stiffness matrix 

% PROM ____________________________________________________________________
PROM_Assembly = ReducedAssembly(MeshNominal, V);
PROM_Assembly.DATA.M = PROM_Assembly.mass_matrix();     % reduced mass matrix 
PROM_Assembly.DATA.C = V.'*Dn*V;                        % reduced damping matrix 
PROM_Assembly.DATA.K = V.'*Kn*V;                        % reduced stiffness matrix 

%% ROM TENSORS - HYDRODYNAMIC FORCES ______________________________________

[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);
vwater = [1;0.3];   % water velocity vector
rho = 997*0.03;

% ROM-n
tensors_hydro_ROMn = reduced_tensors_hydro_ROM(NominalAssembly, elements, Vn, skinElements, skinElementFaces, vwater, rho);
%%
% ROM-sv
tensors_hydro_ROMsv = reduced_tensors_hydro_ROM(svAssembly, elements, Vsv, skinElements, skinElementFaces, vwater, rho);
% PROM
FOURTHORDER = 1;
tensors_hydro_PROM = reduced_tensors_hydro_PROM(NominalAssembly, elements, V, U, FOURTHORDER, skinElements, skinElementFaces, vwater, rho);


%% FOM TENSORS - HYDRODYNAMIC FORCES (optional) ___________________________

FOM = 1;    % FOM=1 for assembling the hydrodynamic tensors at the assembly level (FOM)

if FOM == 1
    tensors_hydro_FOM = unreduced_tensors_hydro_FOM(NominalAssembly, elements, skinElements, skinElementFaces, vwater, rho);
    tensors_hydro_FOMsv = unreduced_tensors_hydro_FOM(svAssembly, elements, skinElements, skinElementFaces, vwater, rho);
end
%% TIME INTEGRATION _______________________________________________________

% Parameters' initialization for all models _______________________________
% time step for integration
h = 0.01;
tmax = 1.0; 

%% ROM-n __________________________________________________________________
% initial condition: equilibrium
tic
q0 = zeros(size(Vn,2),1);
qd0 = zeros(size(Vn,2),1);
qdd0 = zeros(size(Vn,2),1);

% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_ROMn.Tr1) + ...
    double(tensors_hydro_ROMn.Tru2*q) + double(tensors_hydro_ROMn.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_ROMn = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-8);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMn_Assembly,F_ext);

% nonlinear Time Integration
TI_NL_ROMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMn.Solution.u = Vn * TI_NL_ROMn.Solution.q; % get full order solution
toc

%% ROM-sv _________________________________________________________________
% initial condition: equilibrium
q0 = zeros(size(Vsv,2),1);
qd0 = zeros(size(Vsv,2),1);
qdd0 = zeros(size(Vsv,2),1);
tic
% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_ROMsv.Tr1) + ...
    double(tensors_hydro_ROMsv.Tru2*q) + double(tensors_hydro_ROMsv.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_ROMsv.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMsv.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMsv.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration (TI)
TI_NL_ROMsv = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-8);

% nonlinear residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMsv_Assembly,F_ext);

% nonlinear time integration (TI)
TI_NL_ROMsv.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMsv.Solution.u = Vsv * TI_NL_ROMsv.Solution.q; % get full order solution
toc

%% PROM - nominal (defect=0) ______________________________________________
% initial condition: equilibrium
q0 = zeros(size(V,2),1);
qd0 = zeros(size(V,2),1);
qdd0 = zeros(size(V,2),1);
tic
% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_PROM.Tr1) + ...
    double(tensors_hydro_PROM.Tru2*q) + double(tensors_hydro_PROM.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_PROM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-8);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,PROM_Assembly,F_ext);

% nonlinear Time Integration
TI_NL_PROM.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution
toc
%% SENSITIVITY ANALYSIS (PROM) ____________________________________________

% First order sensitivity _________________________________________________
% initial conditions
s0 = zeros(size(V,2),m);
sd0 = zeros(size(V,2),m);
sdd0 = zeros(size(V,2),m);

tic
% evaluate partial derivatives along solution
secondOrderDer = 1; % set to 1 if you want to compute the 2nd order derivative to perform a second order sensitivity analysis
pd_fext_PROM = @(q,qd)DpROM_hydro_derivatives(q,qd,xi,tensors_hydro_PROM,FOURTHORDER,secondOrderDer);
pd_fint_PROM = @(q)DpROM_derivatives(q,tensors_PROM); 

% instantiate object for time integration
TI_sens = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true,'sens',true,'MaxNRit',60,'RelTol',1e-6);

% residual function handle
Residual_sens = @(s,sd,sdd,t,it)residual_linear_sens(s,sd,sdd,t,PROM_Assembly,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, TI_NL_PROM.Solution.qdd,pd_fext_PROM,pd_fint_PROM,h);

% time integration (TI)
TI_sens.Integrate(s0,sd0,sdd0,tmax,Residual_sens);
TI_sens.Solution.s = V * TI_sens.Solution.q; % get full order solution
toc
tic
% Second order sensitivity ________________________________________________
if secondOrderDer == 1
    % initial conditions
    s20 = zeros(size(V,2),m,m);
    s2d0 = zeros(size(V,2),m,m);
    s2dd0 = zeros(size(V,2),m,m);
    
    
    % instantiate object for time integration
    TI_2ndSens = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true,'sens',true,'MaxNRit',60);
    
    % residual function handle
    Residual_2ndSens = @(s2,s2d,s2dd,t,it)residual_linear_2ndSens(s2,s2d,s2dd,t,PROM_Assembly,TI_NL_PROM.Solution, TI_sens.Solution,pd_fext_PROM,pd_fint_PROM,h);
    
    % time integration (TI)
    TI_2ndSens.Integrate(s20,s2d0,s2dd0,tmax,Residual_2ndSens);
    TI_2ndSens.Solution.s = V * TI_2ndSens.Solution.q; % get full order solution
end
toc
%% FOM-n (nominal) ________________________________________________________
% initial condition: equilibrium
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
F_ext = @(t,q,qd) (T1 + Tu2*q + Tudot2*q + double(ttv(ttv(Tuu3,q,3), q,2)) + ...
                    double(ttv(ttv(Tuudot3,qd,3),q,2))+ double(ttv(ttv(Tudotudot3,qd,3), qd,2)));

% instantiate object for nonlinear time integration
TI_NL_FOMn = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-8);

% modal nonlinear Residual evaluation function handle
Residual_NL = @(q,qd,qdd,t)residual_nonlinear_hydro(q,qd,qdd,t,NominalAssembly,F_ext);

% nonlinear Time Integration
TI_NL_FOMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
TI_NL_FOMn.Solution.u = zeros(NominalAssembly.Mesh.nDOFs,size(TI_NL_FOMn.Solution.q,2));
for t=1:size(TI_NL_FOMn.Solution.q,2)
    TI_NL_FOMn.Solution.u(:,t) = NominalAssembly.unconstrain_vector(TI_NL_FOMn.Solution.q(:,t));
end
toc
%% FOM-sv _________________________________________________________________
% initial condition: equilibrium
nUncDOFs = size(svMesh.EBC.unconstrainedDOFs,2);
q0 = zeros(nUncDOFs,1);
qd0 = zeros(nUncDOFs,1);
qdd0 = zeros(nUncDOFs,1);

tic
% hydrodynamic forces
T1 = svAssembly.constrain_vector(double(tensors_hydro_FOMsv.T1));
Tu2 = svAssembly.constrain_matrix(double(tensors_hydro_FOMsv.Tu2));
Tudot2 = svAssembly.constrain_matrix(double(tensors_hydro_FOMsv.Tudot2));
Tuu3 = tensor(svAssembly.constrain_tensor(double(tensors_hydro_FOMsv.Tuu3)));
Tuudot3 = tensor(svAssembly.constrain_tensor(double(tensors_hydro_FOMsv.Tuudot3)));
Tudotudot3 = tensor(svAssembly.constrain_tensor(double(tensors_hydro_FOMsv.Tudotudot3)));
F_ext = @(t,q,qd) (T1 + Tu2*q + Tudot2*q + double(ttv(ttv(Tuu3,q,3), q,2)) + ...
                    double(ttv(ttv(Tuudot3,qd,3),q,2))+ double(ttv(ttv(Tudotudot3,qd,3), qd,2)));

% instantiate object for nonlinear time integration
TI_NL_FOMsv = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-8');

% modal nonlinear Residual evaluation function handle
Residual_NL = @(q,qd,qdd,t)residual_nonlinear_hydro(q,qd,qdd,t,svAssembly,F_ext);

% nonlinear Time Integration
TI_NL_FOMsv.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
TI_NL_FOMsv.Solution.u = zeros(svAssembly.Mesh.nDOFs,size(TI_NL_FOMsv.Solution.q,2));
for t=1:size(TI_NL_FOMsv.Solution.q,2)
    TI_NL_FOMsv.Solution.u(:,t) = svAssembly.unconstrain_vector(TI_NL_FOMsv.Solution.q(:,t));
end
toc

%% FOM (original nonlinear force, small time step needed) _________________
% initial condition: equilibrium
nUncDOFs = size(MeshNominal.EBC.unconstrainedDOFs,2);
q0 = zeros(nUncDOFs,1);
qd0 = zeros(nUncDOFs,1);
qdd0 = zeros(nUncDOFs,1);

tic
% hydrodynamic forces
F_ext = @(t,q,qd) hydro_force_TRI3(NominalAssembly, skinElements, skinElementFaces, vwater, rho,q,qd);

% instantiate object for nonlinear time integration
TI_NL_FOMfull = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',200,'MaxNRit',60,'RelTol',1e-8);

% modal nonlinear Residual evaluation function handle
Residual_NL = @(q,qd,qdd,t)residual_nonlinear_hydro(q,qd,qdd,t,NominalAssembly,F_ext);

% nonlinear Time Integration
TI_NL_FOMfull.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
TI_NL_FOMfull.Solution.u = zeros(NominalAssembly.Mesh.nDOFs,size(TI_NL_FOMfull.Solution.q,2));
for t=1:size(TI_NL_FOMsv.Solution.q,2)
    TI_NL_FOMfull.Solution.u(:,t) = NominalAssembly.unconstrain_vector(TI_NL_FOMfull.Solution.q(:,t));
end
toc

%% 1-DOF PLOT _____________________________________________________________

% find a specific result node and corresponding DOF
tailNodeDOFS = MeshNominal.get_DOF_from_location([Lx, 0]);
tailNodeDOF = tailNodeDOFS(2); % y-direction

% plot
figure('units','normalized','position',[.1 .1 .8 .6],'name','Vertical displacement of the tail node')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
tplot=linspace(0,tmax,tmax/h+1);
%plot(tplot,TI_NL_FOMfull.Solution.u(tailNodeDOF,1:end-1)*100, "--")
hold on
%plot(tplot,TI_NL_FOMn.Solution.u(tailNodeDOF,1:end-1)*100, "--")
plot(tplot,TI_NL_ROMn.Solution.u(tailNodeDOF,1:end)*100)

%plot(tplot,TI_NL_FOMsv.Solution.u(tailNodeDOF,1:end-1)*100, "--")
plot(tplot,TI_NL_ROMsv.Solution.u(tailNodeDOF,1:end)*100)

qSol = TI_NL_PROM.Solution.q(:,1:end);
s1Sol = TI_sens.Solution.q(:,1:end);
s2Sol = TI_2ndSens.Solution.q(:,1:end);
approx = V*(qSol+s1Sol*xi)*100;
plot(tplot,approx(tailNodeDOF,:), "-.")
approx2 = approx;
for t=1:size(approx,2)
    approx2(:,t) = V*(qSol(:,t)+s1Sol(:,t)*xi+ 0.5*double(ttv(ttv(tensor(s2Sol(:,t),[size(s2Sol(:,t)) 1]),xi,3),xi,2)))*100;
end
plot(tplot,approx2(tailNodeDOF,:), "-.")

ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
xlabel('Time [s]')
set(gca,'FontName','ComputerModern');
grid on
%legend({'ROM-n','ROM-d','PROM-n','PROM-d','Approximation of ROM-d using PROM sensitivities'}, 'Location', 'southoutside','Orientation','horizontal')
legend({'FOM-n exact','FOM-n','ROM-n','FOM-sv','ROM-sv','Approx. 1st o. sens.','Approx. 2nd o. sens.'}, 'Location', 'southoutside','Orientation','horizontal')
%legend({'FOM-n','ROM-n','FOM-sv','ROM-sv','Approximation of ROM-sv using PROM','2nd order'}, 'Location', 'southoutside','Orientation','horizontal')
hold off


%% Animations _____________________________________________________________
%% ROM-n __________________________________________________________________

AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_ROMn.Solution.u,'factor',1,'index',1:2,'filename','ROMn','framerate',1/h)

%actuationValues = zeros(size(TI_NL_ROMn.Solution.u,2),1); % no
%AnimateFieldonDeformedMeshActuation(nodes, elementPlot,actuationElements,actuationValues,TI_NL_ROMn.Solution.u,'factor',1,'index',1:2,'filename','result_video','framerate',1/h)

%% ROM-d __________________________________________________________________
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_ROMsv.Solution.u,'factor',100,'index',1:2,'filename','result_video')

