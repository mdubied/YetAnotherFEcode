% ------------------------------------------------------------------------ 
% Script to test the implementation of actuation forces on 3D structures
% with TET4 elements, using a PROM formulation
% 
% Last modified: 11/03/2023, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

FORMULATION = 'N1'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 0;

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 0.6*263824;     % Young's modulus [Pa]
rho     = 1070;     % density [kg/m^3]
nu      = 0.499;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
myElementConstructor = @()Tet4Element(myMaterial);

% MESH_____________________________________________________________________
filename = 'testPartTET4D4';%'testPartTET4D4';
[nodes, elements, nset, elset] = mesh_ABAQUSread(filename);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);
%% PLOT MESH WITH NODES AND ELEMENTS
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 1);
hold off

%% BOUNDARY CONDITIONS
nel = size(elements,1);
nset = {};
Lx=10;
Ly=20;
Lz=20;

for el=1:nel   
    nodePosX = 0;
    for n=1:size(elements,2)
        if  nodes(elements(el,n),2)>Ly/2 && ~any(cat(2, nset{:}) == elements(el,n))
              nset{end+1} = elements(el,n);         
        end
    end   
end

for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)
end


%% ASSEMBLY ________________________________________________________________

% nominal
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;


%% DAMPING 
alfa = 3.1;
beta = 6.3*1e-6;
% nominal
Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);


%% EIGENMODES - VIBRATION MODES (VMs)

n_VMs = 2;

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
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMn(:,mod), 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0n(mod),3) ' Hz']);



%% MODAL DERIVATIVES (MDs)                     

% nominal
[MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);


%% ROM TENSORS (INTERNAL FORCES)                         
% define reduced order basis
Vn = [VMn MDn];     % reduced order basis (ROM-n)

% orthonormalize reduction basis
Vn = orth(Vn);	% ROM-n

% ROM-n: standard reduced order model (no defects in the mesh)
tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, Vn, USEJULIA);


%% REDUCED ASSEMBLIES

% ROM-n ___________________________________________________________________
ROMn_Assembly = ReducedAssembly(MeshNominal, Vn);
ROMn_Assembly.DATA.M = ROMn_Assembly.mass_matrix();  % reduced mass matrix (ROM-n)
ROMn_Assembly.DATA.C = Vn.'*Dn*Vn;    % reduced damping matrix (ROM-n), using C as needed by the residual function of Newmark integration
ROMn_Assembly.DATA.K = Vn.'*Kn*Vn;    % reduced stiffness matrix (ROM-n)




%% ROM TENSORS - ACTUATION FORCES
[skin,allfaces,skinElements, skinElementFaces] = getSkin3D(elements);

nel = size(elements,1);
actuationElements = zeros(nel,1);
for el=1:nel
    nodePosX = 0;
    if skinElements(el)==1
        for n=1:size(elements,2)
            if nodes(elements(el,n),3)>10.5 
                actuationElements(el) = 1;
%                 if nodePosX == 0
%                     nodePosX = 1;
%                 else
%                     actuationElements(el) = 1;
%                     disp(elements(el,n))
%                 end
                
            end
        end
    end
end


actuationDirection = [0;1;0;0;0;0]; % [0;1;0]-->[x^2,y^2,z^2,xy,xz,yz]

% ROM-n
tensors_actuation_ROMn = reduced_tensors_actuation_ROM(NominalAssembly, Vn, actuationElements, actuationDirection);

%%
actuationValue = 1.1;
disp=zeros(1041,3);
PlotFieldonDeformedMeshActuation(nodes,elements,actuationElements,0.8,disp,'factor',1,'color', 'k') ;


%% TIME INTEGRATION

% Parameters' initialization for all models _______________________________
% time step for integration
h = 0.05;

%% ROM-n __________________________________________________________________
% initial condition: equilibrium
q0 = zeros(size(Vn,2),1);
qd0 = zeros(size(Vn,2),1);
qdd0 = zeros(size(Vn,2),1);

% Actuation forces
B1 = tensors_actuation_ROMn.B1;
B2 = tensors_actuation_ROMn.B2;
F_ext = @(t,q) 1/2*(1-(1+1000000000*sin(t*2*pi)))*(B1+B2*q); % q is the reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_ROMn = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_actuation(q,qd,qdd,t,ROMn_Assembly,F_ext);

% nonlinear Time Integration
tmax = 1; 
TI_NL_ROMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMn.Solution.u = Vn * TI_NL_ROMn.Solution.q; % get full order solution

%% ROM-n __________________________________________________________________
actuationValues = zeros(size(TI_NL_ROMn.Solution.u,2),1);
for t=1:size(TI_NL_ROMn.Solution.u,2)
    actuationValues(t) = 1+0.003*sin(t*h*2*pi/10);
end

AnimateFieldonDeformedMeshActuation(nodes, elementPlot,actuationElements,actuationValues,TI_NL_ROMn.Solution.u,'factor',1,'index',1:3,'filename','result_video','framerate',1/h)

%% ROM-d __________________________________________________________________
% initial condition: equilibrium
q0 = zeros(size(Vd,2),1);
qd0 = zeros(size(Vd,2),1);
qdd0 = zeros(size(Vd,2),1);

% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_ROMd.Tr1) + ...
    double(tensors_hydro_ROMd.Tru2*q) + double(tensors_hydro_ROMd.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_ROMd.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMd.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMd.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration (TI)
TI_NL_ROMd = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMd_Assembly,F_ext);

% nonlinear time integration (TI)
TI_NL_ROMd.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMd.Solution.u = Vd * TI_NL_ROMd.Solution.q; % get full order solution

%% PROM - nominal (defect=0) ______________________________________________
% initial condition: equilibrium
q0 = zeros(size(V,2),1);
qd0 = zeros(size(V,2),1);
qdd0 = zeros(size(V,2),1);

% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_PROM.Tr1) + ...
    double(tensors_hydro_PROM.Tru2*q) + double(tensors_hydro_PROM.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_PROM = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,PROM_Assembly,F_ext);

% nonlinear Time Integration
tmax = 4.0; 
TI_NL_PROM.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

%% PROM-d (defect=xi) _____________________________________________________
% % initial condition: equilibrium
% q0 = zeros(size(V,2),1);
% qd0 = zeros(size(V,2),1);
% qdd0 = zeros(size(V,2),1);
% 
% % hydrodynamic forces
% F_ext = @(t,q,qd) (double(tensors_hydro_PROM.Tr1) + double(tensors_hydro_PROM.Tr2)*xi + ...
%     double(tensors_hydro_PROM.Tru2*q) + double(ttv(ttv(tensors_hydro_PROM.Tru3,xi,3), q,2)) + ...
%     double(tensors_hydro_PROM.Trudot2*qd) + double(ttv(ttv(tensors_hydro_PROM.Trudot3,xi,3), qd,2)) + ...
%     double(ttv(ttv(tensors_hydro_PROM.Truu3,q,3), q,2)) + ...
%     double(ttv(ttv(tensors_hydro_PROM.Truudot3,qd,3), q,2)) + ...
%     double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs
% 
% % instantiate object for nonlinear time integration
% TI_NL_PROMd = ImplicitNewmark('timestep',h,'alpha',0.005);
% 
% % modal nonlinear Residual evaluation function handle
% Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,PROM_Assembly,F_ext);
% 
% % nonlinear Time Integration
% tmax = 4.0; 
% TI_NL_PROMd.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
% TI_NL_PROMd.Solution.u = V * TI_NL_PROMd.Solution.q; % get full order solution

%% SENSITIVITY ANALYSIS ___________________________________________________
% initial conditions
s0 = zeros(size(V,2),1);
sd0 = zeros(size(V,2),1);
sdd0 = zeros(size(V,2),1);

% evaluate partial derivatives along solution
pd_fext_PROM = @(q,qd)DpROM_hydro_derivatives(q,qd,xi,tensors_hydro_PROM,FOURTHORDER);
pd_fint_PROM = @(q)DpROM_derivatives(q,tensors_PROM); 

% instantiate object for time integration
TI_sens = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true,'sens',true);

% residual function handle
Residual_sens = @(s,sd,sdd,t,it)residual_linear_sens(s,sd,sdd,t,PROM_Assembly,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, TI_NL_PROM.Solution.qdd,pd_fext_PROM,pd_fint_PROM,h);

% time integration (TI)
TI_sens.Integrate(q0,qd0,qdd0,tmax,Residual_sens);
TI_sens.Solution.s = V * TI_sens.Solution.q; % get full order solution


%% VISUALIZE ______________________________________________________________

%% 1-DOF PLOT _____________________________________________________________

% find a specific result node and corresponding DOF
rNodeDOFS = MeshNominal.get_DOF_from_location([0.9*Lx, 0]);
rNodeDOF = rNodeDOFS(2); % y-direction

% plot
figure 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
tplot=linspace(0,tmax,tmax/h);
plot(tplot,TI_NL_ROMn.Solution.u(rNodeDOF,1:end-2)*100)
hold on
plot(tplot,TI_NL_ROMd.Solution.u(rNodeDOF,1:end-2)*100, "-.")
plot(tplot,TI_NL_PROM.Solution.u(rNodeDOF,1:end-2)*100, "--")
%plot(tplot,TI_NL_PROMd.Solution.u(rNodeDOF,1:end-2)*100, ":")
%plot(tplot,(TI_NL_PROM.Solution.u(rNodeDOF,1:end-2)+TI_sens.Solution.s(rNodeDOF,1:end-2)*xi)*100,":")
approx = V*(TI_NL_PROM.Solution.q(:,1:end-2)+TI_sens.Solution.q(:,1:end-2)*xi)*100;
plot(tplot,approx(rNodeDOF,:))

title('Vertical displacement for a specific node')
ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
xlabel('Time [s]')
set(gca,'FontName','ComputerModern');
grid on
%legend({'ROM-n','ROM-d','PROM-n','PROM-d','Approximation of ROM-d using PROM sensitivities'}, 'Location', 'southoutside','Orientation','horizontal')
legend({'ROM-n','ROM-d','PROM-n','Approximation of ROM-d using PROM sensitivities'}, 'Location', 'southoutside','Orientation','horizontal')
hold off


%% Animations _____________________________________________________________
%% ROM-n __________________________________________________________________
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_ROMn.Solution.u,'factor',1,'index',1:2,'filename','result_video')

%% ROM-d __________________________________________________________________
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_ROMd.Solution.u,'factor',100,'index',1:2,'filename','result_video')


%% OPTIMIZATION ___________________________________________________________
d = [1;0];
%pd_fhydro_PROM = @(q,qd,xi)DpROM_hydro_derivatives(q,qd,xi,tensors_hydro_PROM);
[xiStar,LrEvo] = optimization_pipeline_1(V,d,tensors_PROM, tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd,TI_sens.Solution.q,TI_sens.Solution.qd)
%xi_star = optimization_pipeline_1(Vd,d,tensors_ROMd, tensors_hydro_ROMd,TI_NL_ROMd.Solution.q,TI_NL_ROMd.Solution.qd)

