% ------------------------------------------------------------------------ 
% 2D optimization of a fish.
% 
% Last modified: 02/04/2023, Mathieu Dubied, ETH Zurich
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
E       = 260000;   % Young's modulus [Pa]
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
xi1 = 0.1;  % thinning fish
xi2 = -0.2; % shortening fish
xi3 = 0.6;  % linear tail thinning
xi4 = 0.3;  % long ellipse tail thinning
xi5 = -0.1;  % short ellipse tail thinning
xi6 = 0;    % linear head thinning 
xi7 = 0.0;  % long ellipse head thinning
xi8 = 0.5;  % short ellipse head thinning
xi9 = 0.1;  % tail smoothening
xi10 = 0.3; % head smoothening

% xi1 = 0.2;  % thinning fish
% xi2 = 0.2; % shortening fish
% xi3 = -0.2;  % linear tail thinning
% xi4 = 0.4;  % long ellipse tail thinning
% xi5 = -0.5;  % short ellipse tail thinning
% xi6 = -0.2;    % linear head thinning 
% xi7 = 0.3;  % long ellipse head thinning
% xi8 = 0.1;  % short ellipse head thinning
% xi9 = -0.2;  % tail smoothening
% xi10 = 0.1; % head smoothening




% MESH ____________________________________________________________________

% nominal mesh
switch upper( whichModel )
    case 'ABAQUS'
        filename = 'rectangle192Elements';%'fishNominalTRI3';%'fishNominalTRI3_420El';
        [nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
end

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil

% plot nominal mesh
% elementPlot = elements(:,1:3); 
% figure('units','normalized','position',[.2 .1 .6 .4],'name','Nominal mesh with element and node indexes')
% PlotMesh(nodes, elementPlot, 0);

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

% shape variations 
[thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead,smoothTail,smoothHead] = shape_variations_2D(nodes,Lx,Ly);

% shape variations basis
U = [thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead,smoothTail,smoothHead];   % shape variations basis
xi = [xi1;xi2;xi3;xi4;xi5;xi6;xi7;xi8;xi9;xi10];            % shape variations parameters
%U = [thinFish,shortFish];   % shape variations basis
%xi = [xi1;xi2];            % shape variations parameters
m = length(xi);             % number of shape variations parameters


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

%% FORWARD SIMULATION ON NOMINAL MESH FOR VISUALIZATION PURPOSES __________

% ASSEMBLY ________________________________________________________________
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;

% DAMPING _________________________________________________________________
alfa = 0.912;
beta = 0.002;

Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
NominalAssembly.DATA.C = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);

% EIGENMODES - VIBRATION MODES (VMs) ______________________________________

n_VMs = 4;  % number of vibration modes to include in the ROM

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

% MODAL DERIVATIVES (MDs) _________________________________________________                   

[MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);

% ROM TENSORS (INTERNAL FORCES) ___________________________________________                       
% define reduced order basis
Vn = [VMn MDn];     % reduced order basis (ROM-n) 

% orthonormalize reduction basis
Vn = orth(Vn);	% ROM-n

% ROM-n: reduced order model for nominal mesh
tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, Vn, USEJULIA);

% REDUCED ASSEMBLIES ______________________________________________________

ROMn_Assembly = ReducedAssembly(MeshNominal, Vn);
ROMn_Assembly.DATA.M = ROMn_Assembly.mass_matrix();     % reduced mass matrix 
ROMn_Assembly.DATA.C = Vn.'*Dn*Vn;                      % reduced damping matrix 
ROMn_Assembly.DATA.K = Vn.'*Kn*Vn;                      % reduced stiffness matrix 

% ROM TENSORS - HYDRODYNAMIC FORCES _______________________________________

[~,~,skinElements, skinElementFaces] = getSkin2D(elements);
vwater = [0.7;0.21];%[1;0.3] %[0.7;0.21];   % water velocity vector
rho = 997*0.01;%997*0.001;
tensors_hydro_ROMn = reduced_tensors_hydro_ROM(NominalAssembly, elements, Vn, skinElements, skinElementFaces, vwater, rho);

% ROM TENSORS - ACTUATION FORCES __________________________________________

nel = size(elements,1);
actuationDirection = [1;0;0];%[1;0]-->[1;0;0] (Voigt notation)

% top muscle
topMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY>0.00 &&  elementCenterX > Lx*0.25
        topMuscle(el) = 1;
    end    
end
tensors_topMuscle_ROMn = reduced_tensors_actuation_ROM(NominalAssembly, Vn, topMuscle, actuationDirection);

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY<0.00 &&  elementCenterX > Lx*0.25
        bottomMuscle(el) = 1;
    end    
end
tensors_bottomMuscle_ROMn = reduced_tensors_actuation_ROM(NominalAssembly, Vn, bottomMuscle, actuationDirection);


%% TIME INTEGRATION ________________________________________________________

% parameters' initialization 
h = 0.1;
tmax = 20.0; 

% initial condition: equilibrium
q0 = zeros(size(Vn,2),1);
qd0 = zeros(size(Vn,2),1);
qdd0 = zeros(size(Vn,2),1);

k=1;
B1TopMuscle = tensors_topMuscle_ROMn.B1;
B2TopMuscle = tensors_topMuscle_ROMn.B2;
B1BottomMuscle = tensors_bottomMuscle_ROMn.B1;
B2BottomMuscle = tensors_bottomMuscle_ROMn.B2;
% hydrodynamic and actuation forces
F_ext = @(t,q,qd) (double(tensors_hydro_ROMn.Tr1) + ...
    double(tensors_hydro_ROMn.Tru2*q) + double(tensors_hydro_ROMn.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Trudotudot3,qd,3), qd,2)) + ...
    k/2*(1-(1+0.002*sin(t*2*pi/5)))*(B1TopMuscle+B2TopMuscle*q) + ...
    k/2*(1-(1-0.002*sin(t*2*pi/5)))*(B1BottomMuscle+B2BottomMuscle*q)); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_ROMn = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMn_Assembly,F_ext);

% nonlinear Time Integration
TI_NL_ROMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMn.Solution.u = Vn * TI_NL_ROMn.Solution.q; % get full order solution

% 1-DOF PLOT ______________________________________________________________

% find a specific result node and corresponding DOF
tailNodeDOFS = MeshNominal.get_DOF_from_location([Lx, 0]);
tailNodeDOF = tailNodeDOFS(2); % y-direction

% plot
figure('units','normalized','position',[.1 .1 .8 .6],'name','Vertical displacement of the tail node')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
tplot=linspace(0,tmax,size(TI_NL_ROMn.Solution.u,2));
plot(tplot,TI_NL_ROMn.Solution.u(tailNodeDOF,1:end)*100)
ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
xlabel('Time [s]')
set(gca,'FontName','ComputerModern');
grid on
legend({'ROM-n'}, 'Location', 'southoutside','Orientation','horizontal')
hold off


% Animation _______________________________________________________________

actuationValues = zeros(size(TI_NL_ROMn.Solution.u,2),1);
for t=1:size(TI_NL_ROMn.Solution.u,2)
    actuationValues(t) = 1+0.002*sin(t*h*2*pi/5);
end

actuationValues2 = zeros(size(TI_NL_ROMn.Solution.u,2),1);
for t=1:size(TI_NL_ROMn.Solution.u,2)
    actuationValues2(t) = 1-0.002*sin(t*h*2*pi/5);
end

AnimateFieldonDeformedMeshActuation2Muscles(nodes, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,TI_NL_ROMn.Solution.u, ...
    'factor',1,'index',1:2,'filename','result_video','framerate',1/h)





%% OPTIMIZATION PARAMETERS
dSwim = [-1;0]; %swimming direction
h = 0.05;
tmax = 3.0;
A=[1 0;-1 0;0 1;0 -1];
b=[0.3;0.3;0.3;0.3];

%% OPTIMIZATION PIPELINE P4 _______________________________________________
% U = [thinFish,shortFish,linearTail,longTail];
% A=[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1; 0 0 -1];
% b=[0.3;0.3;0.3;0.3;0.3;0.3];
A=[1 0 0 0;-1 0 0 0;0 1 0 0;0 -1 0 0;0 0 1 0; 0 0 -1 0; 0 0 0 1; 0 0 0 -1];
b=[0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3];
A=[ 1 0 0 0 0;
    -1 0 0 0 0;
    0 1 0 0 0;
    0 -1 0 0 0;
    0 0 1 0 0; 
    0 0 -1 0 0; 
    0 0 0 1 0; 
    0 0 0 -1 0;
    0 0 0 0 1;
    0 0 0 0 -1];
b=[0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3];
% A=[1;-1]
% b=[0.3,0.3]
U = [thinFish,shortFish,linearTail,shortHead,shortTail];
%%
tStart = tic;
[xiStar4,xiEvo4,LrEvo4] = optimization_pipeline_4(myElementConstructor,nset, ...
    nodes,elements,U,dSwim,h,tmax,A,b,'ACTUATION',1,'maxIteration',9,'convCrit',0.0001,'barrierParam',10,'gStepSize',0.001,'nRebuild',5);
tP4 = toc(tStart);
fprintf('Computation time: %.2fs\n',tP4)


%% VISUALIZATION __________________________________________________________

% shape-varied mesh 
df = U*xiStar4;                       % displacement field introduced by shape variations
dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

% plot the two meshes
figure('units','normalized','position',[.2 .3 .6 .4],'name','Shape-varied mesh');
elementPlot = elements(:,1:3); hold on 
PlotMesh(nodes_sv, elementPlot, 0); 
PlotMesh(nodes,elementPlot,0);
v1 = reshape(U*xiStar4, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow

%% COST COMPUTATION ON FINAL ROM __________________________________________
xiFinal = xiStar4;

% shape-varied mesh 
df = U*xiFinal;                       % displacement field introduced by shape variations
dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

% (P)ROM creation
FORMULATION = 'N1';VOLUME = 1; USEJULIA = 0;FOURTHORDER = 0; ACTUATION = 1;
[V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
        build_PROM(svMesh,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);

% %%
% % plot
% mod = 2;
% elementPlot = elements(:,1:3); 
% figure('units','normalized','position',[.2 .1 .6 .4],'name','Vibration mode for shape-varied mesh')
% PlotMesh(nodes, elementPlot, 0);
% v1 = reshape(V(:,mod), 2, []).';
% PlotFieldonDeformedMesh(nodes_sv, elementPlot, v1, 'factor', max(nodes_sv(:,2)));
% %title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0d(mod),3) ' Hz'])

%%        
% solve EoMs
TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax, ...
     'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
eta = TI_NL_PROM.Solution.q;
etad = TI_NL_PROM.Solution.qd;
N = size(eta,2);

%%
% compute cost Lr without barrier functions (no constraints, to obtain the
% the cost stemming from the hydrodynamic force only)
dr = reduced_constant_vector(dSwim,V);
AFinal = [];     % no constraint
bFinal= [];      % no constraint
barrierParam = 10;
Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xiFinal,dr,AFinal,bFinal,barrierParam);
%fprintf('The cost function w/o constraint is: %.4f\n',Lr)

%% PLOTs __________________________________________________________________

% find a specific result node and corresponding DOF
tailNodeDOFS = MeshNominal.get_DOF_from_location([Lx, 0]);
tailNodeDOF = tailNodeDOFS(2); % y-direction
% time axis
tplot=linspace(0,tmax,tmax/h+1);

% plot
figure('units','normalized','position',[.1 .1 .8 .6],'name','Vertical displacement of the tail node')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

plot(tplot,TI_NL_PROM.Solution.u(tailNodeDOF,1:end-1)*100, "--")

ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
xlabel('Time [s]')
set(gca,'FontName','ComputerModern');
grid on
%legend({'FOM','FOM-t','ROM','PROM'}, 'Location', 'eastoutside','Orientation','vertical')
hold off

%%
figure
plot(Lr)


%%
headnode = MeshNominal.get_DOF_from_location([0, 0]);
headnodeY = headnode(2);

% time axis
tplot=linspace(0,tmax,tmax/h+1);

% plot
figure('units','normalized','position',[.1 .1 .8 .6],'name','Vertical displacement of the tail node')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

plot(tplot,TI_NL_PROM.Solution.u(headnodeY,1:end-1)*100, "--")

ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
xlabel('Time [s]')
set(gca,'FontName','ComputerModern');
grid on
%legend({'FOM','FOM-t','ROM','PROM'}, 'Location', 'eastoutside','Orientation','vertical')
hold off

%% Animations _____________________________________________________________
% actuation elements
Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of the nominal fish
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of the nominal fish

nel = size(elements,1);
actuationDirection = [1;0;0];%[1;0]-->[1;0;0] (Voigt notation)

% top muscle
topMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY>0.00 &&  elementCenterX > Lx*0.25
        topMuscle(el) = 1;
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY<0.00 &&  elementCenterX > Lx*0.25
        bottomMuscle(el) = 1;
    end    
end

% actuation values
actuationValues1 = zeros(size(TI_NL_PROM.Solution.u,2),1);
actuationValues2 = zeros(size(TI_NL_PROM.Solution.u,2),1);
for t=1:size(TI_NL_PROM.Solution.u,2)
    actuationValues1(t) = 1+0.04*sin(t*h*2*pi);
    actuationValues2(t) = 1-0.04*sin(t*h*2*pi);
end
%ActuationElements,ActuationValues,ActuationElements2,ActuationValues2
AnimateFieldonDeformedMeshActuation2Muscles(nodes_sv, elementPlot,topMuscle,actuationValues1,bottomMuscle,actuationValues2,TI_NL_PROM.Solution.u,'factor',1,'index',1:2,'filename','optimised_shape','framerate',1/h)


%% COST COMPUTATION ON INITIAL ROM ________________________________________
% 88.57
% xiFinal = zeros(size(U,2),1);
% 
% % shape-varied mesh 
% df = U*xiFinal;                       % displacement field introduced by shape variations
% dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
% nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
% svMesh = Mesh(nodes_sv);
% svMesh.create_elements_table(elements,myElementConstructor);
% svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
% 
% % (P)ROM creation
% FORMULATION = 'N1';VOLUME = 1; USEJULIA = 0;FOURTHORDER = 0; ACTUATION = 1;
% [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
%         build_PROM(svMesh,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);
%         
% % solve EoMs
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
% eta = TI_NL_PROM.Solution.q;
% etad = TI_NL_PROM.Solution.qd;
% N = size(eta,2);
% 
% % compute cost Lr without barrier functions (no constraints, to obtain the
% % the cost stemming from the hydrodynamic force only)
% dr = reduced_constant_vector(dSwim,V);
% AFinal = [];     % no constraint
% bFinal= [];      % no constraint
% barrierParam = 10;
% Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xiFinal,dr,AFinal,bFinal,barrierParam);
% fprintf('The cost function w/o constraint is: %.4f\n',Lr)

