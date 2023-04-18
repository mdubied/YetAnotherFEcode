% ------------------------------------------------------------------------ 
% Script to test and compare the optimization pipeline presented in the
% paper on a simple 2D (non-actuated) structure.
%
% Optimization objective: reduce drag on the airfoil, in the negative
%                         x-direction (foward swimming direction)
%
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

%% PREPARE (NOMINAL) MODEL AND SHAPE VARIATION                                                    

% DATA ____________________________________________________________________
E       = 260000;       %263824;       % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 
thickness = .1;         % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain

% Element
switch elementType
    case 'TRI3'
        myElementConstructor = @()Tri3Element(thickness, myMaterial);
end

% MESH_____________________________________________________________________

% nominal mesh
filename = 'naca0012TRI3_90Elements';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil

% plot nominal mesh
elementPlot = elements(:,1:3); 
figure('units','normalized','position',[.2 .1 .6 .4],'name','Nominal mesh')
PlotMesh(nodes, elementPlot, 0);

% boundary conditions: front end of the airfoil fixed through 2 nodes
frontNode = find_node_2D(0,0,nodes);
frontNode2 = find_node_2D(0.005,0,nodes);
nset = {frontNode, frontNode2};
MeshNominal.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)


% shape variations (only 1 in this example)
% (1) thinner airfoil 
nodes_projected = [nodes(:,1), nodes(:,2)*0];   % projection on x-axis
yDif = nodes_projected(:,2) - nodes(:,2);       % y-difference projection vs nominal
thinAirfoil = zeros(numel(nodes),1);            % create a single long vectors [x1 y1 x2 y2 ...]^T
thinAirfoil(2:2:end) = yDif;                    % fill up all y-positions

% shape variations basis
U = thinAirfoil;    % shape variations basis

%% OPTIMIZATION PARAMETERS ________________________________________________
d = [-1;0];
h = 0.05;
tmax = 1;
A=[1;-1];
b=[0.3;0.3];

%% OPTIMIZATION PIPELINE P1 _______________________________________________

tStart = tic;
[xiStar1,xiEvo1,LrEvo1] = optimization_pipeline_1(myElementConstructor, ...
    nset,nodes,elements,U,d,h,tmax,A,b,'maxIteration',200,'convCrit',0.002,'barrierParam',1e4,'gStepSize',0.1);  
tP1 = toc(tStart);
fprintf('Computation time for P1: %.2fs\n',tP1)

%% OPTIMIZATION PIPELINE P2 _______________________________________________

tStart = tic;
[xiStar2,xiEvo2,LrEvo2] = optimization_pipeline_2(MeshNominal, ...
    nodes,elements,U,d,h,tmax,A,b,'maxIteration',200,'convCrit',0.002,'barrierParam',1e4,'gStepSize',0.1);
tP2 = toc(tStart);
fprintf('Computation time for P2: %.2fs\n',tP2)

%% OPTIMIZATION PIPELINE P3 _______________________________________________

tStart = tic;
[xiStar3,xiEvo3,LrEvo3] = optimization_pipeline_3(MeshNominal, ...
    nodes,elements,U,d,h,tmax,A,b,'maxIteration',200,'convCrit',0.002,'barrierParam',1e4,'gStepSize',0.1);
tP3 = toc(tStart);
fprintf('Computation time for P3: %.2fs\n',tP3)

%% OPTIMIZATION PIPELINE P4 _______________________________________________

tStart = tic;
[xiStar4,xiEvo4,LrEvo4] = optimization_pipeline_4(myElementConstructor, ...
    nset,nodes,elements,U,d,h,tmax,A,b,'maxIteration',200,'convCrit',0.002,'barrierParam',1e4,'gStepSize',0.1,'nRebuild',8);
tP4 = toc(tStart);% nRebuild 20
fprintf('Computation time for P4: %.2fs\n',tP4)

%% OPTIMIZATION PIPELINE P5 _______________________________________________
% PFOM cannot be constructed. The memory will saturate wenn computing the
% PFOM internal force tensors.
% tStart = tic;
% [xiStar5,xiEvo5,LrEvo5] = optimization_pipeline_5(myElementConstructor, ...
%     nset,nodes,elements,U,d,h,tmax,A,b,'maxIteration',2);
% tP5 = toc(tStart);
% fprintf('Computation time for P5: %.2fs\n',tP5)

%% COST COMPUTATION ON FINAL ROM __________________________________________
xiFinal = xiStar4;%xiStar4;


% shape-varied mesh 
df = U*xiFinal;                       % displacement field introduced by shape variations
dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

% (P)ROM creation
FORMULATION = 'N1';VOLUME = 1; USEJULIA = 0;FOURTHORDER = 0; ACTUATION = 0;
[V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
        build_PROM(svMesh,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);
        
% solve EoMs
TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
eta = TI_NL_PROM.Solution.q;
etad = TI_NL_PROM.Solution.qd;
N = size(eta,2);

% compute cost Lr without barrier functions (no constraints, to obtain the
% the cost stemming from the hydrodynamic force only)
dr = reduced_constant_vector(d,V);
AFinal = [];     % no constraint
bFinal= [];      % no constraint
barrierParam = 3000;
Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xiFinal,dr,AFinal,bFinal,barrierParam);
fprintf('The cost function w/o constraint is: %.4f\n',Lr)

%% PLOT COST FUNCTION OVER ITERATIONS _____________________________________

figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(LrEvo2)
grid on
ylabel('$$L_r$$','Interpreter','latex')
xlabel('Iterations')

%% PLOT XI PARAMATER OVER ITERATIONS ______________________________________

figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(xiEvo2)
grid on
ylabel('$$\xi$$','Interpreter','latex')
xlabel('Iterations')



%% VISUALIZATION __________________________________________________________

% defected mesh
xi = xiFinal;       % parameter vector
m = length(xi);     % number of parameters

% update defected mesh nodes
d = U*xi;                       % displacement fields introduced by defects
dd = [d(1:2:end) d(2:2:end)]; 
nodes_defected = nodes + dd;    % nominal + d ---> defected 
DefectedMesh = Mesh(nodes_defected);
DefectedMesh.create_elements_table(elements,myElementConstructor);
DefectedMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

figure('units','normalized','position',[.2 .3 .6 .4],'name','Optimized mesh')
elementPlot = elements(:,1:3); hold on % plot only corners (otherwise it's a mess)
PlotMesh(nodes_defected, elementPlot, 0); 
PlotMesh(nodes,elementPlot,0);
v1 = reshape(U*xi, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow





