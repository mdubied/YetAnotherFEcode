% ------------------------------------------------------------------------ 
% 3D shape optimisation of a fish.
% 
% Last modified: 19/11/2023, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 1;

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 2600000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain

% element
switch elementType
    case 'TET4'
        myElementConstructor = @()Tet4Element(myMaterial);
    case 'TET10'
        myElementConstructor = @()Tet10Element(myMaterial);
end
% MESH ____________________________________________________________________

% nominal mesh
filename = '3d_rectangle_660el';%'fish3_664el';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);

nodes = nodes*0.01;
nodes(:,2)=0.8*nodes(:,2);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil

% plot nominal mesh
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off

% boundary conditions of nominal mesh
nel = size(elements,1);
nset = {};
for el=1:nel   
    for n=1:size(elements,2)
        if  nodes(elements(el,n),1)>-Lx*0.1 && ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
    end   
end

%% SHAPE VARIATIONS _______________________________________________________

[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% shape variations basis
% U = [thinFish,fishTailsv,fishHeadsv];    % shape variations basis
% U = [y_thinFish,z_tail, z_head, y_head,y_ellipseFish];    % shape variations basis
% U = fishEllipseYZ;

% SO1
U = [z_tail,z_head,y_thinFish];

% SO2
U = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];

% SO3

U = [z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_head,y_linLongTail,y_ellipseFish];


% plot the two meshes
% xiPlot = [0.23;-0.39;0.1091];
% xiPlot = [-0.6;0.3;0.5;0.3;0.5];
xiPlot = [0.2;-0.6;0.2;0.1;0.3;0.2;0.2;0.4];
% xiPlot = 0.5;
f1 = figure('units','centimeters','position',[3 3 10 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:4); hold on 
v1 = reshape(U*xiPlot, 3, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);
axis equal; grid on; box on; 
set(f1,'PaperUnits','centimeters');
% set(f1,'PapeyPositionMode','auto');
% set(f1,'PaperSize',[7 3.5]); % Canvas Size
set(f1,'Units','centimeters');


%%
xiTest = [-0.4;0.2;0.4];
% shape-varied mesh 
df = U*xiTest;                       % displacement field introduced by shape variations
dd = [df(1:3:end) df(2:3:end) df(3:3:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);


%%
nNodes = size(nodes,1);
matrix_inp = [linspace(1,nNodes,nNodes)',nodes_sv];
matrix_inp(:,1) = cast(matrix_inp(:,1),"uint8");
disp(matrix_inp)
writematrix(matrix_inp,'M.csv')
 
%% OPTIMIZATION PARAMETERS
h = 0.005;
tmax = 2.0;

%% OPTIMISATION SO1 _______________________________________________________
A = [1 0 0 ;
    -1 0 0;
    0 1 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
b = [0.5;0.5;0.5;0.5;0.5;0.5];

tStart = tic;
[xiStar,xiEvo,LEvo, LwoBEvo] = optimise_shape_3D(myElementConstructor,nset, ...
    nodes,elements,U,h,tmax,A,b,'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, 'maxIteration',25,'convCrit',0.004,'convCritCost',0.8,'barrierParam',10, ...
    'gStepSize',0.002,'nRebuild',6, 'rebuildThreshold',0.15,'USEJULIA',0);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)


%% OPTIMISATION SO2 _______________________________________________________

% [z_tail,z_head,y_thinFish,y_linLongTail,y_head,y_ellipseFish];
% Constraints
nParam = 5;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
yTotConstr = [0 0 1 0 1;0 0  -1 0 -1;
              0 0 0 1 1;0 0  0 -1 -1];
A = [A;yTotConstr];
disp(A);
b = [0.6;0.6;
    0.6;0.6;
    0.5;0.5;
    0.5;0.5;
    0.5;0.5;
    0.8;0.8;
    0.8;0.8];

%%

tStart = tic;
[xiStar,xiEvo,LEvo, LwoBEvo] = optimise_shape_3D(myElementConstructor,nset, ...
    nodes,elements,U,h,tmax,A,b,'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, 'maxIteration',25,'convCrit',0.004,'convCritCost',0.8,'barrierParam',1, ...
    'gStepSize',0.001,'nRebuild',6, 'rebuildThreshold',0.15,'USEJULIA',1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)

%% OPTIMISATION SO5 _______________________________________________________

% [y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
%    y_tail,y_head,y_linLongTail,y_ellipseFish]
% Constraints
nParam = 10;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
yTotConstr = [1 0 0 0 0 0 1 0 1 1;-1 0 0 0 0 0 -1 0 -1 -1;
              1 0 0 0 0 0 0 1 0 1;-1 0 0 0 0 0 0 -1 0 -1];
A = [A;yTotConstr];
disp(A);
b = [0.2;0.2;
    0.3;0.3;
    0.6;0.4;
    0.5;0.5;
    0.3;0.3;
    0.5;0.5;
    0.4;0.4;
    0.4;0.4;
    0.4;0.4;
    0.4;0.4;
    0.8;0.8;
    0.8;0.8];

%% OPTIMISATION SO5 _______________________________________________________

% [,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
%    y_head,y_linLongTail,y_ellipseFish]
% Constraints
nParam = 8;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
yTotConstr = [0 0 0 0 0 0 1 1;0 0 0 0 0 0 -1 -1;
              0 0 0 0 0 1 0 1;0 0 0 0 0 -1 0 -1];
A = [A;yTotConstr];
disp(A);
b = [0.3;0.3;
    0.5;0.5;
    0.5;0.5;
    0.3;0.3;
    0.5;0.5;
    0.4;0.4;
    0.5;0.5;
    0.4;0.4;
    0.9;0.9;
    0.9;0.9];

%%

tStart = tic;
[xiStar,xiEvo,LEvo, LwoBEvo] = optimise_shape_3D(myElementConstructor,nset, ...
    nodes,elements,U,h,tmax,A,b,'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, 'maxIteration',40,'convCrit',0.004,'convCritCost',1,'barrierParam',3, ...
    'gStepSize',0.0006,'nRebuild',12, 'rebuildThreshold',0.15,'USEJULIA',1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)
%% PLOT SHAPE VARIATIONS AND OPTIMAL SHAPE ________________________________
f1 = figure('units','centimeters','position',[3 3 9 7]);
elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
pos1 = [0,0.5,0.5,0.5];
pos2 = [0.5,0.5,0.5,0.5];
pos3 = [0,0,0.5,0.5];
pos4 = [0.5,0,0.5,0.5];
textPosX = -0.2;
textPosY = -0.05;
textPosZ = -0.14;

% shape variation 1
ax1 = subplot(2,2,1,'Position',pos1);
subU = U(:,1);
xiPlot = 0.5;

v1 = reshape(subU*xiPlot, 3, []).';
dm = PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
subplotName = strcat('$$\xi_7=',num2str(xiPlot),'$$');
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')

% shape variation 2
ax2 = subplot(2,2,2,'Position',pos2);
subU = U(:,2);
xiPlot = 0.5;

v2 = reshape(subU*xiPlot, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v2, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
subplotName = strcat('$$\xi_8=',num2str(xiPlot),'$$');
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')

% shape variation 3
ax3 = subplot(2,2,3,'Position',pos3);
subU = U(:,3);
xiPlot = 0.5;

v3 = reshape(subU*xiPlot, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v3, 'factor', 1); 
plotcube(L,O,.05,[0 0 0]);
subplotName = strcat('$$\xi_9=',num2str(xiPlot),'$$');
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')

% optimal shape
ax4 = subplot(2,2,4,'Position',pos4);
xiPlot = xiStar;
xiPlotName = strcat('[',num2str(xiStar(1)),', ',num2str(xiStar(2)),', ',num2str(xiStar(3)),']^\top$$');

v1 = reshape(U*xiPlot, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
subplotName = '$$\mathbf{\xi}^\ast$$';
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')


axis([ax1 ax2 ax3 ax4],[-0.35 0 -0.04 0.04 -0.16 0.16])
% set(ax2, 'box', 'on', 'Visible', 'on')
% set(ax1, 'box', 'on', 'Visible', 'on')

exportgraphics(f1,'SO5_shapes_V0.pdf','Resolution',600)

%% PLOT COST FUNCTION WITH PARAMETERS _____________________________________
f2 = figure('units','centimeters','position',[3 3 9 5]);
t = tiledlayout(1,2);
t.TileSpacing = 'loose';
t.Padding = 'tight';

% cost function
ax1 = nexttile;
plot(LEvo)
grid on
ylabel('$$L$$','Interpreter','latex')
xlabel('Iterations')

% Parameters
ax2 = nexttile;
plot(xiEvo(1,:));%,LineStyle,"-");
hold on
plot(xiEvo(2,:));%,LineStyle,"--");
plot(xiEvo(3,:));%,LineStyle,"-.");
%plot(xiEvo(4,:));%,LineStyle,"-");
%plot(xiEvo(5,:));%,LineStyle,"--");
% plot(xiEvo(6,:),LineStyle="-.");

grid on
ylabel('$$\xi$$','Interpreter','latex')
xlabel('Iterations')
legend('$$\xi_1$$','$$\xi_2$$','$$\xi_3$$','Interpreter','latex', ...
    'Position',[0.75 0.35 0.2 0.2])

legend('$$\xi_1$$','$$\xi_2$$','$$\xi_4$$','$$\xi_5$$','$$\xi_6$$','Interpreter','latex', ...
    'Position',[0.35 0.64 0.2 0.35])
exportgraphics(f2,'SO2_evo_V0.pdf','Resolution',600)


%% PLOT PARAMETERS' EVOLUTION OVER ITERATIONS _____________________________
figure('Position',[100,100,600,200])
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
subplot(1,3,1)
plot(xiEvo(1,:));
ylabel('$$\xi_1$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,2)
plot(xiEvo(2,:));
ylabel('$$\xi_2$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,3)
plot(xiEvo(3,:));
ylabel('$$\xi_3$$','Interpreter','latex')
xlabel('Iterations')
grid on


%% COST COMPUTATION ON FINAL ROM __________________________________________
%xiStar4 = [0.2;0.1;0.5;0.2;0.1]
xiFinal = xiStar;

% tmax=4;
% h=0.0025

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
       
% solve EoMs
TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax, ...
     'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
eta = TI_NL_PROM.Solution.q;
etad = TI_NL_PROM.Solution.qd;
N = size(eta,2);

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
tailNodeDOF = tailNodeDOFS(1); % y-direction
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



%% ANIMATIONS _____________________________________________________________
