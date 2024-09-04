% ------------------------------------------------------------------------ 
% 3D shape optimisation of a fish.
% 
% Last modified: 03/09/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc
set(groot,'defaulttextinterpreter','latex');

%% PREPARE MODELS _________________________________________________________                                                   

% load material parameters
load('parameters.mat') 

% specify and create FE mesh
filename ='3d_rectangle_8086el' ;%'3d_rectangle_8086el'; %'3d_rectangle_8086el'
%'3d_rectangle_1272el';%'3d_rectangle_1272el';%'3d_rectangle_660el'; 
kActu = 1.0;    % multiplicative factor for the actuation forces, dependent on the mesh
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions(nodes);


% plot nominal mesh
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off


%% SHAPE VARIATIONS _______________________________________________________

[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, x_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% SO1
U = [z_tail,z_head,y_thinFish];

% SO2
% U = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
% 
% % SO3

% U = [z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
%     y_head,y_linLongTail,y_ellipseFish];

% test 
% U = [x_concaveTail,z_tail];


% plot the two meshes
% xiPlot = [0.23;-0.39;0.1091];
% xiPlot = [-0.6;0.3;0.5;0.3;0.5];
% xiPlot = [0.2;-0.6;0.2;0.1;0.3;0.2;0.2;0.4];
xiPlot = [-0.5;0.5,;0.5];
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


%% OPTIMIZATION PARAMETERS
h = 0.01;
tmax = 2.0;
kActu = 0.1;    % tested with 3.0 for 660el, 1.3 for 1272el, 0.25 for 4270el, 0.1 for 8086

%% OPTIMISATION SO1 _______________________________________________________

U = [z_tail,z_head,y_thinFish]; 
% U = [z_tail,z_head, x_concaveTail];
A = [1 0 0 ;
    -1 0 0;
    0 1 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
b = [0.5;0.5;0.5;0.5;0.3;0.3];

barrierParam = 10*ones(1,length(b));

tStart = tic;
[xiStar,xiEvo,LEvo, LwoBEvo, nIt] = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,kActu,U,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',35, ...
    'convCrit',0.004, ...
    'convCritCost',1.0, ... % 1.0 for 8086, 0.3 below
    'barrierParam',barrierParam, ...
    'gStepSize',0.0005, ...  % 0.0005 for 8086, 0.002 below
    'nRebuild',6, ...
    'rebuildThreshold',0.15,...
    'USEJULIA',1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)
fprintf('Number of built models and solved EoMs: %5d\n',nIt)
fprintf('Computation time per models/EoMs: %.2f\n',topti/60/nIt)
filename='SO1_results';
save(filename,'xiStar','xiEvo','LEvo','LwoBEvo','topti')

%% OPTIMISATION SO2 _______________________________________________________

U = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
% Constraints
nParam = 5;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
yTotConstr = [0 0 1 0 1;0 0  -1 0 -1;
              0 0 0 1 1;0 0  0 -1 -1];
A = [A;yTotConstr];
b = [0.5;0.5;
    0.5;0.5;
    0.4;0.4;
    0.5;0.5;
    0.5;0.5;
    0.8;0.8;
    0.8;0.8];

barrierParam = ones(1,length(b));

tStart = tic;
[xiStar,xiEvo,LEvo, LwoBEvoSO2,nIt] = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,kActu,U,h,tmax,A,b,...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ...
    'convCrit',0.004, ...
    'convCritCost',0.5, ... % set to 0.5 for 4270 el, 0.1 below
    'barrierParam',barrierParam, ...
    'gStepSize',0.0005,...   % set to 0.001 for 4270 el, 0.03 below, 0.005 for 8086
    'nRebuild',5, ...
    'rebuildThreshold',0.15,...
    'USEJULIA',1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)
fprintf('Number of built models and solved EoMs: %5d\n',nIt)
fprintf('Computation time per models/EoMs: %.2f\n',topti/60/nIt)
filename='SO2_results';
save(filename,'xiStar','xiEvo','LEvo','LwoBEvoSO2','topti')

%% OPTIMISATION SO5 _______________________________________________________

%[0.1968;-0.4929;0.3093;-0.2883;0.2715;0.4913;0.4962;0.3890]

U = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
    z_smallFish, z_notch, x_concaveTail];

% Constraints
nParam = 8;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
% yTotConstr = [0 0 0 0 0 0 1 1;0 0 0 0 0 0 -1 -1;
%               0 0 0 0 0 1 0 1;0 0 0 0 0 -1 0 -1];
% A = [A;yTotConstr];

b = [0.4;0.4;
    0.4;0.4;
    0.4;0.4;
    0.3;0.3;
    0.5;0.5;
    0.3;0.3;
    0.3;0.3;
    0.3;0.01];  % concave tail only in one direction 
barrierParam = 3*ones(1,length(b));

%%

tStart = tic;
[xiStar,xiEvo,LEvo, LwoBEvoSO5,nIt] = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,kActu,U,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',40, ...
    'convCrit',0.015, ...    % set to 0.015 for 8086el, 0.01 below
    'convCritCost',5.0, ...% set to 0.5 for 4270 el, 0.2 below, 5.0 for 8086
    'barrierParam',barrierParam, ...
    'gStepSize',0.0002, ... % set to 0.005 for 4270 el, 0.003 below, 0.0002 for 8086
    'nRebuild',5, ...
    'rebuildThreshold',0.15, ...
    'USEJULIA',1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)
fprintf('Number of built models and solved EoMs: %5d\n',nIt)
fprintf('Computation time per models/EoMs: %.2f\n',topti/60/nIt)

filename='SO3_results';
save(filename,'xiStar','xiEvo','LEvo','LwoBEvoSO5','topti')

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
subU = U(:,6);
xiPlot = 0.5;

v1 = reshape(subU*xiPlot, 3, []).';
dm = PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
subplotName = strcat('$$\xi_7=',num2str(xiPlot),'$$');
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')

% shape variation 2
ax2 = subplot(2,2,2,'Position',pos2);
subU = U(:,7);
xiPlot = 0.5;

v2 = reshape(subU*xiPlot, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v2, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
subplotName = strcat('$$\xi_8=',num2str(xiPlot),'$$');
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')

% shape variation 3
ax3 = subplot(2,2,3,'Position',pos3);
subU = U(:,8);
xiPlot = 0.5;

v3 = reshape(subU*xiPlot, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v3, 'factor', 1); 
plotcube(L,O,.05,[0 0 0]);
subplotName = strcat('$$\xi_9=',num2str(xiPlot),'$$');
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')

% optimal shape
ax4 = subplot(2,2,4,'Position',pos4);
xiPlot = xiStar;% [-0.3;0.3;0.1;0.4;0.4;0.1;0.2;0.1];%xiStar;
% [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
%     z_smallFish, z_notch, x_concaveTail];
% xiPlot(1) = 0.2;
xiPlotName = strcat('[',num2str(xiStar(1)),', ',num2str(xiStar(2)),', ',num2str(xiStar(3)),']^\top$$');

v1 = reshape(U*xiPlot, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
subplotName = '$$\mathbf{\xi}^\ast$$';
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')


axis([ax1 ax2 ax3 ax4],[-0.40 0 -0.04 0.04 -0.16 0.16])
% set(ax2, 'box', 'on', 'Visible', 'on')
% set(ax1, 'box', 'on', 'Visible', 'on')

% exportgraphics(f1,'SO3_shapes_V1.pdf','Resolution',1200)
exportgraphics(f1,'SO3_shapes_V1.jpg','Resolution',600)

%% PLOT COST FUNCTION WITH PARAMETERS _____________________________________
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
f2 = figure('units','centimeters','position',[3 3 9 5]);
t = tiledlayout(1,2);
t.TileSpacing = 'loose';
t.Padding = 'tight';

% cost function
ax1 = nexttile;
plot(LwoBEvo)
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
    'Position',[0.8 0.44 0.2 0.2])

% legend('$$\xi_1$$','$$\xi_2$$','$$\xi_4$$','$$\xi_5$$','$$\xi_6$$','Interpreter','latex', ...
%     'Position',[0.35 0.64 0.2 0.35])
% exportgraphics(f2,'SO3_evo_V1.pdf','Resolution',1200)
exportgraphics(f2,'SO3_evo_V1.jpg','Resolution',600)

%% PLOT COST FUNCTIONS, SO2 + SO3 _________________________________________
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
f2 = figure('units','centimeters','position',[3 3 9 5]);
t = tiledlayout(1,2);
t.TileSpacing = 'loose';
t.Padding = 'tight';

% cost function
ax1 = nexttile;
plot(LwoBEvoSO2)
hold on
annotation('textbox',[.35 .6 .3 .3],'String','SO2','EdgeColor','None','Interpreter','latex');
grid on
ylabel('$$L$$','Interpreter','latex')
xlabel('Iterations','Interpreter','latex')
ylim([-1200,-18])
hold off


% Parameters
ax2 = nexttile;
plot(LwoBEvoSO5)
hold on
annotation('textbox',[.87 .6 .3 .3],'String','SO3','EdgeColor','None','Interpreter','latex');
grid on
ylabel('$$L$$','Interpreter','latex')
xlabel('Iterations','Interpreter','latex')
ylim([-1200,-18])
hold off

% exportgraphics(f2,'SO2-3_cost_V1.pdf','Resolution',600)
exportgraphics(f2,'SO2-3_cost_V1.jpg','Resolution',600)

%%
xiTest = [0;0;0];%[-0.4;0.2;0.4];
xiTest=[0.5;-0.5]
% shape-varied mesh 
df = U*xiTest;                       % displacement field introduced by shape variations
dd = [df(1:3:end) df(2:3:end) df(3:3:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);

elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes_sv, elementPlot, 0);
hold off


%%
nNodes = size(nodes,1);
matrix_inp = [linspace(1,nNodes,nNodes)',nodes_sv];
matrix_inp(:,1) = cast(matrix_inp(:,1),"uint8");
disp(matrix_inp)
writematrix(matrix_inp,'M.csv')
 

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

%%
xiTest=[0.1968;-0.4929;0.3093;-0.2883;0.2715;0.4913;0.4962;0.3890;0.0]
% shape-varied mesh 
df = U*xiTest;                       % displacement field introduced by shape variations
dd = [df(1:3:end) df(2:3:end) df(3:3:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);

elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes_sv, elementPlot, 0);
hold off



