% ------------------------------------------------------------------------ 
% 3D shape optimisation of a fish.
% 
% Last modified: 16/01/2025, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

%% PREPARE MODELS _________________________________________________________                                                   

% load material parameters
load('parameters.mat') 

% specify and create FE mesh
filename ='InputFiles/3d_rectangle_4270el'; % 24822el' ;%'3d_rectangle_8086el';


%'3d_rectangle_1272el';%'3d_rectangle_1272el';%'3d_rectangle_660el'; 
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
n_elements = size(elements,1);

% plot nominal mesh
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off


%% SHAPE VARIATIONS _______________________________________________________

[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% SO1
U = [z_tail,z_head,y_thinFish];
% U = [xz_concaveTail];

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
S = 1.0;
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
h = 0.02;
tmax = 2.0;
kActu = 4.0e5;    % multiplicative factor for the actuation forces

%% OPTIMISATION SO1 _______________________________________________________
% Typical results:  -0.4952, 0.3272, 0.2974

% Preparation (shape variation and opt. constraints)
U_1 = [z_tail,z_head,y_thinFish]; 
A = [1 0 0 ;
    -1 0 0;
    0 1 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
b = [0.5;0.5;0.5;0.5;0.3;0.3];
barrierParam = 10*ones(1,length(b));

% Optimisation
tStart = tic;
[xiStar_1,xiEvo_1,LEvo_1, LwoBEvo_1, nIt_1] = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_1,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',35, ...
    'convCrit',0.004, ...
    'convCritCost',1.0, ... % 1.0 for 8086, 0.3 below, 5.0 for 24822
    'barrierParam',barrierParam, ...
    'gStepSize',0.0003, ...  % 0.0005 for 8086, 0.002 below, 0.0001 for 24822
    'nRebuild',6, ...
    'rebuildThreshold',0.15,...
    'USEJULIA',1);
tOpti_1 = toc(tStart)/60;   % unit is minute
tPerIt_1 = tOpti_1/nIt_1;

% Print stats
fprintf('Computation time: %.2fmin\n',tOpti_1)
fprintf('Number of built models and solved EoMs: %5d\n',nIt_1)
fprintf('Computation time per models/EoMs: %.2f\n',tPerIt_1)

% Save results
filename = sprintf('Results/Data/SO1_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
save(filename,'xiStar_1','xiEvo_1','LEvo_1','LwoBEvo_1','tOpti_1','nIt_1','tPerIt_1')

%% OPTIMISATION SO2 _______________________________________________________
% Typical results: -0.4767, 0.3332, 0.3694, 0.3752, 0.4036

% Preparation (shape variation and opt. constraints)
U_2 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
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
    0.5;0.5;    % set to 0.4 for 24822, otherwise 0.5
    0.8;0.8;
    0.8;0.8];
barrierParam = ones(1,length(b));

% Optimisation
tStart = tic;
[xiStar_2,xiEvo_2,LEvo_2, LwoBEvo_2,nIt_2] = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_2,h,tmax,A,b,...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ...
    'convCrit',0.01, ...   % set to 0.01 for 8086el, 0.004 below
    'convCritCost',0.5, ... % set to 0.5 for 4270 el, 0.1 below, 5 for 24822
    'barrierParam',barrierParam, ...
    'gStepSize',0.0005,...   % set to 0.001 for 4270 el, 0.003 below, 0.0005 for 8086, 0.0002 for 24822
    'nRebuild',5, ...
    'rebuildThreshold',0.15,...
    'USEJULIA',1);
tOpti_2 = toc(tStart)/60;   % unit: minutes
tPerIt_2 = tOpti_2/nIt_2;  

% Print stats
fprintf('Computation time: %.2fmin\n',tOpti_2)
fprintf('Number of built models and solved EoMs: %5d\n',nIt_2)
fprintf('Computation time per models/EoMs: %.2f\n',tPerIt_2)
%%
% Save results
filename = sprintf('Results/Data/SO2_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
save(filename,'xiStar_2','xiEvo_2','LEvo_2','LwoBEvo_2','tOpti_2','nIt_2','tPerIt_2')

%% OPTIMISATION SO3 _______________________________________________________
% Typical results: -0.3584,0.3080,0.3893,0.2906,0.4838,0.2701,0.2271,0.2975

% Preparation (shape variation and opt. constraints)
U_3 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
    z_smallFish, z_notch, xz_concaveTail];
nParam = 8;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
b = [0.4;0.4;
    0.4;0.4;
    0.4;0.4;
    0.3;0.3;
    0.4;0.4;    % set to 0.4 for 24822 elements, 0.5 else
    0.3;0.3;
    0.3;0.3;
    0.2; 0.01];  % concave tail only in one direction 
barrierParam = 3*ones(1,length(b));

% Optimisation
tStart = tic;
[xiStar_3,xiEvo_3,LEvo_3, LwoBEvo_3,nIt_3] = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_3,h,tmax,A,b, ...
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
tOpti_3 = toc(tStart)/60; % unit: minutes
tPerIt_3 = tOpti_3/nIt_3;
xiStar = xiStar_3;  

% Print stats
fprintf('Computation time: %.2fmin\n',tOpti_3)
fprintf('Number of built models and solved EoMs: %5d\n',nIt_3)
fprintf('Computation time per models/EoMs: %.2f\n',tPerIt_3)

% Save results
filename = sprintf('Results/Data/SO3_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
save(filename,'xiStar_3','xiEvo_3','LEvo_3','LwoBEvo_3','tOpti_3','nIt_3','tPerIt_3')


%% PLOT OPTIMAL SHAPE _____________________________________________________

SOIdx = 3;
switch SOIdx
    case 1
        xiStar = xiStar_1;
        U = U_1;   
    case 2
        xiStar = xiStar_2;
        U = U_2;
    case 3
        xiStar = xiStar_3;
        U = U_3;
end

f_opt_shape = figure('units','centimeters','position',[3 3 9 4]);
elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];

t = tiledlayout(1, 2, 'TileSpacing', 'None', 'Padding', 'tight', 'InnerPosition',[0.01, 0.01 ,0.98,0.98]);
v = reshape(U*xiStar, 3, []).';

% View 1
ax1 = nexttile(t);
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);

% View 2
ax2 = nexttile(t);
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
view(ax2,[23.28863587057058,35.063216347984103])

axis([ax1 ax2],[-0.23 0 -0.02 0.02 -0.1 0.1])

fig_title = sprintf('Results/Figures/SO%d_opt_shape_%d_el_kActu_%.3f', SOIdx, n_elements, kActu);
exportgraphics(f_opt_shape,strcat(fig_title,'.pdf'),'Resolution',1200)

%% PLOT COST FUNCTION FOR THE 3 EXPERIMENTS _______________________________
f_cost = figure('units','centimeters','position',[3 3 9 6]);
hold on
plot(LwoBEvo_1,'LineWidth',1)
plot(LwoBEvo_2,'--','LineWidth',1)
plot(LwoBEvo_3,'-.','LineWidth',1)
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('SO1','SO2','SO3')
hold off

exportgraphics(f_cost,'results/figures/SO_cost_evolution.pdf','Resolution',1200)
% exportgraphics(f_cost,'results/figures/SO_cost_evolution.jpg','Resolution',600)



%% OTHER/PREVIOUS PLOTS ___________________________________________________








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
