% -------------------------------------------------------------------------
% Creation of a figure for the fish position
%
% Last modified: 30/09/2024, Mathieu Dubied, ETH Zurich
% -------------------------------------------------------------------------
clear; 
close all; 
clc
set(groot,'defaulttextinterpreter','latex');

%% PREPARE MODEL   

% load material parameters
load('parameters.mat') 

% specify and create FE mesh
filename ='input_files/3d_rectangle_1272el'; % 24822el' ;%'3d_rectangle_8086el'; 
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions(nodes);
n_elements = size(elements,1);


%% SHAPE VARIATIONS _______________________________________________________
% gather all shape variations
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, x_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% create shape variation basis U
U = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
    z_smallFish, z_notch, x_concaveTail];

% set xi
% xi =  [-0.3584,0.3080,0.3893,0.2906,0.4838,0.2701,0.2271,0.297].';
xi = zeros(8,1);

% plot
f0 = figure('units','centimeters','position',[3 3 10 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:4); hold on 
v1 = reshape(U*xi, 3, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);
axis equal; grid on; box on; 

%% CREATE FIGURE __________________________________________________________

elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
t = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'loose', 'InnerPosition',[0.04, 0.01 ,0.92,0.98]);

v1 = reshape(U*xi, 3, []).';
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);

%%
% Initialize figure
f_abstract = figure('units','centimeters','position',[3 3 9.5 9]);

% Arrow from top to bottom
annotation('arrow', [0.15, 0.15], [0.75, 0.3], 'LineWidth', 2);

% 1. Nominal shape
annotation('textbox', [0.3, 0.9, 0.6, 0.05], 'String', '1. Load nominal shape', ...
               'EdgeColor', 'none', 'FontSize', 12, 'HorizontalAlignment', 'left', 'Interpreter','latex');
axes('Position', [0.02,0.73, 0.26, 0.26]);  
hf = PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
hold off
% annotation('textbox', [0.3, 0.35-0.02,0.6,0.3], 'String','','BackgroundColor', 'yellow')

% axes('Position', [0.3,0.52, 0.7, 0.25]); 
% rectangle('Position',[0,1.0,1.0,0.5])
% axis equal
% axis off

% 2. Shape variations
annotation('textbox', [0.3, 0.8, 0.55, 0.05], 'String', '2. Define shape variations', ...
               'EdgeColor', 'none', 'FontSize', 12,  'Interpreter','latex', ...
               'BackgroundColor','White', 'VerticalAlignment', 'middle');
ax = axes('Position', [0.3,0.57, 0.2, 0.2]);  
v = reshape(U(:,7)*0.9, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
axes('Position', [0.425,0.57, 0.2, 0.2]);
v = reshape(U(:,3)*0.9, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
axes('Position', [0.55,0.57, 0.2, 0.2]);
v = reshape(U(:,8)*0.5, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
axes('Position', [0.675,0.57, 0.2, 0.2]);
v = reshape(U(:,5)*0.9, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
annotation('textbox', [0.79,0.6, 0.6, 0.05], 'String', '...', ...
               'EdgeColor', 'none', 'FontSize', 16, 'Interpreter','latex');
           
axes('Position', [0.33,0.55, 0.62, 0.275]); 
rectangle('Position',[0,0.0,1,1])
xlim([0 1])
ylim([0 1])
axis off

% 3. Optimisation pipeline
annotation('textbox', [0.3, 0.45, 0.44, 0.05], 'String', '3. Run optimisation', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex', ...
               'BackgroundColor','White', 'VerticalAlignment', 'middle');
annotation('textbox', [0.35, 0.4, 0.6, 0.05], 'String', '- Create PROM', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');
annotation('textbox', [0.35, 0.35, 0.6, 0.05], 'String', '- Solve EoMs and sensitivity', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');
annotation('textbox', [0.35, 0.3, 0.6, 0.05], 'String', '- Take gradient step', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');

% box
% annotation('textbox', [0.325, 0.22, 0.65, 0.20], 'String','');
axes('Position', [0.33,0.275, 0.62, 0.2]); 
rectangle('Position',[0,0.0,1,1])
xlim([0 1])
ylim([0 1])
axis off



% 4. Optimal shape
annotation('textbox', [0.3, 0.15, 0.6, 0.05], 'String', '4. Obtain optimal shape', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');
axes('Position', [0.,0.0, 0.28, 0.28]); 
xi =  [-0.3584,0.3080,0.3893,0.2906,0.4838,0.2701,0.2271,0.297].';
v = reshape(U*xi, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', S);
plotcube(L, O, .05, [0 0 0]);

exportgraphics(f_abstract,'figure_abstract.pdf','Resolution',1200)



% L = [Lx,Ly,Lz];
% positions = [0.65, 0.85; 0.65, 0.7; 0.65, 0.1];  % x and y positions for the mesh plots
% 
% for i = [1, 2, 6]
%     axes('Position', [positions(i == [1, 2, 6], 1), positions(i == [1, 2, 6], 2), 0.1, 0.1]);  % Define small axes for each mesh
%     plotcube(L, O, .05, [0 0 0]);  % Plot the cube mesh using the custom plotcube function
% end


%% BUILD ROM AND SOLVE EOMS: Nominal Shape ________________________________
% Parameters
h = 0.01;
tmax = 16.0;
kActu = 1.0;

%%
% Mesh        
df = U*xi;                    % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_sv = nodes + ddf;   % nominal + d ---> defected 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nsetBC)
    svMesh.set_essential_boundary_condition([nsetBC{l}],1:3,0)   
end

% build ROM
fprintf('____________________\n')
fprintf('Building PROM ... \n')
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM_3D(svMesh,nodes_sv,elements,U,USEJULIA,VOLUME,FORMULATION); 

% store dorsal nodes for future use
dorsalNodesStructFromUser.matchedDorsalNodesIdx = spineProperties.dorsalNodeIdx;
dorsalNodesStructFromUser.dorsalNodesElementsVec = spineProperties.dorsalNodesElementVec;
dorsalNodesStructFromUser.matchedDorsalNodesZPos = spineProperties.zPos;

%% 
% Solve EoMs
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax); 
% TI_NL_ROM = solve_EoMs(V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax); 

% Retrieving solutions    
solution = TI_NL_PROM.Solution.u;

%% CHECK SOLUTION _________________________________________________________
uTail = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for a=1:tmax/h
    uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
end

figure
subplot(2,1,1);
plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
xlabel('Time [s]')
ylabel('x-position tail node')
legend('Location','northwest')
subplot(2,1,2);
plot(timePlot,uTail(2,:),'DisplayName','k=0')
xlabel('Time [s]')
ylabel('y-position tail node')
legend('Location','southwest')
drawnow

%% BUILD ROM AND SOLVE EOMS: Optimal Shape ________________________________
xi =  [-0.3584,0.3080,0.3893,0.2906,0.4838,0.2701,0.2271,0.297].';
% Mesh        
df = U*xi;                    % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_sv = nodes + ddf;   % nominal + d ---> defected 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nsetBC)
    svMesh.set_essential_boundary_condition([nsetBC{l}],1:3,0)   
end

% build PROM
fprintf('____________________\n')
fprintf('Building PROM ... \n')
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM_3D(svMesh,nodes_sv,elements,U,USEJULIA,VOLUME,FORMULATION,'dorsalNodes',dorsalNodesStructFromUser); 

%%
% Solve EoMs
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM, ...
    tailProperties,spineProperties,dragProperties,actuTop,actuBottom,...
    kActu,h,tmax); 

% Retrieving solutions    
solution = TI_NL_PROM.Solution.u;

%% CHECK SOLUTION _________________________________________________________
uTail = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for a=1:tmax/h
    uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
end

figure
subplot(2,1,1);
plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
xlabel('Time [s]')
ylabel('x-position tail node')
legend('Location','northwest')
subplot(2,1,2);
plot(timePlot,uTail(2,:),'DisplayName','k=0')
xlabel('Time [s]')
ylabel('y-position tail node')
legend('Location','southwest')
drawnow
%% CREATE SHAPE FIGURE ____________________________________________________
elementPlot = elements(:,1:4); 
nel = size(elements,1);
nDOFperNode = 3;

timesToPlot = [0,5.88,13.99];
scalefactor = 1;
opacityVec = [0.7,0.5,0.3];
colors = [[0.2,0.2,0.2];[0.1,0.1,0.1];[0,0,0]];

nt = size(solution,2);
S = {solution};
ns = 1;

f2 = figure('units','centimeters','position',[3 3 9 4]);

for idx = 1:length(timesToPlot)
    hold on
    timeStep = int16(timesToPlot(idx)/h+1);
    reshapedSol = reshape(solution(:,timeStep),nDOFperNode,[]).';
    yShift = zeros(size(reshapedSol,1),nDOFperNode);
    yShift(:,2) = -ones(size(reshapedSol,1),1);
    displ = reshapedSol(:,1:3) + idx/20*yShift;
    PlotFieldonDeformedMeshSpecificColor(nodes_sv,elements,displ, ...
        brighten([0 0.4470 0.7410],opacityVec(idx)), ...
        'color', colors(idx,:) ,...
        'factor',scalefactor) ;
    
end

% print(f2,'figure1.svg','-dsvg','-r800');

%% -----------












%% CREATE ACTUATION FIGURE
f2 = figure('units','centimeters','position',[3 3 10 2]);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
timePlot = linspace(0,tmax-h,tmax/h);
pPlot = [0.3;2.4;1];
k=300;
actuPlot = zeros(1,length(timePlot));
actuStar = zeros(1,length(timePlot));
for i = 1:length(timePlot)
    actuPlot(i) = actuation_signal_6(k,timePlot(i),pPlot);
end
plot(timePlot,actuPlot)
hold on
grid on
ylabel('Actuation','Interpreter','latex')
xlabel('Time [s]')
hold off

print(f2,'graphical_abstract_actu_optimal_V0.svg','-dsvg','-r800');





