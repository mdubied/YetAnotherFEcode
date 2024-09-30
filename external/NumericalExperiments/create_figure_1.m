% -------------------------------------------------------------------------
% Creation of a figure for the fish position
%
% Last modified: 30/09/2024, Mathieu Dubied, ETH Zurich
% -------------------------------------------------------------------------
clear; 
close all; 
clc

%% PREPARE MODEL   

% load material parameters
load('parameters.mat') 

% specify and create FE mesh
filename ='3d_rectangle_8086el'; % 24822el' ;%'3d_rectangle_8086el'; 
kActu = 1.0;    % multiplicative factor for the actuation forces, dependent on the mesh
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions(nodes);
n_elements = size(elements,1);


%% SHAPE VARIATIONS _______________________________________________________
% gather all shape variations
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% create shape variation basis U
U = [z_tail,y_head];
U = [z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_head,y_linLongTail,y_ellipseFish];

% set xi
% xi = [0.16;-0.4831;0.3;-0.2801;0.2;0.2;0.475;0.45];
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

%% BUILD ROM AND SOLVE EOMS _______________________________________________
% Parameters
h = 0.005;
tmax = 4.0;
pActu = [0.2;2.0;0.3];

% Mesh        
df = U*xi;                    % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_sv = nodes + ddf;   % nominal + d ---> defected 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nset)
    svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
end

% build ROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')

[V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_ROM_3D(svMesh,nodes_sv,elements,USEJULIA);      

% Solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_ROM = solve_EoMs_and_sensitivities_actu(V,ROM_Assembly, ...
    tensors_ROM,tailProperties,spineProperties,dragProperties, ...
    actuTop,actuBottom,h,tmax,pActu);
toc

% Retrieving solutions    
solution = TI_NL_ROM.Solution.u;

%% CHECK SOLUTION
uTail = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for a=1:tmax/h
    uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_ROM.Solution.q(:,a);
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



%% CREATE SHAPE FIGURE
elementPlot = elements(:,1:4); 
nel = size(elements,1);
nDOFperNode = 3;

timesToPlot = [0,1.88,3.99];
scalefactor = 1;
opacityVec = [0.7,0.5,0.3];
colors = [[0.2,0.2,0.2];[0.1,0.1,0.1];[0,0,0]];

nt = size(solution,2);
S = {solution};
ns = 1;

f2 = figure('units','centimeters','position',[3 3 14 4]);

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

print(f2,'figure1.svg','-dsvg','-r800');












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





