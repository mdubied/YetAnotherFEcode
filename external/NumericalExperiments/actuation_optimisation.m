% ------------------------------------------------------------------------ 
% 3D actuation signal optimisation of a fish.
% 
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 0;

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
end
% MESH ____________________________________________________________________

% nominal mesh
filename = '3d_fish_for_mike';%'3d_rectangle_660el';%'fish3_664el';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);

% nodes = nodes*0.01;
% nodes(:,2)=0.8*nodes(:,2);

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
% for l=1:length(nset)
%     MeshNominal.set_essential_boundary_condition([nset{l}],2,0)   %
%     fixed head
%     MeshNominal.set_essential_boundary_condition([nset{l}],2,0)     % head on "rail"
% end




%% OPTIMIZATION PARAMETERS
dSwim = [1;0;0]; %swimming direction
h = 0.005;
tmax = 2.0;

%% OPTIMISATION ___________________________________________________________
% The actuation force needs to be chosen in the following script:
%   solve_EoMs_and_sensitivities_actu
% The initial conditions should be changed in the following script:
%   optimise_actuation_3D

% AO 1
A = [1 0 0 0;
    -1 0 0 0;
    0 1 0 0;
    0 -1 0 0;
    0 0 1 0;
    0 0 -1 0;
    0 0 0 1;
    0 0 0 -1];
b = [1.15;-0.75;0.4;0.4;0.2;0.2;1.3;-0.7];

tStart = tic;
[xiStar,xiEvo,LEvo, ~] = optimise_actuation_3D(myElementConstructor,nset, ...
    nodes,elements,dSwim,h,tmax,A,b, ...
    'maxIteration',30, ...
    'convCrit',0.0005, ...
    'convCritCost',0.002, ...
    'barrierParam',600, ...
    'gStepSize',0.05, ...
    'nResolve',8, ...
    'resolveThreshold',0.1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)

%%
AActu = [1 0 0;
        -1 0 0;
        0 1 0;
        0 -1 0;
        0 0 1;
        0 0 -1];
bActu = [0.25;-0.15;2.6;-1.4;1.05;-0.05];

tStart = tic;
[xiStar,xiEvo,LEvo, ~] = optimise_actuation_3D(myElementConstructor,nset, ...
    nodes,elements,dSwim,h,tmax,AActu,bActu, ...
    'maxIteration',30, ...
    'convCrit',0.0005, ...
    'convCritCost',0.002, ...
    'barrierParam',600, ...
    'gStepSize',0.001, ...
    'nResolve',8, ...
    'resolveThreshold',0.1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)

%%
% AO 2 (Change the acutation functions if you change the time step h)
% nParam = int32(tmax/h*0.1+2);
% A = zeros(2 * nParam, nParam);
% for i = 1:nParam
%     A(2*i-1:2*i,i) =[1;-1];
% end
% b = 1.3*ones(nParam*2,1);
% 
% tStart = tic;
% [xiStar,xiEvo,LEvo, ~] = optimise_actuation_3D(myElementConstructor,nset, ...
%     nodes,elements,dSwim,h,tmax,A,b, ...
%     'maxIteration',30, ...
%     'convCrit',0.0005, ...
%     'convCritCost',0.01, ...
%     'barrierParam',10, ...
%     'gStepSize',0.04, ...
%     'nResolve',8, ...
%     'resolveThreshold',0.1);
% topti = toc(tStart);
% fprintf('Computation time: %.2fmin\n',topti/60)


%% PLOT COST FUNCTION OVER ITERATIONS _____________________________________

fCost = figure;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(LEvo)
grid on
ylabel('$$L$$','Interpreter','latex')
xlabel('Iterations')
%%
exportgraphics(fCost,'AO1_cost_V0.pdf','Resolution',600)


%% PLOT PARAMETERS' EVOLUTION OVER ITERATIONS _____________________________
figure('Position',[100,100,600,200])
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
subplot(1,3,1)
plot(xiEvo(1,:));
ylabel('$$p_1$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,2)
plot(xiEvo(2,:));
ylabel('$$p_2$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,3)
plot(xiEvo(3,:));
ylabel('$$p_3$$','Interpreter','latex')
xlabel('Iterations')
grid on

%% PLOT ACTUATION SIGNAL _____________________________________

f1 = figure('units','centimeters','position',[3 3 9 4]);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
timePlot = linspace(0,tmax-h,tmax/h);
p0 = [1;0;0;1];
k=300;
actu0 = zeros(1,length(timePlot));
actuStar = zeros(1,length(timePlot));
for i = 1:length(timePlot)
    actu0(i) = actuation_signal_4(k,timePlot(i),p0);
    actuStar(i) = actuation_signal_4(k,timePlot(i),xiStar);
end
plot(timePlot,actu0,'--')
hold on
plot(timePlot,actuStar)
grid on
ylabel('Actuation signal','Interpreter','latex')
xlabel('Time [s]')
legend('Initial','Optimised','Interpreter','latex')%Location='northoutside',Orientation='horizontal'
hold off
exportgraphics(f1,'AO1_signal_V0.pdf','Resolution',600)
