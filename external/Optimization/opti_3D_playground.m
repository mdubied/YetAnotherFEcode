% ------------------------------------------------------------------------ 
% 3D optimization of a fish.
% 
% Last modified: 22/10/2023, Mathieu Dubied, ETH Zurich
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
% for l=1:length(nset)
%     MeshNominal.set_essential_boundary_condition([nset{l}],2,0)   %
%     fixed head
%     MeshNominal.set_essential_boundary_condition([nset{l}],2,0)     % head on "rail"
% end


%% SHAPE VARIATIONS _______________________________________________________

% (1) thinner fish (y direction)
nodes_projected = [nodes(:,1), nodes(:,2)*0, nodes(:,3)];   % projection on x-axis
yDif = nodes_projected(:,2) - nodes(:,2);       % y-difference projection vs nominal
thinFish = zeros(numel(nodes),1);            % create a single long vectors [x1 y1 x2 y2 ...]^T
thinFish(2:3:end) = yDif;                    % fill up all y-positions

% (2) smaller fish (z direction)
nodes_projected = [nodes(:,1), nodes(:,2), nodes(:,3)*0];   % projection on x-axis
zDif = nodes_projected(:,3) - nodes(:,3);       % y-difference projection vs nominal
smallFish = zeros(numel(nodes),1);            % create a single long vectors [x1 y1 x2 y2 ...]^T
smallFish(3:3:end) = zDif;                    % fill up all y-positions

% (3) smaller z tail fish (z direction)
elCenter = [-Lx*0.7;0];
a = Lx*0.3;
b = Lz*0.5;
nodes_ellipse = nodes;
zDif = zeros(size(nodes,1),1);
for n = 1:size(nodes,1) 
    if nodes(n,1) <= elCenter(1)
        nodes_ellipse(n,:) = [nodes(n,1), nodes(n,2), sign(nodes(n,3)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
        zDif(n) = (nodes_ellipse(n,3) - max(nodes(:,3)).*sign(nodes(n,3))).*abs(nodes(n,3))/(Lz*0.5);       % z-difference projection vs nominal
    end
end
fishTailsv = zeros(numel(nodes),1);
fishTailsv(3:3:end) = real(zDif);

% (4) smaller z head fish (z direction)
elCenter = [-Lx*0.3;0];
a = Lx*0.3;
b = Lz*0.5;
nodes_ellipse = nodes;
zDif = zeros(size(nodes,1),1);
for n = 1:size(nodes,1)
    if nodes(n,1) >= elCenter(1)
        nodes_ellipse(n,:) = [nodes(n,1), nodes(n,2), sign(nodes(n,3)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
        zDif(n) = (nodes_ellipse(n,3) - max(nodes(:,3)).*sign(nodes(n,3))).*abs(nodes(n,3))/(Lz*0.5);       % z-difference projection vs nominal
    end
end
fishHeadsv = zeros(numel(nodes),1);
fishHeadsv(3:3:end) = real(zDif);


% shape variations basis
U = [thinFish,fishTailsv,fishHeadsv];    % shape variations basis

% plot the two meshes
xiPlot = [0.23;-0.39;0.1091];
f1 = figure('units','centimeters','position',[3 3 10 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:4); hold on 
% PlotMesh(nodes, elementPlot, 0); 
v1 = reshape(U*xiPlot, 3, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);
axis equal; grid on; box on; 
% set(hf{1},'FaceAlpha',.7); drawnow
set(f1,'PaperUnits','centimeters');
set(f1,'PaperPositionMode','auto');
% set(f1,'PaperSize',[7 3.5]); % Canvas Size
set(f1,'Units','centimeters');
% 
%%
xiTest = [0.23;-0.39;0.1091];
% shape-varied mesh 
df = U*xiTest;                       % displacement field introduced by shape variations
dd = [df(1:3:end) df(2:3:end) df(3:3:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);


%% OPTIMIZATION PARAMETERS
dSwim = [1;0;0]; %swimming direction
h = 0.01;
tmax = 1.5;

%% OPTIMISATION TEST 1 ____________________________________________________
% A = [1 0 0 0;
%     -1 0 0 0;
%     0 1 0 0;
%     0 -1 0 0;
%     0 0 1 0;
%     0 0 -1 0;
%     0 0 0 1;
%     0 0 0 -1];
% b = [0.6;0.6;0.6;0.6;0.6;0.6;0.6;0.6];
A = [1 0 0 ;
    -1 0 0;
    0 1 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
b = [0.4;0.4;0.4;0.4;0.4;0.4];
% A = [1 0;
%     -1 0;
%     0 1;
%     0 -1];
% b = [0.6;0.6;0.6;0.6];
% A = [1;
%     -1];
% b = [0.2;0.2];

tStart = tic;
[xiStar,xiEvo,LrEvo] = optimise_shape_3D(myElementConstructor,nset, ...
    nodes,elements,U,dSwim,h,tmax,A,b,'FORMULATION',FORMULATION,'VOLUME',VOLUME, ...
    'maxIteration',22,'convCrit',0.002,'barrierParam',4,'gStepSize',0.0001,'nRebuild',5);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)


%% VISUALIZATION __________________________________________________________

% plot the two meshes
xiPlot = xiStar;
f1 = figure('units','centimeters','position',[3 3 10 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:4); hold on 
% PlotMesh(nodes, elementPlot, 0); 
v1 = reshape(U*xiPlot, 3, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);
axis equal; grid on; box on; drawnow
set(f1,'PaperUnits','centimeters');
set(f1,'PaperPositionMode','auto');
% set(f1,'PaperSize',[10 3.5]); % Canvas Size
set(f1,'Units','centimeters');
% 

%% PLOT COST FUNCTION OVER ITERATIONS _____________________________________

figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(LrEvo)
grid on
ylabel('$$L_r$$','Interpreter','latex')
xlabel('Iterations')

%% PLOT XI PARAMATER OVER ITERATIONS ______________________________________

figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot3(xiEvo(1,:),xiEvo(2,:),linspace(1,size(xiEvo,2),size(xiEvo,2)));
grid on
xlabel('$$\xi_1$$','Interpreter','latex')
ylabel('$$\xi_2$$','Interpreter','latex')
zlabel('Iterations')

%%
figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(xiEvo(2,:));
grid on
xlabel('$$\xi_1$$','Interpreter','latex')
ylabel('Iterations')


%%
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

TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

% ROM TENSORS - ACTUATION FORCES __________________________________________

elementPlot = elements(:,1:3); 
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


actuationValues = zeros(size(TI_NL_PROM.Solution.u,2),1);
for t=1:size(TI_NL_PROM.Solution.u,2)
    actuationValues(t) = 1+0.005*sin(t*h*2*pi);
end

actuationValues2 = zeros(size(TI_NL_PROM.Solution.u,2),1);
for t=1:size(TI_NL_PROM.Solution.u,2)
    actuationValues2(t) = 1-0.005*sin(t*h*2*pi);
end

AnimateFieldonDeformedMeshActuation2Muscles(nodes_sv, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,TI_NL_PROM.Solution.u, ...
    'factor',1,'index',1:2,'filename','result_video','framerate',1/h)

%% COST COMPUTATION ON INITIAL ROM ________________________________________
% % 88.57
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
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax,...
%                     'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
% 
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

% %% PLOTs __________________________________________________________________
% 
% % find a specific result node and corresponding DOF
% tailNodeDOFS = MeshNominal.get_DOF_from_location([Lx, 0]);
% tailNodeDOF = tailNodeDOFS(1); % y-direction
% % time axis
% tplot=linspace(0,tmax,tmax/h+1);
% 
% % plot
% figure('units','normalized','position',[.1 .1 .8 .6],'name','Vertical displacement of the tail node')
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% set(groot,'defaultAxesTickLabelInterpreter','latex'); 
% 
% plot(tplot,TI_NL_PROM.Solution.u(tailNodeDOF,1:end-1)*100, "--")
% 
% ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
% xlabel('Time [s]')
% set(gca,'FontName','ComputerModern');
% grid on
% %legend({'FOM','FOM-t','ROM','PROM'}, 'Location', 'eastoutside','Orientation','vertical')
% hold off
