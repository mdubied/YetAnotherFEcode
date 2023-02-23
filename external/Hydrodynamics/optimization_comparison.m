% ------------------------------------------------------------------------ 
% Script to test and compare the optimization pipeline presented in the
% paper.
% 
% Last modified: 19/12/2022, Mathieu Dubied, ETH Zurich
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
E       = 0.6*263824;       %263824;       % Young's modulus [Pa]
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
filename = 'naca0012TRI_medium_mesh';%'naca0012TRI';
[nodes, elements, nset, elset] = mesh_ABAQUSread(filename);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

% boundary conditions: front end of the airfoil fixed through 2 nodes
frontNode = find_node_2D(0,0,nodes);
Lx = 0.15;
backNode = find_node_2D(Lx,0,nodes);
frontNode2 = find_node_2D(0.005,0,nodes);
nset = {frontNode, frontNode2};
MeshNominal.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)


% shape variation
% (1) thinning airfoil 
nodes_translated = [nodes(:,1), nodes(:,2)*0];
vdA = nodes_translated(:,2) - nodes(:,2);
thinAirfoil = zeros(numel(nodes),1);
thinAirfoil(2:2:end) = vdA;
U = thinAirfoil;   % defect basis

%% OPTIMIZATION PIPELINE P1
d = [-1;0];
h = 0.05;
tmax = 0.5;
FOURTHORDER = 0;
tStart = tic;
[xiStar1,xiEvo1,LrEvo1] = optimization_pipeline_1(myElementConstructor,nset,nodes,elements,U,d,h,tmax,...
    FORMULATION,VOLUME,USEJULIA,FOURTHORDER);
tP1 = toc(tStart)

%% OPTIMIZATION PIPELINE P2

tStart = tic;
[xiStar2,xiEvo2,LrEvo2] = optimization_pipeline_2(MeshNominal,nodes,elements,U,d,h,tmax,...
    FORMULATION,VOLUME,USEJULIA,FOURTHORDER);
tP2 = toc(tStart)

%% OPTIMIZATION PIPELINE P3
tStart = tic;
[xiStar3,xiEvo3,LrEvo3] = optimization_pipeline_3(MeshNominal,nodes,elements,U,d,h,tmax,...
    FORMULATION,VOLUME,USEJULIA,FOURTHORDER);
tP3 = toc(tStart)

%% VISUALIZATION

% defected mesh
xi = xiStar1;           % parameter vector
m = length(xi);     % number of parameters

% update defected mesh nodes
d = U*xi;                       % displacement fields introduced by defects
dd = [d(1:2:end) d(2:2:end)]; 
nodes_defected = nodes + dd;    % nominal + d ---> defected 
DefectedMesh = Mesh(nodes_defected);
DefectedMesh.create_elements_table(elements,myElementConstructor);
DefectedMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

figure('units','normalized','position',[.2 .3 .6 .4])
elementPlot = elements(:,1:3); hold on % plot only corners (otherwise it's a mess)
PlotMesh(nodes_defected, elementPlot, 0); 
PlotMesh(nodes,elementPlot,0);
v1 = reshape(U*xi, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
title(sprintf('Defect, \\xi=[%.1f], S=%.1f\\times',...
    xi, S))
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow





