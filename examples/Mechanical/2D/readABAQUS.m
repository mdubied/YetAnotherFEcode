% Test how to read meshes coming from abaqus inp files
clear; 
close all; 
clc

whichModel = 'ABAQUS'; 
elementType = 'TRI3';
%elementType = 'QUAD4';


%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
switch elementType
    case 'QUAD4'
        myElementConstructor = @()Quad4Element(thickness, myMaterial);
    case 'QUAD8'
        myElementConstructor = @()Quad8Element(thickness, myMaterial);
    case 'TRI3'
        myElementConstructor = @()Tri3Element(thickness, myMaterial);
end

% MESH_____________________________________________________________________
Lx = 3;
Ly = .2;
nx = 30;
ny = 3;
switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,elementType);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'naca0012TRI';  % triangle mesh
        %filename = 'naca0012QUAD';  % quad mesh
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
end

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

elementPlot = elements(:,1:3); % for TRI elements
PlotMesh(nodes, elementPlot, 0);
