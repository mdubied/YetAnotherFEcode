% ------------------------------------------------------------------------ 
% Script to test the implementation of hydrodynamic forces on 2D structures
% with TRI3 elements
% 
% Last modified: 21/10/2022, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

%whichModel = 'CUSTOM'; % or "ABAQUS"
whichModel = 'ABAQUS';
%elementType = 'QUAD4';
elementType = 'TRI3';
% elementType = 'QUAD8'; % only QUAD4 is implemented for now

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
nx = 5;
ny = 6;
switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,elementType);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'naca0012TRI';
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
end

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

% boundary conditions of nominal mesh: front and back nodes fixed
frontNode = find_node_2D(0,0,nodes);
Lx = 0.15;
backNode = find_node_2D(Lx,0,nodes);
nset = {frontNode,backNode};
MeshNominal.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)


% ASSEMBLY ________________________________________________________________
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);

% store matrices
NominalAssembly.DATA.K = Kn;
NominalAssembly.DATA.M = Mn;



%% PLOT MESH WITH NODES AND ELEMENTS
elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 1);

%% TESTING BASIC FORCE ORIENTATION (PRELIMINARY STUDY - NOT NEEDED)
[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);
%f = BeamAssembly.vector('drag_force', u, ud);
vwater = [1;0];
rho = 1;
%f = AirfoilAssembly.skin_force('force_length_prop_skin_normal', 'weights',
%skinElements, skinElementFaces); % change in assembly -> skin_force needed
%for this function
f = NominalAssembly.skin_force('drag_force', 'weights', skinElements, skinElementFaces, vwater, rho);
% plot forces
PlotMeshandForce(nodes, elementPlot, 1,30*f);

%% HYDRODYNAMIC FORCES - TENSORIAL APPROACH
[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);

vwater = [1;0.1];
rho = 1;
tensorsHydroFOM = tensors_hydro_FOM(NominalAssembly, elements, skinElements, skinElementFaces, vwater, rho);




