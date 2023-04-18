% ------------------------------------------------------------------------ 
% 2D optimization of a fish: plot shape variations.
%
% Things to change in PlotMesh:
%           p = patch(X,Y,'w','DisplayName','Mesh','LineWidth',0.01,'EdgeAlpha',0.2); 
%       else % line
%           %p=plot(Nodes(:,1),Nodes(:,2),'.-k', 'Markersize',10,'LineWidth',0.1);% no linewidth
%
% Things to change in PlotFieldonDeformedMesh
%           h{1} = patch(defoX,defoY,profile,'EdgeColor',meshcolor,'EdgeAlpha',0.2);
%           %h{2} = plot(defoX(:),defoY(:),'.','Color', meshcolor, 'Markersize',5); %5 %10 
%       else
%           %h = plot(Nodes(:,1)+factor*ux,Nodes(:,2)+factor*uy,'.-','Color', meshcolor, 'Markersize',5); %10
% 
% Last modified: 17/04/2023, Mathieu Dubied, ETH Zurich
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

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 260000;   % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 
thickness = .1;         % [m] out-of-plane thickness

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain

% element
switch elementType
    case 'TRI3'
        myElementConstructor = @()Tri3Element(thickness, myMaterial);
end

% PROM parameters
xi1 = 0.1;  % thinning fish
xi2 = -0.2; % shortening fish
xi3 = 0.6;  % linear tail thinning
xi4 = 0.3;  % long ellipse tail thinning
xi5 = -0.1;  % short ellipse tail thinning
xi6 = 0;    % linear head thinning 
xi7 = 0.0;  % long ellipse head thinning
xi8 = 0.5;  % short ellipse head thinning
xi9 = 0.1;  % tail smoothening
xi10 = 0.3; % head smoothening

% xi1 = 0.2;  % thinning fish
% xi2 = 0.2; % shortening fish
% xi3 = -0.2;  % linear tail thinning
% xi4 = 0.4;  % long ellipse tail thinning
% xi5 = -0.5;  % short ellipse tail thinning
% xi6 = -0.2;    % linear head thinning 
% xi7 = 0.3;  % long ellipse head thinning
% xi8 = 0.1;  % short ellipse head thinning
% xi9 = -0.2;  % tail smoothening
% xi10 = 0.1; % head smoothening




% MESH ____________________________________________________________________

% nominal mesh
switch upper( whichModel )
    case 'ABAQUS'
        filename = 'fishNominalTRI3_420El';%'rectangle192Elements';%'fishNominalTRI3';%'fishNominalTRI3_420El';
        [nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
end

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil

% plot nominal mesh
% elementPlot = elements(:,1:3); 
% figure('units','normalized','position',[.2 .1 .6 .4],'name','Nominal mesh with element and node indexes')
% PlotMesh(nodes, elementPlot, 0);

% boundary conditions of nominal mesh
nel = size(elements,1);
nset = {};
for el=1:nel   
    for n=1:size(elements,2)
        if  nodes(elements(el,n),1)<Lx*0.15 && ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
    end   
end
for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],1:2,0)
end

%% PLOTS
% shape variations 
[thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead,smoothTail,smoothHead] = shape_variations_2D(nodes,Lx,Ly);
Utot = [thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead,smoothTail,smoothHead];
for i=1:10
    U = Utot(:,i);
    if i==1 || i==2
        xi = 0.5;
    else
        xi = 1;
    end
 
    % positive direction
    % shape-varied mesh 
    d = U*xi;                       % displacement field introduced by shape variations
    dd = [d(1:2:end) d(2:2:end)];   % rearrange as two columns matrix
    nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
    svMesh = Mesh(nodes_sv);
    svMesh.create_elements_table(elements,myElementConstructor);
    svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
    
    % plot the two meshes
    f1 = figure('units','centimeters','position',[3 3 10 3.5],'name','Shape-varied mesh');
    elementPlot = elements(:,1:3); hold on 
    %PlotMesh(nodes_sv, elementPlot, 0); 
    PlotMesh(nodes,elementPlot,0);
    v1 = reshape(U*xi, 2, []).';
    S = 1;
    hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
    axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow
    set(f1,'PaperUnits','centimeters');
    set(f1,'PaperPositionMode','auto');
    set(f1,'PaperSize',[10 3.5]); % Canvas Size
    set(f1,'Units','centimeters');
    
    % save image as pdf
    f1=gcf;
    name1 = "sv";
    name2 = "1.pdf";
    caseID = num2str(i);
    str = [name1,caseID,name2];
    fileName = join(str,"_");
    saveas(gcf,fileName); %Set desired file name

    % negative direction
    xi = -xi;
    % shape-varied mesh 
    d = U*xi;                       % displacement field introduced by shape variations
    dd = [d(1:2:end) d(2:2:end)];   % rearrange as two columns matrix
    nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
    svMesh = Mesh(nodes_sv);
    svMesh.create_elements_table(elements,myElementConstructor);
    svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
    
    % plot the two meshes
    f2 = figure('units','centimeters','position',[3 3 10 3.5],'name','Shape-varied mesh');
    elementPlot = elements(:,1:3); hold on 
    PlotMesh(nodes_sv, elementPlot, 0); 
    PlotMesh(nodes,elementPlot,0);
    v1 = reshape(U*xi, 2, []).';
    S = 1;
    hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
    axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow
    set(f2,'PaperUnits','centimeters');
    set(f2,'PaperPositionMode','auto');
    set(f2,'PaperSize',[10 3.5]); % Canvas Size
    set(f2,'Units','centimeters');
    
    
    % save image as pdf
    f2=gcf;
    name1 = "sv";
    name2 = "2.pdf";
    caseID = num2str(i);
    str = [name1,caseID,name2];
    fileName = join(str,"_");
    saveas(gcf,fileName); %Set desired file name
end
