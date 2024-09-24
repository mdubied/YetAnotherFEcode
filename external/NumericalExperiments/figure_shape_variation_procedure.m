% ------------------------------------------------------------------------ 
% Creation of a figure showing the procedure to create a shape variation
% 
% Last modified: 13/09/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc
set(groot,'defaulttextinterpreter','latex');

%% PREPARE NOMINAL MESH ___________________________________________________                                                   

% load material parameters
load('parameters.mat') 

% specify and create FE mesh
filename ='3d_rectangle_8086el' ;%'3d_rectangle_8086el'; %'3d_rectangle_8086el'
%'3d_rectangle_1272el';%'3d_rectangle_1272el';%'3d_rectangle_660el'; 
kActu = 1.0;    % multiplicative factor for the actuation forces, dependent on the mesh
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions(nodes);


% plot nominal mesh [Can be removed]
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off

%% DEFINE SHAPE VARIATIONS ________________________________________________
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, x_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

U = z_notch;

%% CREATE FIGURE __________________________________________________________
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


% shape variation 
ax1 = subplot(2,2,1,'Position',pos1);
variation = reshape(U*1, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, variation, 'factor', 1);
view(ax1,[-7.138557064977 29.9999999681973]);

% nominal mesh
ax2 = subplot(2,2,2,'Position',pos2);
variationNominal = reshape(U*0, 3, []).';
dm = PlotFieldonDeformedMesh(nodes, elementPlot, variationNominal, 'factor', 1);
view(ax2,[-7.138557064977 29.9999999681973]);

% compute vectors to show variation field
vecShapeVariation = [];
UReshaped = reshape(U*1, 3, []).';
for n=1:size(nodes,1)
    if (nodes(n,3)==Lz/2 || nodes(n,3)==-Lz/2) && (nodes(n,2)==0 ||nodes(n,2)==Ly/2 || nodes(n,2)==-Ly/2)
        dif = UReshaped(n,:);
        vecForPlot = [nodes(n,:),dif].';
        vecShapeVariation = [vecShapeVariation,vecForPlot];
    end
end



% shape variation 3
ax3 = subplot(2,2,3,'Position',pos3);
quiver3(vecShapeVariation(1,:),vecShapeVariation(2,:),vecShapeVariation(3,:),...
        vecShapeVariation(4,:),vecShapeVariation(5,:),vecShapeVariation(6,:),...
        0,'r', 'LineWidth', 0.8)
axis equal
axis off;
plotcube(L,O,.05,[0 0 0]);
view(ax3,[-7.138557064977 29.9999999681973]);



% optimal shape
ax4 = subplot(2,2,4,'Position',pos4);
xiPlot = 0.5;% [-0.3;0.3;0.1;0.4;0.4;0.1;0.2;0.1];%xiStar;
% [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
%     z_smallFish, z_notch, x_concaveTail];
% xiPlot(1) = 0.2;
% xiPlotName = strcat('[',num2str(xiStar(1)),', ',num2str(xiStar(2)),', ',num2str(xiStar(3)),']^\top$$');

v1 = reshape(U*xiPlot, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', 1);
plotcube(L,O,.05,[0 0 0]);
subplotName = '$$\mathbf{\xi}^\ast$$';
text(textPosX, textPosY, textPosZ, subplotName,'Interpreter','latex')


annotation('textbox',[0,0,0.5,0.5],'String','Variation ()','EdgeColor','none','Interpreter','latex');
annotation('arrow',[0.5,0.5],[0.5,0.0])

axis([ax1 ax2 ax3 ax4],[-0.40 0 -0.05 0.05 -0.16 0.16])
% set(ax2, 'box', 'on', 'Visible', 'on')
% set(ax1, 'box', 'on', 'Visible', 'on')

% exportgraphics(f1,'SO3_shapes_V1.pdf','Resolution',1200)
% exportgraphics(f1,'SO3_shapes_V1.jpg','Resolution',600)