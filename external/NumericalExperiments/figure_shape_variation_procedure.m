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
propRigid = 0.6;
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
f1 = figure('units','centimeters','position',[3 3 11 3]);
elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
t = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'loose', 'InnerPosition',[0.04, 0.01 ,0.92,0.98]);

% Shape variation 
ax1 = nexttile(t);
variation = reshape(U*1, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, variation, 'factor', 1);
view(ax1,[0.219559017020941,1.058624487082271]);

% Nominal mesh
ax2 =  nexttile(t);
variationNominal = reshape(U*0, 3, []).';
dm = PlotFieldonDeformedMesh(nodes, elementPlot, variationNominal, 'factor', 1);
view(ax2,[0.219559017020941,1.058624487082271]);
% view(ax2,[-7.138557064977 29.9999999681973]);

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


% Vector field
ax3 =  nexttile(t);
quiver3(vecShapeVariation(1,:),vecShapeVariation(2,:),vecShapeVariation(3,:),...
        vecShapeVariation(4,:),vecShapeVariation(5,:),vecShapeVariation(6,:),...
        0,'r', 'LineWidth', 0.8)
axis equal
axis off;
plotcube(L,O,.05,[0 0 0]);
view(ax3,[0.219559017020941,1.058624487082271]);

% Axis limits
axis([ax1 ax2 ax3],[-0.35 0 -0.045 0.045 -0.1 0.1])

% Add annotation
annotation('textbox',[0.291,0.57,0.1,0.1], 'String','$$-$$', ...
    'FontSize', 20,'Interpreter','latex',...
    'FitBoxToText','on', 'LineStyle','none')
annotation('textbox',[0.625,0.57,0.1,0.1], 'String','$$=$$', ...
    'FontSize', 20,'Interpreter','latex',...
    'FitBoxToText','on', 'LineStyle','none')
annotation('textbox',[0.0,0.07,0.2,0.1], 'String','Shape-varied mesh', ...
    'FontSize', 11,'Interpreter','latex',...
    'FitBoxToText','on', 'LineStyle','none')
annotation('textbox',[0.37,0.07,0.1,0.1], 'String','Nominal mesh', ...
    'FontSize', 11,'Interpreter','latex',...
    'FitBoxToText','on', 'LineStyle','none')
annotation('textbox',[0.725,0.07,0.1,0.1], 'String','Vector field', ...
    'FontSize', 11,'Interpreter','latex',...
    'FitBoxToText','on', 'LineStyle','none')



exportgraphics(f1,'shape_variation_procedure.pdf','Resolution',1200)
% exportgraphics(f1,'shape_variation_procedure.jpg','Resolution',600)