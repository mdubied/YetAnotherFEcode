% ------------------------------------------------------------------------ 
% Creation of a figure containing all shape variation used
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

U = [z_tail,z_head, y_thinFish,z_smallFish,z_notch, y_head, y_linLongTail,y_ellipseFish, y_ellipseFish];

%% PLOT SHAPE VARIATIONS __________________________________________________
f1 = figure('units','centimeters','position',[3 0 18 19]);
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
textPos = [textPosX, textPosY, textPosZ];

n_subplot = 2*9;

axesHandles = gobjects(1, 14); % Preallocate an array for the axes handles

for i = 1:9
    % Shape variation i, positive
    axesHandles((i-1)*2 + 1) = subplot(5, 4, (i-1)*2 + 1);
    create_subplot(U(:,i), 0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Shape variation i, negative
    axesHandles((i-1)*2 + 2) = subplot(5, 4, (i-1)*2 + 2);
    create_subplot(U(:,i), -0.5, num2str(i), nodes, elementPlot, L, O, textPos);
end

% Set axis limits for all the axes
axis(axesHandles, [-0.40 0 -0.061 0.061 -0.16 0.16]);

%%
%% PLOT SHAPE VARIATIONS __________________________________________________
f1 = figure('units','centimeters','position',[3 0 18 19]);
elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
textPosX = -0.2;
textPosY = -0.05;
textPosZ = -0.14;
textPos = [textPosX, textPosY, textPosZ];

% Define tiled layout with tighter spacing
t = tiledlayout(5, 4, 'TileSpacing', 'none', 'Padding', 'none');

n_subplot = 2 * 9;

axesHandles = gobjects(1, n_subplot); % Preallocate an array for the axes handles

for i = 1:9
    % Shape variation i, positive
    axesHandles((i-1)*2 + 1) = nexttile(t, (i-1)*2 + 1);
    create_subplot(U(:,i), 0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Shape variation i, negative
    axesHandles((i-1)*2 + 2) = nexttile(t, (i-1)*2 + 2);
    create_subplot(U(:,i), -0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Get the positions of the two subplots
    pos1 = get(axesHandles((i-1)*2 + 1), 'Position');
    pos2 = get(axesHandles((i-1)*2 + 2), 'Position');
    
    % Calculate the rectangle position that encloses both subplots
    x = min(pos1(1), pos2(1)); % Minimum x position
    y = min(pos1(2), pos2(2)); % Minimum y position
    width = max(pos1(1) + pos1(3), pos2(1) + pos2(3)) - x; % Total width
    height = max(pos1(2) + pos1(4), pos2(2) + pos2(4)) - y; % Total height
    
    % Draw a rectangle around the two subplots
    annotation('rectangle', [x, y, width, height], ...
        'Color', 'k', 'LineWidth', 1.0, 'LineStyle', '--');
end

% Set axis limits for all the axes
axis(axesHandles, [-0.40 0 -0.061 0.061 -0.16 0.16]);


%%
%% PLOT SHAPE VARIATIONS AND OPTIMAL SHAPE ________________________________
f1 = figure('units','centimeters','position',[3 0 18 19], 'Visible', 'off');  % Create invisible figure to suppress default axes
elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
textPosX = -0.2;
textPosY = -0.05;
textPosZ = -0.14;
textPos = [textPosX, textPosY, textPosZ];

% Define grid parameters
nRows = 5;   % Number of rows in the layout (for 9 variations, 5 rows is more than enough)
nCols = 4;   % Number of columns (positive/negative for each variation)
subplot_width = 0.2;  % Width of each subplot in normalized units
subplot_height = 0.15; % Height of each subplot in normalized units
x_offset = 0.05;      % Horizontal offset between subplots
y_offset = 0.05;      % Vertical offset between rows

axesHandles = gobjects(1, 18); % Preallocate an array for 18 axes handles (9 pairs)

for i = 1:9
    % Calculate row and column for each pair of subplots
    row = floor((i-1) / 2); % Determine row number (0 to nRows-1)
    col1 = mod((i-1) * 2, 4); % First column position for positive shape variation
    col2 = col1 + 1; % Second column position for negative shape variation
    
    % Calculate subplot positions (positive and negative variations)
    pos1 = [x_offset + col1*(subplot_width + 0.02), 1 - (row + 1)*(subplot_height + y_offset), subplot_width, subplot_height];
    pos2 = [x_offset + col2*(subplot_width + 0.02), 1 - (row + 1)*(subplot_height + y_offset), subplot_width, subplot_height];
    
    % Shape variation i, positive
    axesHandles((i-1)*2 + 1) = axes('Position', pos1);
    create_subplot(U(:,i), 0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Shape variation i, negative
    axesHandles((i-1)*2 + 2) = axes('Position', pos2);
    create_subplot(U(:,i), -0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Draw a rectangle around the two subplots
    rect_x = pos1(1);
    rect_y = pos1(2);
    rect_w = pos2(1) - pos1(1) + subplot_width;
    rect_h = subplot_height;
    annotation('rectangle', [rect_x, rect_y, rect_w, rect_h], 'Color', 'k', 'LineWidth', 1.0, 'LineStyle', '--');
end

% Set axis limits for all the axes
axis(axesHandles, [-0.40 0 -0.061 0.061 -0.16 0.16]);

% Make the figure visible after all axes are set
set(f1, 'Visible', 'on');

%%

%% PLOT SHAPE VARIATIONS AND OPTIMAL SHAPE ________________________________
f1 = figure('units','centimeters','position',[3 0 18 19]);

% Create a tiled layout without any padding or extra space
t = tiledlayout(5, 4, 'TileSpacing', 'tight', 'Padding', 'tight');

elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
textPosX = -0.2;
textPosY = -0.05;
textPosZ = -0.14;
textPos = [textPosX, textPosY, textPosZ];

axesHandles = gobjects(1, 18); % Preallocate an array for 18 axes handles (9 pairs)

for i = 1:9
    % Shape variation i, positive
    axesHandles((i-1)*2 + 1) = nexttile(t);  % Assign the next tile in the layout
    create_subplot(U(:,i), 0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Shape variation i, negative
    axesHandles((i-1)*2 + 2) = nexttile(t);  % Assign the next tile in the layout
    create_subplot(U(:,i), -0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Draw a rectangle around the two subplots
    pos1 = get(axesHandles((i-1)*2 + 1), 'Position');
    pos2 = get(axesHandles((i-1)*2 + 2), 'Position');
    
    % Get the position for the bounding rectangle
    rect_x = min(pos1(1), pos2(1));
    rect_y = min(pos1(2), pos2(2));
    rect_w = max(pos1(1) + pos1(3), pos2(1) + pos2(3)) - rect_x;
    rect_h = max(pos1(2) + pos1(4), pos2(2) + pos2(4)) - rect_y;
    
    annotation('rectangle', [rect_x, rect_y, rect_w, rect_h], 'Color', 'k', 'LineWidth', 1.0, 'LineStyle', '--');
end

% Set axis limits for all the axes
axis(axesHandles, [-0.40 0 -0.061 0.061 -0.16 0.16]);

% Set figure visibility to 'on' (in case it's off)
set(f1, 'Visible', 'on');

%%
%% PLOT SHAPE VARIATIONS AND OPTIMAL SHAPE ________________________________
%% PLOT SHAPE VARIATIONS AND OPTIMAL SHAPE ________________________________
f1 = figure('units','centimeters','position',[3 0 18 19]);

elementPlot = elements(:,1:4); hold on
L = [Lx, Ly, Lz];
O = [-Lx, -Ly/2, -Lz/2];
textPosX = -0.2;
textPosY = -0.05;
textPosZ = -0.14;
textPos = [textPosX, textPosY, textPosZ];

nRows = 5;
nCols = 4;
n_subplot = 2 * 9;

axesHandles = gobjects(1, n_subplot); % Preallocate an array for the axes handles

% Set the desired spacing
horizontalSpacing = -0.00; % Horizontal space between subplots
verticalSpacing = -0.00;   % Vertical space between subplots

% Loop through and create subplots
for i = 1:9
    % Shape variation i, positive
    axesHandles((i-1)*2 + 1) = subplot(nRows, nCols, (i-1)*2 + 1);
    create_subplot(U(:,i), 0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Shape variation i, negative
    axesHandles((i-1)*2 + 2) = subplot(nRows, nCols, (i-1)*2 + 2);
    create_subplot(U(:,i), -0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Get the position of the current subplots and adjust for spacing
    pos1 = get(axesHandles((i-1)*2 + 1), 'Position');
    pos2 = get(axesHandles((i-1)*2 + 2), 'Position');
    
    % Adjust the positions to reduce spacing
    pos1(3) = pos1(3) - horizontalSpacing;  % Reduce width (to create space)
    pos1(4) = pos1(4) - verticalSpacing;    % Reduce height (to create space)
    pos2(3) = pos2(3) - horizontalSpacing;
    pos2(4) = pos2(4) - verticalSpacing;
    
    set(axesHandles((i-1)*2 + 1), 'box', 'on', 'Visible', 'on')
    set(axesHandles((i-1)*2 + 2), 'box', 'on', 'Visible', 'on')

%     % Apply new positions
%     set(axesHandles((i-1)*2 + 1), 'Position', pos1);
%     set(axesHandles((i-1)*2 + 2), 'Position', pos2);
    
    % Draw a rectangle around the two subplots
    rect_x = min(pos1(1), pos2(1));
    rect_y = min(pos1(2), pos2(2));
    rect_w = max(pos1(1) + pos1(3), pos2(1) + pos2(3)) - rect_x;
    rect_h = max(pos1(2) + pos1(4), pos2(2) + pos2(4)) - rect_y;
    
    annotation('rectangle', [rect_x, rect_y, rect_w, rect_h], 'Color', 'k', 'LineWidth', 1.0, 'LineStyle', '--');
end

% Set axis limits for all the axes
axis(axesHandles, [-0.40 0 -0.061 0.061 -0.16 0.16]);




%%
function create_subplot(subU,xiPlot,nr_sv,nodes,elementPlot,L,O, textPos)
    v = reshape(subU*xiPlot, 3, []).';
    dm = PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
    plotcube(L,O,.05,[0 0 0]);
    subplotName = strcat('$$\xi_',nr_sv,'=',num2str(xiPlot),'$$');
%     text(textPos(1), textPos(2), textPos(3), subplotName,'Interpreter','latex')
end

