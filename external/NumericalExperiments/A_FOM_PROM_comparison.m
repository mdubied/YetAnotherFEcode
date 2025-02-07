% ------------------------------------------------------------------------ 
% A_FOM_PROM_comparison.m
%
% Description: Accuracy and computational speed comparison between the FOM
% and the (P)ROM formulations.
%
% Last modified: 06/02/2025, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
clear; 
close all; 
clc
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');

%% PREPARE MODELS _________________________________________________________                                                   

% load material and mesh parameters
load('parameters.mat') 

% parameters to test: propRigid, muscleBoundaries
% propRigid = 0.55;
muscleBoundaries = [0.83,0.6]; % 0.82 for the plot
%%
% specify and create FE mesh
filename = 'InputFiles/3d_rectangle_8086el';%_24822el'; %'3d_rectangle_8086el'
%'3d_rectangle_1272el';%'3d_rectangle_1272el';%'3d_rectangle_660el'; 4270,
%16009 3d_rectangle_47423el


[Mesh_ROM, ~, ~, ~, ~] = create_mesh(filename, myElementConstructor, propRigid);
[Mesh_FOM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);

% set boundary conditions
for l=1:length(nsetBC)
    Mesh_ROM.set_essential_boundary_condition([nsetBC{l}],1:3,0)  % all DOFs constrained to get VMs. Rigid body modes are added in build_ROM
    Mesh_FOM.set_essential_boundary_condition([nsetBC{l}],2:3,0)
end  

% shape variations for PROM
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

U_3 = [z_tail,z_head,y_thinFish];
U_5 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
U_8 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
        z_smallFish, z_notch, xz_concaveTail];
%% FIGURE HIGHLIGHTING THE SPINE AND TAIL ELEMENT _________________________
f_spine = create_fig_spine(elements, nodes, 'cyan', 'r', [1 2 16 8]);

%% FIGURE A1 (muscles, rigid part, VM) ____________________________________
f_A1 = create_fig_muscle_placement_VM(Mesh_ROM, nodes, elements,muscleBoundaries, esetBC);
fig_filename = sprintf('Setup/Figures/A_muscles_placement_VM_%del.pdf', Mesh_ROM.nElements);
exportgraphics(f_A1, fig_filename, 'Resolution', 1400);

%% SIMULATION PARAMETERS __________________________________________________
h = 0.01;
tmax = 0.3;
kActu = 6.5e5;  % [0.8,0.6] and 2.5e5 --> reaches 1.53 for 24822

%% FOM ____________________________________________________________________
tStartFOM = tic;

% build FOM
fprintf('____________________\n')
fprintf('Building FOM ... \n')
[Assembly,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
build_FOM_3D(Mesh_FOM,nodes,elements,muscleBoundaries);  

% %%
% % static solution
% Fext = tip_actuation_force_FOM_1(0.25,kActu,Assembly,tailProperties);
% [ u_lin, u ] = static_equilibrium(Assembly, Fext);
% 
% % y displacement of tail
% disp(u_lin(tailProperties.tailNode*3-1))
% disp(u(tailProperties.tailNode*3-1))   
% 
% %% Plot
% elementPlot = elements(:,1:4);
% u_plot = reshape(u, 3, []).';
% figure
% PlotFieldonDeformedMesh(nodes, elementPlot, u_plot, 'factor', 1);
% hold off



            
%%
% solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_FOM = solve_EoMs_FOM(Assembly,elements,tailProperties,spineProperties,dragProperties,actuLeft,actuRight,kActu,h,tmax); 
toc

fprintf('Time needed to solve the problem using FOM: %.2fsec\n',toc(tStartFOM))
timeFOM = toc(tStartFOM);
sol_FOM = TI_NL_FOM.Solution.u;
fprintf('Tail node lateral displacement, dynamic: %.4f\n', max(sol_FOM(tailProperties.tailNode*3-1,:))*100)

%% Plot
elementPlot = elements(:,1:4);
sol_FOM = TI_NL_FOM.Solution.u;
u_plot = reshape(sol_FOM(:,0.25/h), 3, []).';
figure
PlotFieldonDeformedMesh(nodes, elementPlot, u_plot, 'factor', 1);

% hold off

%%
sol_FOM = TI_NL_FOM.Solution.u;
fprintf('\n')
fprintf('______________________\n')
fprintf('Tail node lateral displacement, static : %.4f\n', u(tailProperties.tailNode*3-1)*100)
% fprintf('Tail node lateral displacement, dynamic: %.4f\n', max(sol_FOM(tailProperties.tailNode*3-1,:))*100)
fprintf('-------\n')

left_node = find_node(-0.2,0.02,0,nodes);
fprintf('Left node lateral displacement, static : %.4f\n', u(left_node*3-1)*100)
% fprintf('Left node lateral displacement, dynamic: %.4f\n', max(sol_FOM(left_node*3-1,:))*100)
fprintf('-------\n')

top_left_node = find_node(-0.2,0.02,0.05,nodes);
fprintf('Top left node lateral displacement, static : %.4f\n', u(top_left_node*3-1)*100)
% fprintf('Top left node lateral displacement, dynamic: %.4f\n', max(sol_FOM(top_left_node*3-1,:))*100)
fprintf('-------\n')

fprintf('Max lateral displacement, static : %.4f\n', max(u(2:3:end))*100)
% fprintf('Max lateral displacement, dynamic: %.4f\n', max(sol_FOM(2:3:end,:),[],'all')*100)
fprintf('-------\n')

%%
xdata = [1272,4270,8086,16009,24822];
ydata = [2.35,3.24,3.81,4.02,4.19];
ydata = [1.21,1.71,1.98,2.07,2.18];
ydata = [1.26,1.80,2.18,2.37,2.50];
figure
plot(xdata,ydata)
xlabel('Number of FEs')
ylabel('Lateral displacement (tail node)')
xticks([1000,5000,10000,15000,20000, 25000])
xticklabels({1000,5000,10000,15000,20000,25000})
grid on

%% ROM ____________________________________________________________________
tStartROM = tic;
kActu = kActu;
% build PROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')
[V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
build_ROM_3D(Mesh_ROM,nodes,elements,muscleBoundaries,USEJULIA);  
%%
% solve EoMs 
tic
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_ROM = solve_EoMs(V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight,kActu,h,tmax); 
toc

fprintf('Time needed to solve the problem using ROM: %.2fsec\n',toc(tStartROM))
timeROM = toc(tStartROM);

%% PROM ___________________________________________________________________
tStartPROM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building PROM ... \n')
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
build_PROM_3D(Mesh_ROM,nodes,elements,muscleBoundaries,U,USEJULIA,VOLUME,FORMULATION);      

% solve EoMs (with sensitivities for the PROM)
tic 
fprintf('____________________\n')
fprintf('Solving EoMs and sensitivities...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight,kActu,h,tmax); 
toc

fprintf('Time needed to solve the problem using PROM: %.2fsec\n',toc(tStartPROM))
timePROM = toc(tStartPROM);


%% PLOT ___________________________________________________________________
uTail_FOM = zeros(3,tmax/h);
uTail_ROM = zeros(3,tmax/h);
uTail_PROM = zeros(3,tmax/h);
uHead_FOM = zeros(3,tmax/h);
uHead_ROM = zeros(3,tmax/h);
uHead_PROM = zeros(3,tmax/h);

modelToPlot = ['FOM'];%,'ROM'];

if contains(modelToPlot,'FOM')
    sol_FOM = TI_NL_FOM.Solution.u;
end
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));

headNode = find_node(0,0,0,nodes);

for t=1:tmax/h
     if contains(modelToPlot,'FOM')
        uTail_FOM(:,t) = sol_FOM(tailProperties.tailNode*3-2:tailProperties.tailNode*3,t);
        uHead_FOM(:,t) = sol_FOM(headNode*3-2:headNode*3,t);
     end
     if contains(modelToPlot,'ROM')
        uTail_ROM(:,t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_ROM.Solution.q(:,t);  
        uHead_ROM(:,t) = V(headNode*3-2:headNode*3,:)*TI_NL_ROM.Solution.q(:,t); 
     end
end

f_A2 = figure('units','centimeters','position',[3 3 9 6]);
% x-position
subplot(2,1,1);
hold on
if contains(modelToPlot,'FOM')
    plot(timePlot,x0Tail+uHead_FOM(1,:),'--','DisplayName','FOM')
end
if contains(modelToPlot,'ROM')
    plot(timePlot,x0Tail+uHead_ROM(1,:),'DisplayName','ROM')
end
% plot(timePlot,x0Tail+uTail_PROM(1,:),'DisplayName','PROM')
hold on
grid on
xlabel('Time [s]')
ylabel('head x-position [m]')
legend('Location','northwest', 'interpreter','latex')

% y-position
subplot(2,1,2);
hold on
if contains(modelToPlot,'FOM')
    plot(timePlot,uTail_FOM(2,:),'--','DisplayName','FOM')
end
if contains(modelToPlot,'ROM')
    plot(timePlot,uTail_ROM(2,:),'DisplayName','ROM')
end

% plot(timePlot,uTail_PROM(2,:),'DisplayName','PROM')
grid on
% ylim([-0.03,0.03])
xlabel('Time [s]')
ylabel('tail y-position [m]')
% legend('Location','southwest', 'interpreter','latex')
% exportgraphics(f_A2,'A_FOM_vs_ROM_1272el.pdf','Resolution',1400)
% exportgraphics(f_A2,strcat('A_FOM_vs_ROM_',num2str(size(elements,1)),'.jpg'),'Resolution',600)

%% COMPUTE AVERAGE ERROR
error_x_vector = zeros(1,tmax/h);
error_y_vector = zeros(1,tmax/h);
for t=2:tmax/h
    error_x_vector(t) = (uHead_FOM(1,t) - uHead_ROM(1,t))/uHead_FOM(1,t);
    error_y_vector(t) = (uTail_FOM(2,t) - uTail_ROM(2,t))/uTail_FOM(2,t);
end
% mean(error_x_vector)
% mean(error_y_vector)
error_x_vector(end)

%%
nElVec = [660,1272,4270,8086,24822];
timeROM = [41.09,22.35,25.39,29.06,50.30];
timeFOM = [49.0,83.12,327.68,601.84,2421];

f_comp_time = figure('units','centimeters','position',[3 3 9 6]);
plot(nElVec,(timeROM),'-x')
hold on
plot(nElVec,(timeFOM),'-*')
grid on
xlabel('Number of FE')
ylabel('Time [s]')
legend('ROM', 'FOM','Location','northwest', 'interpreter','latex')
%% ________________________________________________________________________
%
%
%
%% ITERATIVE TESTING OF MULTIPLE CASES ____________________________________

% Vector of element counts
elements_vec = [4270]; % Number of elements for each input file
kActu_values = [3.0e5];    % Actuation values 1.0e5,2.0e5,3.0e5,4.0e5

% 4270el:   5.0e5 limit (doesn't make it)
% 8086el:   4.0e5 just make it
% 16009el:  3.5e5 just ok, 4.0 for 
% 16009el:  4.0e5 just ok
% 24822el:  3.0e5 not ok, 2.75 ok


% Set simulation parameters
h = 0.02;
tmax = 2.0;

% Pre-allocate matrix to store results. 16 columns
% Columns:  [num_elements, kActu, ...
%             max_uTail_FOM, max_uTail_ROM, abs_diff_Tail, rel_error_Tail, ...
%             max_uHead_FOM, max_uHead_ROM, abs_diff_Head, rel_error_Head, ...
%             rel_error_Tail_ROM_PROM_3, rel_error_Head_ROM_PROM_3, ...
%             rel_error_Tail_ROM_PROM_5, rel_error_Head_ROM_PROM_5, ...
%             rel_error_Tail_ROM_PROM_8, rel_error_Head_ROM_PROM_8, ...
%             timeFOMBuild, timeROMBuild, timeFOMSolve, timeROMSolve, timeFOM, timeROM, ...
%             timePROMBuild_3, timePROMSolve_3, timePROM_3,...
%             timePROMBuild_5, timePROMSolve_5, timePROM_5,
%             timePROMBuild_8, timePROMSolve_8, timePROM_8];
results_matrix = zeros(length(elements_vec) * numel(kActu_values), 31); 

% Initialize row index for results_matrix
row_idx = 1;

% Loop over each element count and each actuation value
for elem_idx = 1:length(elements_vec)
    num_elements = elements_vec(elem_idx); % Get current number of elements
    
    % Load the FE mesh for FOM and ROM
    filename = strcat('InputFiles/3d_rectangle_', num2str(num_elements), 'el');  % Construct filena
    [Mesh_ROM, ~, ~, ~, ~] = create_mesh(filename, myElementConstructor, propRigid);
    [Mesh_FOM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
    [Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
    
    % shape variations for PROM
    [y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
        y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail] = ...
            shape_variations_3D(nodes,Lx,Ly,Lz);
    
    U_3 = [z_tail,z_head,y_thinFish];
    U_5 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
    U_8 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
            z_smallFish, z_notch, xz_concaveTail];

    % Set boundary conditions for FOM and ROM
    for l = 1:length(nsetBC)
        Mesh_ROM.set_essential_boundary_condition([nsetBC{l}], 1:3, 0);
        Mesh_FOM.set_essential_boundary_condition([nsetBC{l}], 2:3, 0);
    end
    
    for k_idx = 1:size(kActu_values,2)
        kActu = kActu_values(elem_idx,k_idx);  % Get current actuation value
    
        % FOM Simulation __________________________________
        tStartFOM = tic;
        fprintf('Building and solving FOM for %d elements, kActu: %.3f...\n', num_elements, kActu);
        [Assembly, tailProperties, spineProperties, dragProperties, actuLeft, actuRight] = ...
            build_FOM_3D(Mesh_FOM, nodes, elements, muscleBoundaries);
        timeFOMBuild = toc(tStartFOM);
        tStartFOMSolve = tic;
        TI_NL_FOM = solve_EoMs_FOM(Assembly, elements, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax);
        timeFOMSolve = toc(tStartFOMSolve);
        timeFOM = toc(tStartFOM);
           
        % ROM Simulation __________________________________
        tStartROM = tic;
        fprintf('Building and solving ROM for %d elements, kActu: %.3f...\n', num_elements, kActu);
        [V, ROM_Assembly, tensors_ROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight] = ...
            build_ROM_3D(Mesh_ROM, nodes, elements, muscleBoundaries, USEJULIA);
        timeROMBuild = toc(tStartROM);
        tStartROMSolve = tic;
        TI_NL_ROM = solve_EoMs(V, ROM_Assembly, tensors_ROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax);
        timeROMSolve = toc(tStartROMSolve);
        timeROM = toc(tStartROM);
        
        % PROM Simulation, 3 parameters ___________________
        tStartPROM_3 = tic;
        fprintf('Building and solving PROM (3 param., w/ sens.) for %d elements, kActu: %.3f...\n', num_elements, kActu)      
        [V_3,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
            build_PROM_3D(Mesh_ROM, nodes, elements, muscleBoundaries, U_3, USEJULIA, VOLUME, FORMULATION); 
        
        timePROMBuild_3 = toc(tStartPROM_3);
        tStartPROMSolve_3 = tic;
        TI_NL_PROM_3 = solve_EoMs_and_sensitivities(V_3, PROM_Assembly, tensors_PROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax); 
        timePROMSolve_3 = toc(tStartPROMSolve_3);
        timePROM_3 = toc(tStartPROM_3);
        
        % PROM Simulation, 5 parameters ___________________
        tStartPROM_5 = tic;
        fprintf('Building and solving PROM (5 param., w/ sens.) for %d elements, kActu: %.3f...\n', num_elements, kActu)      
        [V_5,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
            build_PROM_3D(Mesh_ROM, nodes, elements, muscleBoundaries, U_5, USEJULIA, VOLUME, FORMULATION); 
        
        timePROMBuild_5 = toc(tStartPROM_5);
        tStartPROMSolve_5 = tic;
        TI_NL_PROM_5 = solve_EoMs_and_sensitivities(V_5, PROM_Assembly, tensors_PROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax); 
        timePROMSolve_5 = toc(tStartPROMSolve_5);
        timePROM_5 = toc(tStartPROM_5);
        
        % PROM Simulation, 8 parameters ___________________
        tStartPROM_8 = tic;
        fprintf('Building and solving PROM (8 param., w/ sens.) for %d elements, kActu: %.3f...\n', num_elements, kActu)      
        [V_8,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
            build_PROM_3D(Mesh_ROM, nodes, elements, muscleBoundaries, U_8, USEJULIA, VOLUME, FORMULATION); 
        
        timePROMBuild_8 = toc(tStartPROM_8);
        tStartPROMSolve_8 = tic;
        TI_NL_PROM_8 = solve_EoMs_and_sensitivities(V_8, PROM_Assembly, tensors_PROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax); 
        timePROMSolve_8 = toc(tStartPROMSolve_8);
        timePROM_8 = toc(tStartPROM_8);
        
        
        % Data Analysis ___________________________________
        fprintf('Post-processing results for %d elements, kActu: %.3f...\n\n\n', num_elements, kActu)
        sol_FOM = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);
        timePlot = linspace(0, tmax-h, tmax/h);
        x0Tail = min(nodes(:, 1));
        
        headNode = find_node(0, 0, 0, nodes);
        
        uTail_FOM = zeros(3, tmax/h); 
        uTail_ROM = zeros(3, tmax/h);
        uTail_PROM_3 = zeros(3, tmax/h); 
        uTail_PROM_5 = zeros(3, tmax/h);
        uTail_PROM_8 = zeros(3, tmax/h);
        
        uHead_FOM = zeros(3, tmax/h);
        uHead_ROM = zeros(3, tmax/h);
        uHead_PROM_3 = zeros(3, tmax/h);
        uHead_PROM_5 = zeros(3, tmax/h);
        uHead_PROM_8 = zeros(3, tmax/h);
        
        for t = 1:tmax/h
            % FOM displacement, convert to cm
            uTail_FOM(:, t) = sol_FOM(tailProperties.tailNode*3-2:tailProperties.tailNode*3, t)*100;
            uHead_FOM(:, t) = sol_FOM(headNode*3-2:headNode*3, t)*100;
            
            % ROM displacement, convert to cm
            uTail_ROM(:, t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3, :) * TI_NL_ROM.Solution.q(:, t)*100;
            uHead_ROM(:, t) = V(headNode*3-2:headNode*3, :) * TI_NL_ROM.Solution.q(:, t)*100;
            
            % PROM (3 params) displacement, convert to cm
            uTail_PROM_3(:, t) = V_3(tailProperties.tailNode*3-2:tailProperties.tailNode*3, :) * TI_NL_PROM_3.Solution.q(:, t)*100;
            uHead_PROM_3(:, t) = V_3(headNode*3-2:headNode*3, :) * TI_NL_PROM_3.Solution.q(:, t)*100;
            
            % PROM (5 params) displacement, convert to cm
            uTail_PROM_5(:, t) = V_5(tailProperties.tailNode*3-2:tailProperties.tailNode*3, :) * TI_NL_PROM_5.Solution.q(:, t)*100;
            uHead_PROM_5(:, t) = V_5(headNode*3-2:headNode*3, :) * TI_NL_PROM_5.Solution.q(:, t)*100;
            
            % PROM (8 params) displacement, convert to cm
            uTail_PROM_8(:, t) = V_8(tailProperties.tailNode*3-2:tailProperties.tailNode*3, :) * TI_NL_PROM_8.Solution.q(:, t)*100;
            uHead_PROM_8(:, t) = V_8(headNode*3-2:headNode*3, :) * TI_NL_PROM_8.Solution.q(:, t)*100;
        end
        
        % Calculate maximum displacements (lateral, 2, and horizontal, 1)
        max_uTail_FOM = max(uTail_FOM(2, :));   
        max_uHead_FOM = max(uHead_FOM(1, :));   
        max_uTail_ROM = max(uTail_ROM(2, :));   
        max_uHead_ROM = max(uHead_ROM(1, :));  
        max_uTail_PROM_3 = max(uTail_PROM_3(2, :));   
        max_uHead_PROM_3 = max(uHead_PROM_3(1, :));
        max_uTail_PROM_5 = max(uTail_PROM_5(2, :));   
        max_uHead_PROM_5 = max(uHead_PROM_5(1, :));
        max_uTail_PROM_8 = max(uTail_PROM_8(2, :));   
        max_uHead_PROM_8 = max(uHead_PROM_8(1, :));
        
        % Calculate absolute differences and relative errors
        abs_diff_Tail = max_uTail_ROM - max_uTail_FOM;
        abs_diff_Head = max_uHead_ROM - max_uHead_FOM;
        rel_error_Tail = abs_diff_Tail / max_uTail_FOM;
        rel_error_Head = abs_diff_Head / max_uHead_FOM;
        
        abs_diff_Tail_ROM_PROM_3 = max_uTail_PROM_3 - max_uTail_ROM;
        abs_diff_Head_ROM_PROM_3 = max_uHead_PROM_3 - max_uHead_ROM;
        abs_diff_Tail_ROM_PROM_5 = max_uTail_PROM_5 - max_uTail_ROM;
        abs_diff_Head_ROM_PROM_5 = max_uHead_PROM_5 - max_uHead_ROM;
        abs_diff_Tail_ROM_PROM_8 = max_uTail_PROM_8 - max_uTail_ROM;
        abs_diff_Head_ROM_PROM_8 = max_uHead_PROM_8 - max_uHead_ROM;
        rel_error_Tail_ROM_PROM_3 = abs_diff_Tail_ROM_PROM_3 / max_uTail_ROM;
        rel_error_Head_ROM_PROM_3 = abs_diff_Head_ROM_PROM_3 / max_uHead_ROM;
        rel_error_Tail_ROM_PROM_5 = abs_diff_Tail_ROM_PROM_5 / max_uTail_ROM;
        rel_error_Head_ROM_PROM_5 = abs_diff_Head_ROM_PROM_5 / max_uHead_ROM;
        rel_error_Tail_ROM_PROM_8 = abs_diff_Tail_ROM_PROM_8 / max_uTail_ROM;
        rel_error_Head_ROM_PROM_8 = abs_diff_Head_ROM_PROM_8 / max_uHead_ROM;
        
        % Store results in the results_matrix
        results_matrix(row_idx, :) = [num_elements, kActu, ...
            max_uTail_FOM, max_uTail_ROM, abs_diff_Tail, rel_error_Tail, ...
            max_uHead_FOM, max_uHead_ROM, abs_diff_Head, rel_error_Head, ...
            rel_error_Tail_ROM_PROM_3, rel_error_Head_ROM_PROM_3, ...
            rel_error_Tail_ROM_PROM_5, rel_error_Head_ROM_PROM_5, ...
            rel_error_Tail_ROM_PROM_8, rel_error_Head_ROM_PROM_8, ...
            timeFOMBuild, timeROMBuild, timeFOMSolve, timeROMSolve, timeFOM, timeROM, ...
            timePROMBuild_3, timePROMSolve_3, timePROM_3, ...
            timePROMBuild_5, timePROMSolve_5, timePROM_5, ...
            timePROMBuild_8, timePROMSolve_8, timePROM_8];
        row_idx = row_idx + 1;
        
        % Save results to .mat file
        mat_filename = sprintf('Results/Data/temp_A_results_%del_kActu_%.3f.mat', num_elements, kActu);
        save(mat_filename, 'timeFOM', 'timeROM', ...
            'timeFOMBuild', 'timeROMBuild','timeFOMSolve', 'timeROMSolve', ...
            'uTail_FOM', 'uTail_ROM', 'uHead_FOM', 'uHead_ROM', ...
            'max_uTail_FOM', 'max_uHead_FOM', 'max_uTail_ROM', 'max_uHead_ROM', ...
            'abs_diff_Tail', 'abs_diff_Head', 'rel_error_Tail', 'rel_error_Head');
                
        % Create and save figure comparing FOM and ROM displacements
        f = figure('units','centimeters','position',[3 3 9 6.5]);

        % x-position (Head)
        subplot(2,1,1);
        hold on;
        plot(timePlot, uHead_FOM(1,:), '--', 'DisplayName', 'FOM');
        plot(timePlot, uHead_ROM(1,:), 'DisplayName', 'ROM');
        grid on;
        xlabel('Time [s]');
        ylTop = ylabel('Head x-position [cm]');
        legend('Location','northwest', 'interpreter', 'latex');
    
        % y-position (Tail)
        subplot(2,1,2);
        hold on;
        plot(timePlot, uTail_FOM(2,:), '--', 'DisplayName', 'FOM');
        plot(timePlot, uTail_ROM(2,:), 'DisplayName', 'ROM');
        grid on;
        xlabel('Time [s]');
        ylBottom = ylabel('Tail y-position [cm]');
        
        % Align y-axis labels
        if ylTop<ylBottom
            set(ylBottom,'Pos',[ylTop.Position(1) ylBottom.Position(2) ylBottom.Position(3)]);
        else
            set(ylTop,'Pos',[ylBottom.Position(1) ylTop.Position(2) ylTop.Position(3)]);
        end

        % Save figure
        fig_filename = sprintf('Results/Figures/PDF/temp_A_Displacement_Comparison_%del_kActu_%.3f.pdf', num_elements, kActu);
        exportgraphics(f, fig_filename, 'Resolution', 600);
        fig_filename = sprintf('Results/Figures/JPG/temp_A_Displacement_Comparison_%del_kActu_%.3f.jpg', num_elements, kActu);
        exportgraphics(f, fig_filename, 'Resolution', 600);

        % Close the figure to avoid memory issues during iterations
        close(f);
        
    end
    
    % Save summary table for each mesh
    results_table_filename = sprintf('Results/Data/temp_A_results_new_%del.csv', num_elements);
    csvwrite(results_table_filename,results_matrix);
end
%%
%% ________________________________________________________________________
%
%
%
%% RESULTS ANALYSIS AND PLOTS _____________________________________________

%% Read results from csv
nElementsForResult = [1272, 4272, 8086, 16009, 24822]; % Number of elements for each input file
% nElementsForResult = [1272, 4270];
for i=1:length(nElementsForResult)
    nElements = nElementsForResult(i);
    filename = sprintf('Results/Data/A_results_%del.csv', nElements); % sprintf('Results/Data/temp_A_results_new_%del.csv', nElements);
    if i==1
        resultSummary = readmatrix(filename);
    else
        resultSummary = [resultSummary;readmatrix(filename)];
    end
    
end

% Extract the unique cases from the first column
cases = unique(resultSummary(:, 1));

% Create an average matrix
avgResults = zeros(length(cases), size(resultSummary, 2));

% Loop through each case and compute the averages
for i = 1:length(cases)
    % Get rows corresponding to the current case
    rows = resultSummary(resultSummary(:, 1) == cases(i), :);
    
    % Compute the average for each column
    avgResults(i, :) = mean(rows, 1);
end

%% Create and save figure comparing FOM and ROM computational time
% Columns: [num_elements, kActu, ...
%  3-6      max_uTail_FOM, max_uTail_ROM, abs_diff_Tail, rel_error_Tail, ...
%  7-10     max_uHead_FOM, max_uHead_ROM, abs_diff_Head, rel_error_Head, ...
%  11-12    rel_error_Tail_ROM_PROM_3, rel_error_Head_ROM_PROM_3, ...
%  13-14    rel_error_Tail_ROM_PROM_5, rel_error_Head_ROM_PROM_5, ...
%  15-16    rel_error_Tail_ROM_PROM_8, rel_error_Head_ROM_PROM_8, ...
%  17-22    timeFOMBuild, timeROMBuild, timeFOMSolve, timeROMSolve, timeFOM, timeROM, ...
%  23-25    timePROMBuild_3, timePROMSolve_3, timePROM_3, ...
%  26-28    timePROMBuild_5, timePROMSolve_5, timePROM_5, ...
%  29-31    timePROMBuild_8, timePROMSolve_8, timePROM_8];
f = figure('units','centimeters','position',[3 3 9 15.0]);

subplot(3,1,1);
hold on;
p1=plot(avgResults(:,1), avgResults(:,21), '--', 'DisplayName', 'FOM', 'LineWidth', 1.0);
p2=plot(avgResults(:,1), avgResults(:,22), 'LineWidth', 1.0);
p3=plot(avgResults(:,1), avgResults(:,25), '-.', 'DisplayName', 'PROM (3 param.)', 'LineWidth', 1.0);
p4=plot(avgResults(:,1), avgResults(:,28), '-.', 'DisplayName', 'PROM (5 param.)', 'LineWidth', 1.0);
p5=plot(avgResults(:,1), avgResults(:,31), '-.', 'DisplayName', 'PROM (8 param.)', 'LineWidth', 1.0);
p6=plot(avgResults(:,1), avgResults(:,31),'LineStyle' ,"none");   % dummy lines for legend alignment
grid on;
xlim([1272 24822])
xlabel('Number of finite elements');
set(gca, 'YScale', 'log')
yticks([100,1000]);
ylTop = ylabel('Time [s]');
lgd = legend([p1,p2,p3,p4,p5],{'FOM', 'ROM', 'PROM 3p', 'PROM 5p', 'PROM 8p'},'Location','northwest', 'NumColumns',1, 'interpreter', 'latex');
lgd.Position(2) = lgd.Position(2) + 0.075; 

% Plot for computational time ROM only
subplot(3,1,2);

hold on;
fill([avgResults(:,1); flipud(avgResults(:,1))], [zeros(size(avgResults(:,1))); flipud(avgResults(:,18))],[0.8500 0.3250 0.0980], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Build ROM'); % Fill for build
fill([avgResults(:,1); flipud(avgResults(:,1))], [avgResults(:,18); flipud(avgResults(:,22))],[0.8500 0.3250 0.0980], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Simulate ROM'); % Fill for simulate

grid on
xlim([1272 24822])
xlabel('Number of finite elements');
ylMiddle = ylabel('Time [s]');
legend('Location','northwest', 'interpreter', 'latex');

% Plot for computational time ROM only
subplot(3,1,3);
hold on;
fill([avgResults(:,1); flipud(avgResults(:,1))], [zeros(size(avgResults(:,1))); flipud(avgResults(:,26))],[0.4940 0.1840 0.5560], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Build PROM 5p'); % Fill for build
fill([avgResults(:,1); flipud(avgResults(:,1))], [avgResults(:,26); flipud(avgResults(:,28))],[0.4940 0.1840 0.5560], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Simulate  PROM 5p'); % Fill for simulate
grid on
xlim([1272 24822])
xlabel('Number of finite elements');
ylBottom = ylabel('Time [s]');
legend('Location','northwest', 'interpreter', 'latex');



% Align y-axis labels
if ylTop<ylBottom
    set(ylBottom,'Pos',[ylTop.Position(1) ylBottom.Position(2) ylBottom.Position(3)]);
    set(ylMiddle,'Pos',[ylTop.Position(1) ylMiddle.Position(2) ylMiddle.Position(3)]);
else
    set(ylTop,'Pos',[ylBottom.Position(1) ylTop.Position(2) ylTop.Position(3)]);
    set(ylMiddle,'Pos',[ylBottom.Position(1) ylMiddle.Position(2) ylMiddle.Position(3)]);
end

fig_filename = 'Results/Figures/PDF/08_A_results_computational_time.pdf';
exportgraphics(f, fig_filename, 'Resolution', 600);

%% Plot accuracy as heat map
% Columns: [num_elements, kActu, ...
%  3-6      max_uTail_FOM, max_uTail_ROM, abs_diff_Tail, rel_error_Tail, ...
%  7-10     max_uHead_FOM, max_uHead_ROM, abs_diff_Head, rel_error_Head, ...
%  11-12    rel_error_Tail_ROM_PROM_3, rel_error_Head_ROM_PROM_3, ...
%  13-14    rel_error_Tail_ROM_PROM_5, rel_error_Head_ROM_PROM_5, ...
%  15-16    rel_error_Tail_ROM_PROM_8, rel_error_Head_ROM_PROM_8, ...
%  17-22    timeFOMBuild, timeROMBuild, timeFOMSolve, timeROMSolve, timeFOM, timeROM, ...
%  23-25    timePROMBuild_3, timePROMSolve_3, timePROM_3, ...
%  26-28    timePROMBuild_5, timePROMSolve_5, timePROM_5, ...
%  29-31    timePROMBuild_8, timePROMSolve_8, timePROM_8];

f = figure('units','centimeters','position',[3 3 9 7]);
hold on
X = resultSummary(:, 1);            % nr. elements
Y = resultSummary(:, 3);            % max horizontal displacement
values = resultSummary(:, 10)*100;  % relative error 

% Number of grid points
gridres = 500 ;
% Create a uniform vector on X, from min to max value, made of "gridres" points
xs = linspace(min(X),max(X),gridres) ;
% Create a uniform vector on Y, from min to max value, made of "gridres" points
ys = linspace(min(Y),max(Y),gridres) ;
% Generate 2D grid coordinates from xs and ys
[xq,yq]=meshgrid(xs,ys) ;

% Interpolate the values over the new grid
InterpolatedValues = griddata(X,Y,values,xq,yq,'linear') ;

hmap_above = pcolor(xq,yq,InterpolatedValues);
hmap_above.EdgeColor = [.5 .5 .5] ; % cosmetic adjustment
cb = colorbar('TickLabelInterpreter','latex');
colormap jet
shading interp

% Draw zero line
contour(xq,yq,InterpolatedValues,[0, 0],'LineColor','k','LineWidth',1)

% Draw zero line on colorbar
cbar_limits = caxis();
zero_pos = (0 - cbar_limits(1)) / (cbar_limits(2) - cbar_limits(1));
cbar_pos = cb.Position;
y_zero = cbar_pos(2) + zero_pos * cbar_pos(4);
annotation('line', [cbar_pos(1), cbar_pos(1) + cbar_pos(3)], [y_zero, y_zero], 'Color', 'k', 'LineWidth', 1.0);



% Labels and appearancy
xlabel('Number of finite elements')
ylabel('Tail $y$-oscillation [cm]')
xlim([1700,24000])
ylim([0.5,2.0])
title('Relative error in $x$-position after $t=2$s', 'interpreter', 'latex')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

cb.Label.String = '[$x$(ROM)  - $x$(FOM)] / $x$(FOM) [\%]';
cb.Label.Interpreter = 'latex'; % If you want LaTeX formatting





hold off



fig_filename = 'Results/Figures/PDF/09_A_results_rel_error_forward.pdf';
exportgraphics(f, fig_filename, 'Resolution', 600);























%% ANIMATION ______________________________________________________________
elementPlot = elements(:,1:4); 
nel = size(elements,1);

% top muscle
topMuscle = zeros(nel,1);

for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*propRigid && elementCenterX > -Lx*0.9
        topMuscle(el) = 1;
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*propRigid && elementCenterX > -Lx*0.9
        bottomMuscle(el) = 1;
    end    

end

actuationValues = zeros(size(TI_NL_FOM.Solution.u,2),1);
for t=1:size(TI_NL_FOM.Solution.u,2)
    actuationValues(t) = 0;
end

actuationValues2 = zeros(size(TI_NL_FOM.Solution.u,2),1);
for t=1:size(TI_NL_FOM.Solution.u,2)
    actuationValues2(t) = 0;
end
sol = TI_NL_FOM.Solution.u(:,1:end);
AnimateFieldonDeformedMeshActuation2Muscles(nodes, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,sol, ...
    'factor',1,'index',1:3,'filename','result_video','framerate',1/h)






