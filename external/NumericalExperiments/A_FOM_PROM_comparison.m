% ------------------------------------------------------------------------ 
% A_FOM_PROM_comparison.m
%
% Description: Accuracy and computational speed comparison between the FOM
% and the (P)ROM formulations.
%
% Last modified: 13/01/2025, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
clear; 
close all; 
clc
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');

%% PREPARE MODELS _________________________________________________________                                                   

% load material and mesh parameters
load('parameters.mat') 


% parameters to test: propRigid, muscleBoundaries
% propRigid = 0.55;
%%
% specify and create FE mesh
filename = 'InputFiles/3d_rectangle_1272el';%_24822el'; %'3d_rectangle_8086el'
%'3d_rectangle_1272el';%'3d_rectangle_1272el';%'3d_rectangle_660el'; 4270


[Mesh_ROM, ~, ~, ~, ~] = create_mesh(filename, myElementConstructor, propRigid);
[Mesh_FOM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);

% plot mesh
% elementPlot = elements(:,1:4);
% figure('units','normalized','position',[.2 .1 .6 .8])
% PlotMeshAxis(nodes, elementPlot, 0);
% hold off

% set boundary conditions
for l=1:length(nsetBC)
    Mesh_ROM.set_essential_boundary_condition([nsetBC{l}],1:3,0)  % all DOFs constrained to get VMs. Rigid body modes are added in build_ROM
    Mesh_FOM.set_essential_boundary_condition([nsetBC{l}],2:3,0)
end  

% shape variations for PROM
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);
U = [z_tail,y_head,y_thinFish,z_head,z_linLongTail];

%% FIGURE HIGHLIGHTING THE SPINE AND TAIL ELEMENT _________________________
f_spine = create_fig_spine(elements, nodes, 'cyan', 'r', [1 2 16 8]);

%% FIGURE A1 (muscles, rigid part, VM) ____________________________________
% Note: the position of the muscle is defined in the build_ROM/FOM/PROM
% functions. The number of VMs used also. Only the constraints for the
% rigid par of the fish is defined above (boundary conditions)
muscleBoundaries = [0.9,0.60];
f_A1 = create_fig_muscle_placement_VM(Mesh_ROM, nodes, elements,muscleBoundaries, esetBC);
% exportgraphics(f_A1,'A_muscles_placement_VM.pdf','Resolution',1400)

%% SIMULATION PARAMETERS __________________________________________________
h = 0.01;
tmax = 1.0;
kActu = 5e5;

%% FOM ____________________________________________________________________
tStartFOM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building FOM ... \n')
[Assembly,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
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

%% TEST ACTUATION FORCE
% quiver3(X,Y,Z,U,V,W) plots arrows with directional components U, V, and W at the Cartesian coordinates specified by X, Y, and Z
% actuation force
B1T = actuTop.B1;
B1B = actuBottom.B1;
B2T = actuTop.B2;
B2B = actuBottom.B2; 


fActu = @(t,q)  kActu/2*(-0.2*sin(t*2*pi))*(B1T+B2T*q) + ...
                kActu/2*(0.2*sin(t*2*pi))*(B1B+B2B*q);
            
nUncDOFs = size(Assembly.Mesh.EBC.unconstrainedDOFs,2);
nDOFs = Assembly.Mesh.nDOFs;
q0 = zeros(nUncDOFs,1);
u0 = Assembly.unconstrain_vector(q0);
            
fActuTest = fActu(0.2,u0);
fActuPlot = reshape(fActuTest, 3, []).';
figure
variation = reshape(u0*1, 3, []).';
% PlotFieldonDeformedMesh(nodes, elements, variation, 'factor', 1);
quiver3(nodes(:,1),nodes(:,2),nodes(:,3),fActuPlot(:,1),fActuPlot(:,2),fActuPlot(:,3),5)
axis equal

            
%%
% solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_FOM = solve_EoMs_FOM(Assembly,elements,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax); 
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
[V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_ROM_3D(Mesh_ROM,nodes,elements,muscleBoundaries,USEJULIA);  

% solve EoMs 
tic
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_ROM = solve_EoMs(V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax); 
toc

fprintf('Time needed to solve the problem using ROM: %.2fsec\n',toc(tStartROM))
timeROM = toc(tStartROM);

%% PROM ___________________________________________________________________
tStartPROM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building PROM ... \n')
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_PROM_3D(Mesh_ROM,nodes,elements,muscleBoundaries,U,USEJULIA,VOLUME,FORMULATION);      

% solve EoMs (with sensitivities for the PROM)
tic 
fprintf('____________________\n')
fprintf('Solving EoMs and sensitivities...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax); 
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

modelToPlot = ['ROM'];%,'ROM'];

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


%% ITERATIVE TESTING OF MULTIPLE CASES ____________________________________

% Vector of element counts
elements_vec = [1272]; % Number of elements for each input file
kActu_values = [0.4,0.8,1.3,1.6,2.0];    % Actuation values
            
            % 1.0, 2.0, 3.0, 4.0;

% Set simulation parameters
h = 0.02;
tmax = 2.0;

% Pre-allocate matrix to store results
% Columns: [num_elements, kActu, max_uTail_FOM, relative_error_Tail, relative_error_Head,
% timeFOMBuild, timeROMBuild, timeFOMSolve, timeROMSolve, timeFOM, timeROM]
results_matrix = zeros(length(elements_vec) * numel(kActu_values), 11); 

% Initialize row index for results_matrix
row_idx = 1;

% Loop over each element count and each actuation value
for elem_idx = 1:length(elements_vec)
    num_elements = elements_vec(elem_idx); % Get current number of elements
    filename = strcat('3d_rectangle_', num2str(num_elements), 'el');  % Construct filename
    
    for k_idx = 1:size(kActu_values,2)
        kActu = kActu_values(elem_idx,k_idx);  % Get current actuation value
        
        % Simulation of the models
        % Load the FE mesh for FOM and ROM
        [Mesh_ROM, ~, ~, ~, ~] = create_mesh(filename, myElementConstructor, propRigid);
        [Mesh_FOM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
        [Lx, Ly, Lz] = mesh_dimensions(nodes);
        
        % Set boundary conditions for FOM and ROM
        for l = 1:length(nsetBC)
            Mesh_ROM.set_essential_boundary_condition([nsetBC{l}], 1:3, 0);
            Mesh_FOM.set_essential_boundary_condition([nsetBC{l}], 2:3, 0);
        end
        
        % Shape variations for PROM
        [y_thinFish, z_smallFish, z_tail, z_head, z_linLongTail, z_notch, ...
            y_tail, y_head, y_linLongTail, y_ellipseFish] = shape_variations_3D(nodes, Lx, Ly, Lz);
        U = [z_tail, y_head, y_thinFish, z_head, z_linLongTail];
        
        % FOM Simulation
        tStartFOM = tic;
        fprintf('Building and solving FOM for %d elements, kActu: %.3f...\n', num_elements, kActu);
        [Assembly, tailProperties, spineProperties, dragProperties, actuTop, actuBottom] = ...
            build_FOM_3D(Mesh_FOM, nodes, elements);
        timeFOMBuild = toc(tStartFOM);
        tStartFOMSolve = tic;
        TI_NL_FOM = solve_EoMs_FOM(Assembly, elements, tailProperties, spineProperties, dragProperties, actuTop, actuBottom, kActu, h, tmax);
        timeFOMSolve = toc(tStartFOMSolve);
        timeFOM = toc(tStartFOM);
        
        
        % ROM Simulation
        tStartROM = tic;
        fprintf('Building and solving ROM for %d elements, kActu: %.3f...\n', num_elements, kActu);
        [V, ROM_Assembly, tensors_ROM, tailProperties, spineProperties, dragProperties, actuTop, actuBottom] = ...
            build_ROM_3D(Mesh_ROM, nodes, elements, USEJULIA);
        timeROMBuild = toc(tStartROM);
        tStartROMSolve = tic;
        TI_NL_ROM = solve_EoMs(V, ROM_Assembly, tensors_ROM, tailProperties, spineProperties, dragProperties, actuTop, actuBottom, kActu, h, tmax);
        timeROMSolve = toc(tStartROMSolve);
        timeROM = toc(tStartROM);
        
        % Data Analysis
        sol_FOM = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);
        timePlot = linspace(0, tmax-h, tmax/h);
        x0Tail = min(nodes(:, 1));
        
        headNode = find_node(0, 0, 0, nodes);
        
        uTail_FOM = zeros(3, tmax/h);
        uTail_ROM = zeros(3, tmax/h);
        uHead_FOM = zeros(3, tmax/h);
        uHead_ROM = zeros(3, tmax/h);
        
        for t = 1:tmax/h
            % FOM displacement
            uTail_FOM(:, t) = sol_FOM(tailProperties.tailNode*3-2:tailProperties.tailNode*3, t);
            uHead_FOM(:, t) = sol_FOM(headNode*3-2:headNode*3, t);
            
            % ROM displacement
            uTail_ROM(:, t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3, :) * TI_NL_ROM.Solution.q(:, t);
            uHead_ROM(:, t) = V(headNode*3-2:headNode*3, :) * TI_NL_ROM.Solution.q(:, t);
        end
        
        % Calculate maximum displacements
        max_uTail_FOM = max(uTail_FOM(2, :));
        max_uHead_FOM = max(uHead_FOM(1, :));
        max_uTail_ROM = max(uTail_ROM(2, :));
        max_uHead_ROM = max(uHead_ROM(1, :));
        
        % Calculate absolute differences and relative errors
        abs_diff_Tail = max_uTail_ROM - max_uTail_FOM;
        abs_diff_Head = max_uHead_ROM - max_uHead_FOM;
        rel_error_Tail = abs_diff_Tail / max_uTail_FOM;
        rel_error_Head = abs_diff_Head / max_uHead_FOM;
        
        % Store results in the results_matrix
        results_matrix(row_idx, :) = [num_elements, kActu, ...
            max_uTail_FOM, rel_error_Tail, rel_error_Head, ...
            timeFOMBuild, timeROMBuild, timeFOMSolve, timeROMSolve, timeFOM, timeROM];
        row_idx = row_idx + 1;
        
        % Save results to .mat file
        mat_filename = sprintf('A_results_%del_kActu_%.3f.mat', num_elements, kActu);
        save(mat_filename, 'timeFOM', 'timeROM', ...
            'timeFOMBuild', 'timeROMBuild','timeFOMSolve', 'timeROMSolve', ...
            'uTail_FOM', 'uTail_ROM', 'uHead_FOM', 'uHead_ROM', ...
            'max_uTail_FOM', 'max_uHead_FOM', 'max_uTail_ROM', 'max_uHead_ROM', ...
            'abs_diff_Tail', 'abs_diff_Head', 'rel_error_Tail', 'rel_error_Head');
        
        fprintf('Results saved to: %s\n', mat_filename);
        
        % Create and save figure comparing FOM and ROM displacements
        f = figure('units','centimeters','position',[3 3 9 6]);

        % x-position (Head)
        subplot(2,1,1);
        hold on;
        plot(timePlot, uHead_FOM(1,:), '--', 'DisplayName', 'FOM');
        plot(timePlot, uHead_ROM(1,:), 'DisplayName', 'ROM');
        grid on;
        xlabel('Time [s]');
        ylabel('Head x-position [m]');
        legend('Location','northwest', 'interpreter', 'latex');
    
        % y-position (Tail)
        subplot(2,1,2);
        hold on;
        plot(timePlot, uTail_FOM(2,:), '--', 'DisplayName', 'FOM');
        plot(timePlot, uTail_ROM(2,:), 'DisplayName', 'ROM');
        grid on;
        xlabel('Time [s]');
        ylabel('Tail y-position [m]');

        % Save figure
        fig_filename = sprintf('Displacement_Comparison_%del_kActu_%.3f.jpg', num_elements, kActu);
        exportgraphics(f, fig_filename, 'Resolution', 600);

        % Close the figure to avoid memory issues during iterations
        close(f);
        
    end
end

% The results_matrix now contains all results for all configurations.


%%
nElVec = [660,1272,4270,8086,24822];
timePerIt5 = [0.42,0.46,0.65,0.95,2.09];
timePerIt8 = [0.56,0.74,1.61,2.28,6.49];
f3 = figure('units','centimeters','position',[3 3 9 5]);
plot(nElVec,timePerIt5,'-x')
hold on
plot(nElVec,timePerIt8,'-*')
xlabel('Number of FE')
ylabel('Time per Build + EoMs')
legend('5 parameters', '8 parameters','Interpreter','latex', 'Location','North West')
grid on
hold off
exportgraphics(f3,'comp_time_analysis_FE.jpg','Resolution',600)

%%
nParamVec = [3,5,8];
timePerIt24822 = [1.81,2.09,6.49];
timePerIt8086 = [0.78,0.95,2.28];
timePerIt4270= [0.59,0.65,1.61];

f4 = figure('units','centimeters','position',[3 3 9 5]);
plot(nParamVec,timePerIt4270,'-x')
hold on
plot(nParamVec,timePerIt8086,'-*')
plot(nParamVec,timePerIt24822,'-o')
xlabel('Number of parameters')
ylabel('Time per Build + EoMs')
legend('4270 elements', '8086 elements', '24822 elements', 'Interpreter','latex', 'Location','North West')
grid on
hold off
exportgraphics(f4,'comp_time_analysis_param.jpg','Resolution',600)


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






