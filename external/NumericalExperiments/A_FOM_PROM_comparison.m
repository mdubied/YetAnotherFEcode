% ------------------------------------------------------------------------ 
% Accuracy and computational speed comparison.
% 
% Last modified: 05/05/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 1;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');

%% PREPARE MODEL __________________________________________________________                                                   

% DATA ____________________________________________________________________
E       = 260000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.4;        % Poisson's ratio 

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain
myElementConstructor = @()Tet4Element(myMaterial);

% MESH ____________________________________________________________________

filename = %3d_rectangle_660el'; %3d_fish_for_mike';'3d_rectangle_2385el' % need to set 0.2*k
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
nel = size(elements,1);

% convert to cm to m and reduce the initial y dimension
nodes = nodes*0.01;
nodes(:,2)=0.8*nodes(:,2);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

MeshNominal_FOM = Mesh(nodes);
MeshNominal_FOM.create_elements_table(elements,myElementConstructor);

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
fixedPortion = 0.6;
nset1 = {};
fixedElements = zeros(nel,1);
for el=1:nel   
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterX >= -Lx*fixedPortion
        for n=1:size(elements,2) 
            if  ~any(cat(2, nset1{:}) == elements(el,n))
                nset1{end+1} = elements(el,n); 
            end
        
    end   
        fixedElements(el)=1;
    end
    
%     for n=1:size(elements,2)
%         
%         if  nodes(elements(el,n),1)>-Lx*fixedPortion && ~any(cat(2, nset1{:}) == elements(el,n))
%             nset1{end+1} = elements(el,n); 
%             fixedElements(el)=1;
%         end
%         
%     end   
end

for l=1:length(nset1)
    MeshNominal.set_essential_boundary_condition([nset1{l}],1:3,0)  % all DOFs constrained to get VMs. Rigid body modes are added in build_ROM
    MeshNominal_FOM.set_essential_boundary_condition([nset1{l}],2:3,0)
end  

 


% shape variations for PROM
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);
U = [z_tail,y_head,y_thinFish];


%% FIGURE A1 (muscles, rigid part, VM) ____________________________________
% Note: the position of the muscle is defined in the build_ROM/FOM/PROM
% functions. The number of VMs used also. Only the constraints for the
% rigid par of the fish is defined above (boundary conditions)

% MUSCLES _________________________________________________________________
% left muscle
leftMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*0.6 && elementCenterX > -Lx
        leftMuscle(el) = 1;
    end    
end

% right muscle
rightMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*0.6 && elementCenterX > -Lx
        rightMuscle(el) = 1;
    end    
end

% VIBRATION MODES _________________________________________________________
% get first vibration mode
NominalAssemblyForPlot = Assembly(MeshNominal);
Mn = NominalAssemblyForPlot.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssemblyForPlot.tangent_stiffness_and_force(u0);
% store matrices
NominalAssemblyForPlot.DATA.K = Kn;
NominalAssemblyForPlot.DATA.M = Mn;

% vibration modes
n_VMs = 1;
Kc = NominalAssemblyForPlot.constrain_matrix(Kn);
Mc = NominalAssemblyForPlot.constrain_matrix(Mn);
[VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssemblyForPlot.unconstrain_vector(VMn);

% FIGURE __________________________________________________________________
f_A1 = figure('units','centimeters','position',[3 3 9 3.5]);
pos1 = [0,0,0.5,1];
pos2 = [0.5,0,0.5,1];

% subplot 1: muscles' placement
ax1 = subplot(1,2,1,'Position',pos1);

Plot2MusclesAndConstraints(nodes,elements, ...
    leftMuscle,'green',rightMuscle,'blue', ...
    fixedElements,'red');

% subplot2: VM1
ax2 = subplot(1,2,2,'Position',pos2);
elementPlot = elements(:,1:4);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);
v1 = reshape(-VMn(:,1), 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));

axis([ax1 ax2],[-0.39 0 -0.08 0.08 -0.11 0.11])
exportgraphics(f_A1,'A_muscles_placement_VM.pdf','Resolution',1400)

%% SIMULATION PARAMETERS __________________________________________________
h = 0.01;
tmax = 2;


%% FOM ____________________________________________________________________
tStartFOM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building FOM ... \n')
[Assembly,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_FOM_3D(MeshNominal_FOM,nodes,elements);   

% solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_FOM = solve_EoMs_FOM(Assembly,elements,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax); 
toc

fprintf('Time needed to solve the problem using FOM: %.2fsec\n',toc(tStartFOM))
timeFOM = toc(tStartFOM);

%% ROM ____________________________________________________________________
tStartROM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')
[V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_ROM_3D(MeshNominal,nodes,elements,USEJULIA);  

% solve EoMs 
tic
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_ROM = solve_EoMs(V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax); 
toc

fprintf('Time needed to solve the problem using ROM: %.2fsec\n',toc(tStartROM))
timeROM = toc(tStartROM);

%% PROM ___________________________________________________________________
tStartPROM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building PROM ... \n')
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_PROM_3D(MeshNominal,nodes,elements,U,USEJULIA,VOLUME,FORMULATION);      

% solve EoMs (with sensitivities for the PROM)
tic 
fprintf('____________________\n')
fprintf('Solving EoMs and sensitivities...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax); 
toc

fprintf('Time needed to solve the problem using PROM: %.2fsec\n',toc(tStartPROM))
timePROM = toc(tStartPROM);


%% PLOT ___________________________________________________________________
uTail_FOM = zeros(3,tmax/h);
uTail_ROM = zeros(3,tmax/h);
uTail_PROM = zeros(3,tmax/h);
sol_FOM = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));

for t=1:tmax/h
    uTail_FOM(:,t) = sol_FOM(tailProperties.tailNode*3-2:tailProperties.tailNode*3,t);
    uTail_ROM(:,t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_ROM.Solution.q(:,t);  
%     uTail_PROM(:,t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,t);
end

f_A2 = figure('units','centimeters','position',[3 3 9 6]);
% x-position
subplot(2,1,1);
hold on
plot(timePlot,x0Tail+uTail_FOM(1,:),'--','DisplayName','FOM')
plot(timePlot,x0Tail+uTail_ROM(1,:),'DisplayName','ROM')
% plot(timePlot,x0Tail+uTail_PROM(1,:),'DisplayName','PROM')
hold on
grid on
xlabel('Time [s]')
ylabel('x-position [m]')
legend('Location','northwest', 'interpreter','latex')

% y-position
subplot(2,1,2);
hold on
plot(timePlot,uTail_FOM(2,:),'--','DisplayName','FOM')
plot(timePlot,uTail_ROM(2,:),'DisplayName','ROM')

% plot(timePlot,uTail_PROM(2,:),'DisplayName','PROM')
grid on
ylim([-0.04,0.04])
xlabel('Time [s]')
ylabel('y-position [m]')
% legend('Location','southwest', 'interpreter','latex')
exportgraphics(f_A2,'A_FOM_vs_ROM_2385el.pdf','Resolution',1400)

%% ANIMATION ______________________________________________________________
elementPlot = elements(:,1:4); 
nel = size(elements,1);

% top muscle
topMuscle = zeros(nel,1);

for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        topMuscle(el) = 1;
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        bottomMuscle(el) = 1;
    end    

end

actuationValues = zeros(size(TI_NL_ROM.Solution.u,2),1);
for t=1:size(TI_NL_ROM.Solution.u,2)
    actuationValues(t) = 0;
end

actuationValues2 = zeros(size(TI_NL_ROM.Solution.u,2),1);
for t=1:size(TI_NL_ROM.Solution.u,2)
    actuationValues2(t) = 0;
end
sol = TI_NL_ROM.Solution.u(:,1:2:end);
AnimateFieldonDeformedMeshActuation2Muscles(nodes, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,sol, ...
    'factor',1,'index',1:3,'filename','result_video','framerate',1/h*0.5)






