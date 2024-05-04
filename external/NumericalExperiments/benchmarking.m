% ------------------------------------------------------------------------ 
% Accuracy and computational speed comparison.
% 
% Last modified: 04/02/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 1;

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

filename = '3d_fish_for_mike';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
nel = size(elements,1);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil

% plot nominal mesh
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off

% boundary conditions of nominal mesh
nset = {};
for el=1:nel   
    for n=1:size(elements,2)
        if  ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
%         if  nodes(elements(el,n),1)> fixedPoint-marginFixedPoint && ...
%                 nodes(elements(el,n),1)< fixedPoint+marginFixedPoint && ...
%                 ~any(cat(2, nset{:}) == elements(el,n))
%             nset{end+1} = elements(el,n);
%         end
    end   
end

for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],3,0)
end   



[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% for testing
U = [z_tail,y_head];

%% SIMULATION PARAMETERS __________________________________________________
h = 0.005;
tmax = 4;

%% ROM ____________________________________________________________________
tStartROM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_ROM_3D(MeshNominal,nodes,elements,USEJULIA);  

% solve EoMs 
tic 
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_ROM = solve_EoMs(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax); 
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

%% FOM ____________________________________________________________________
tStartFOM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building FOM ... \n')
[Assembly,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_FOM_3D(MeshNominal,nodes,elements);   
%%
% solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs and sensitivities...\n') 
TI_NL_FOM = solve_EoMs_FOM(Assembly,elements,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax); 
toc

fprintf('Time needed to solve the problem using FOM: %.2fsec\n',toc(tStartFOM))
timeFOM = toc(tStartFOM);

%%



%% PLOT ___________________________________________________________________
uTail_ROM = zeros(3,tmax/h);
uTail_PROM = zeros(3,tmax/h);
uTail_FOM = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for t=1:tmax/h
    uTail_ROM(:,t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_ROM.Solution.q(:,t);
    uTail_PROM(:,t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,t);
    uTail_FOM(:,t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_FOM.Solution.q(:,t);

end

figure
subplot(2,1,1);
plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
hold on
xlabel('Time [s]')
ylabel('x-position tail node')
legend('Location','northwest')
subplot(2,1,2);
plot(timePlot,uTail(2,:),'DisplayName','k=0')
hold on
xlabel('Time [s]')
ylabel('y-position tail node')
legend('Location','southwest')
drawnow

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






