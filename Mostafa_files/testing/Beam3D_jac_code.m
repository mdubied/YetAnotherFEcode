% EXAMPLE: MEMS resonator dynamic response

clear ; close all; clc

addpath MODELs FE_routines MDs_routines DDs_routines tensors_routines
addpath tensor_toolbox utilities

% Load input file (stored in .mat as "model_filename.mat") ________________
meshfile='MeshBeam3D_HEX20.dat';
model=Ansys2Matlab_mesh(meshfile,'HEX20');
% BCs _____________________________________________________________________
% Nset = 1;
% fixed_nodes = model.sets{1};    % node set coming from ABAQUS-defined sets
% fixdofs = model.idb(fixed_nodes,2:4);
% fixdofs = fixdofs(:);
% model.dofs.fixed = fixdofs;
% model.dofs.free = setdiff(model.dofs.all,model.dofs.fixed);
% plotNodeSet(model,Nset,'r*')
fixed_nodes = model.fixedNodes(:)';   
fixdofs = model.idb(fixed_nodes,2:end);               %if 3D put 2:4 2D 2:3 all this is useless in YETANOTHERFEMCODE
fixdofs = fixdofs(:);
model.dofs.fixed = fixdofs;
model.dofs.free = setdiff(model.dofs.all,model.dofs.fixed);

% K and M _________________________________________________________________
ngauss = 2;
NLgeom = 1;
U0 = 0;
K = FE_assemblyKF(model,U0,ngauss,NLgeom);
M = FE_assemblyM(model,ngauss);

%%
% VMs and MDs _____________________________________________________________
Nm = 4;
[VMs,f0,time_vm] = VMs_compute(model,K,M,Nm,1);
NMDs = 3;
[MDs,MDs_names,time_md] = MDs_compute(model,K,VMs,NMDs,ngauss);

% normalize bases
VMs = ROBnormalization(VMs,'displacement');
MDs = ROBnormalization(MDs,'displacement');

% plot VMs and MDs
plotVMs = 1;
plotMDs = 1;
scale = 25;
vis = 'U';      % U, U1, U2, U3
flat = 'no';    % flat
AXvm = plot_VMs(model,VMs,f0,plotVMs,scale,vis,flat);
AXmd = plot_MDs(model,MDs,NMDs,plotMDs,scale,vis,flat);

% compute ROM _____________________________________________________________
V = [VMs MDs];
Q2 = V'*K(model.dofs.free,model.dofs.free)*V;
[Q3,Q4,time_Qs] = tensors_Q34_assembly_classic(model,V,ngauss);
%%
% external force vector ___________________________________________________
P1 = [0 -100 10];       % forced node
fdir = 2;               % direction (1/2/3)
AMP = 750;              % force amplitude
omega0 = 2*pi*f0(1);    % force circular frequency
ncycles = 3;            % force cycles

nf = findnode(P1(1),P1(2),P1(3),model.nodes);
fext = zeros(model.Ndofs,1);
fdof = nf*3-(3-fdir);
fext(fdof) = -1;

% plot force
figure(figmesh)
P2 = P1;
P2(fdir) = P2(fdir)-max(ylim)/3;
arrow3(P1,P2,'k-3',2,2,2,.9)
axis tight; drawnow
%%
% Time integration (ROM) __________________________________________________
u0 = zeros(model.Ndofs,1);
    matrices.forcing.fvector    	= fext;
    matrices.forcing.loadfunction	= @(t) AMP*sin(omega0*t);
    matrices.IC.displacements    	= u0;
    matrices.IC.velocities       	= u0;
    matrices.Qtensors.Q2         	= Q2;
    matrices.Qtensors.Q3          	= Q3;
    matrices.Qtensors.Q4        	= Q4;
    matrices.Mass                	= M;
    matrices.ROB                	= V;
    settings.damping.alpha     = 0;
    settings.damping.beta      = 0;
    settings.time.ntimesteps   = ncycles*50;
    settings.time.Tend         = ncycles*2*pi/omega0;
    settings.solver.tolerance  = 1e-6;
    settings.solver.nmax       = 15;
    settings.PRINTinfo         = 1;
    settings.Uout              = 1;
output_rom = Newmark_ROM(model,matrices,settings);

% linear
    matrices_lin = matrices;
    matrices_lin.Qtensors.Q3 = Q3*0;
    matrices_lin.Qtensors.Q4 = Q4*0;
output_rom_lin = Newmark_ROM(model,matrices_lin,settings);

% Time integration (full) _________________________________________________
    matrices.K0 = K;
    matrices.backup = [];
    settings.solver.ngauss = ngauss;
    settings.monitor = 0;
    settings.interp = 1;
    settings.sim_number = 01;
    settings.savepath = [cd '/SCRATCH/'];
output_full = Newmark_FULL_c(model,matrices,settings);

% Plot response ___________________________________________________________

figure('units','normalized','position',[.2 .2 .6 .6])
plot(output_full.ti,output_full.u(fdof,:),'b-',...
    'linewidth',2,'DisplayName','full')
hold on
plot(output_rom.t,output_rom.u(fdof,:),'r--',...
    'linewidth',2,'DisplayName','ROM')
plot(output_rom_lin.t,output_rom_lin.u(fdof,:),'k:',...
    'linewidth',2,'DisplayName','ROM linear')
grid on;
xlabel('time [s]','interpreter','latex','fontsize',20)
ylabel('y-displacement $\mu$','interpreter','latex','fontsize',20)
hl = legend;
hl.Interpreter = 'latex';
hl.FontSize = 20;




