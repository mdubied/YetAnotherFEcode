% EXAMPLE: beam meshed with 3D element
clear; 
close all; 
clc

whichModel = 'CUSTOM'; % or "ABAQUS"


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
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 3;
Ly = .3;
nx = 30;
ny = 3;
switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'Job-BeamQuad';
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
end

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
switch upper( whichModel )
    case 'CUSTOM'
        myMesh.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
    case 'ABAQUS'
        myMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
end

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(myMesh);
M = BeamAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;
C= 0*M+0*K;

%% EXAMPLE 1                                                        

% % Eigenvalue problem_______________________________________________________
% n_VMs = 2; % first n_VMs modes with lowest frequency calculated 
% Kc = BeamAssembly.constrain_matrix(K);
% Mc = BeamAssembly.constrain_matrix(M);
% [V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
% [f0,ind] = sort(sqrt(diag(om))/2/pi);
% V0 = V0(:,ind);
% for ii = 1:n_VMs
%     V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
% end
% V0 = BeamAssembly.unconstrain_vector(V0);
% 
% % PLOT
% mod = 2;
% elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
% figure('units','normalized','position',[.2 .1 .6 .8])
% PlotMesh(nodes, elementPlot, 1);
% v1 = reshape(V0(:,mod), 2, []).';
% PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', Ly*1.1);
% title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])


%% VMs
 Kc = BeamAssembly.constrain_matrix(K);
 Mc = BeamAssembly.constrain_matrix(M);
 
 n_VMs = 15;
[VMs,f0,time_vm]=BeamAssembly.VMs_compute(n_VMs,1);
%figure('units','normalized','position',[.2 .1 .6 .8])

PNx=ceil(n_VMs/2); %Plot number
PNy=n_VMs-PNx;
if PNx==1
    PNx=2;
    PNy=1;
elseif PNx==2
    PNx=3;
    PNy=1;
end

% normalization
% for ii = 1:size(VMs,2)
%     VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
% end

% PLOT
for mod=3
% mod = 4;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
%subplot(4,3,mod)
figure('units','normalized','position',[.2 .1 .6 .8])

PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMs(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor',Ly*1.1);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])

end
close all
%% MDs
NMDs =15;
[MDs,MDs_names,time_md] = BeamAssembly.MDs_compute( n_VMs,NMDs,u0);

%   figure('units','normalized','position',[.2 .1 .6 .8])
PNxx=ceil(size(MDs,2)/2); %Plot number
PNyy=size(MDs,2)-PNxx;
if PNxx==1
    PNxx=2;
    PNyy=1;
elseif PNxx==2
    PNxx=3;
    PNyy=1;
end
% for ii = 1:size(MDs,2)
%    MDs(:,ii) = MDs(:,ii)/max(sqrt(sum(MDs(:,ii).^2,2)));
% end
m=3
for mod=1:m*(m+1)/2%size(MDs,2)
% mod = 4;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
  figure('units','normalized','position',[.2 .1 .6 .8])
%  subplot(2,3,mod)

 
PlotMesh(nodes, elementPlot, 0);
d1 = reshape(MDs(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, d1, 'factor', Ly*0.3 );
title(['MD' num2str(MDs_names{mod})])

end
 % Q3=BeamAssembly.tensor('tensor_Q3',[myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs],[2 3]);
% Q4=BeamAssembly.tensor('tensor_Q4',[myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs],[2 3 4]);

%deri1=BeamAssembly.stiffness_derivative(u0,VMs);
close all
%% EXAMPLE 2                                                        

% Define external force:
% % Body force
% Pressure = 1e6;
% F = Pressure*BeamAssembly.uniform_body_force();

% Nodal force
F = zeros(myMesh.nDOFs,1);
%F = zeros(length(myMesh.EBC.unconstrainedDOFs),1);
nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(2)) = 10e7;

u_lin = BeamAssembly.solve_system(K, F);
ULIN = reshape(u_lin,2,[]).';	% Linear response
u = static_equilibrium(BeamAssembly, F, 'display', 'iter-detailed');
UNL = reshape(u,2,[]).';        % Nonlinear response

fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ... 
    '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 5;
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,UNL,'factor',scale,'color','k');
colormap jet
title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])

figure('units','normalized','position',[.2 .1 .6 .8])
scale = 5;
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,ULIN,'factor',scale,'color','k');
colormap jet
title(['LINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])



%% Dynamic response using Implicit Newmark
% forcing frequency of the average of first two natural frequencies
omega_ext = 0.5*2*pi*f0(2)+0.5*2*pi*f0(1)+0.5*2*pi*f0(3); 
 %omega_ext=2*pi*f0(1);
T =  2*pi/omega_ext; % time period of forcing

% load amplification factor
amplification_factor = 1;

% forcing function
F_ext = @(t) amplification_factor * F * sin(omega_ext * t);

% Initial condition: equilibrium
u0 = zeros(BeamAssembly.Mesh.nDOFs, 1);
v0 = zeros(BeamAssembly.Mesh.nDOFs, 1);
a0 = zeros(BeamAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0)) 

q0 = BeamAssembly.constrain_vector(u0);
qd0 = BeamAssembly.constrain_vector(v0);
qdd0 = BeamAssembly.constrain_vector(a0);

% time step for integration
h = T/100;

% Precompute data for Assembly object
BeamAssembly.DATA.M = M;
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.C = C; %rayleigh

% Instantiate object for linear time integration
TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,BeamAssembly,F_ext);

% Linearized Time Integration
tmax = 50*T; 
%%
%tmax=0.002;
tic
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin.Solution.u = BeamAssembly.unconstrain_vector(TI_lin.Solution.q);
FullLDuration=toc
% Animate solution on Mesh (very slow)
%AnimateFieldonDeformedMesh(myMesh.nodes,myMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:2,'filename','lineardisp')
%%
% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamAssembly,F_ext);

% Nonlinear Time Integration

tic
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = BeamAssembly.unconstrain_vector(TI_NL.Solution.q);
FullNLduration=toc
%  save('TI_NL.mat','TI_NL');
% Generalized alpha scheme
% % linear
% TI_lin_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7, 'linear',true);
% TI_lin_alpha.Integrate(q0,qd0,qdd0,tmax,residual_lin);
% TI_lin_alpha.Solution.u = BeamAssembly.unconstrain_vector(TI_lin_alpha.Solution.q);
% 
% % nonlinear
% TI_NL_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
% TI_NL_alpha.Integrate(q0,qd0,qdd0,tmax,residual);
% TI_NL_alpha.Solution.u = BeamAssembly.unconstrain_vector(TI_NL_alpha.Solution.q);

%% Reduced solution 
LROM=[];
Ltime=[];
Lduration=[];
NLROM=[];
NLtime=[];
NLduration=[];

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

count=0;
VbMode=[];
MDeri=[];
for m=3%m=[1 3 7 10]
    for n=m*(m+1)/2 % n=[0 m*(m+1)/2]
% m= 3; % use the first five VMs in reduction
% n=m*(m+1)/2;
% n=0;
t=m+n;
% normalization
% for ii = 1:m
%     VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
% end
% for ii = 1:n
%    MDs(:,ii) = MDs(:,ii)/max(sqrt(sum(MDs(:,ii).^2,2)));
% end
% VMs = self.unconstrain_vector(VMs);
            
if n==0
 V = [VMs(:,1:m)];
else
   V = [VMs(:,1:m) MDs(:,1:n)];
end
%V=[ MDs(:,1:6)];

BeamReducedAssembly  = ReducedAssembly(myMesh,V);


BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix(0,0,u0);
BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
% 
q0 = zeros(t,1);
qd0 = zeros(t,1);
qdd0 = zeros(t,1);


TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Modal linear Residual evaluation function handle
Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,BeamReducedAssembly,F_ext);

% time integration
tic
TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
durationL=toc
Lduration=[Lduration durationL];
Ltime=[Ltime TI_lin_red.Solution.time'];
LROM=[LROM TI_lin_red.Solution.u(dof,:)'];

%% Reduced solution Noninear
% For demonstration purposes, we simply reduce the nonlinear system using
% out-of-plane bending modes. This is expected to produce bad results when 
% in-plane stretching is involved in the response.


%TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
TI_NL_alpha_red = ImplicitNewmark('timestep',h,'alpha',0.005);

% Modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,BeamReducedAssembly,F_ext);

% time integration
tic
TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
durationNL=toc

NLduration=[NLduration durationNL];

NLtime=[NLtime TI_NL_alpha_red.Solution.time'];
NLROM=[NLROM TI_NL_alpha_red.Solution.u(dof,:)'];

count=count+1;
VbMode=[VbMode m];
MDeri=[MDeri n];
    end
end
 
%% Reduced solution  Tensor Approach
LROM=[];
Ltime=[];
Lduration=[];
NLROM=[];
NLtime=[];
NLduration=[];
ASEMduration=[];

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

count=0;
VbMode=[];
MDeri=[];
for m=5
    for  n=m*(m+1)/2%n=[0 m*(m+1)/2]
% m= 3; % use the first five VMs in reduction
% n=m*(m+1)/2;
% n=0;
t=m+n;
% normalization
%  for ii = 1:m
%     VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
% end
% for ii = 1:n
%    MDs(:,ii) = MDs(:,ii)/max(sqrt(sum(MDs(:,ii).^2,2)));
% end
% VMs = self.unconstrain_vector(VMs);
            
if n==0
 V = [VMs(:,1:m)];
else
   V = [VMs(:,1:m) MDs(:,1:n)];
end
%V=[ MDs(:,1:6)];
 
 %V=BeamAssembly.unconstrain_vector(V);
 
 BeamReducedAssembly  = ReducedAssembly(myMesh,V);
 
 tensors_ROM = ROM_reduced_tensors(BeamAssembly, elements, V);
% Vr=V(BeamReducedAssembly.Mesh.EBC.unconstrainedDOFs(:)',:);
% Vr=V;

% 
  BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
  BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix(0,0,u0);
  BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
% Mr=BeamAssembly.mass_matrix();
% Kr=BeamAssembly.stiffness_matrix(u0);
% Cr=BeamAssembly.damping_matrix(0,0,u0);
% % 
% Mrr= BeamAssembly.constrain_matrix(Mr);
% Krr=BeamAssembly.constrain_matrix(Kr);
% Crr=BeamAssembly.constrain_matrix(Cr);
% % % 
% Mrrr = Vr'*Mrr*Vr;
% Krrr = Vr'*Kc*Vr;
% Crrr =  Vr'*Crr*Vr;

tic


Q2=BeamReducedAssembly.DATA.K;
 % Q2=Krrr;
Q3r=BeamReducedAssembly.tensor('tensor_Q3',[n+m n+m n+m],[]);
Q4r=BeamReducedAssembly.tensor('tensor_Q4',[n+m n+m n+m n+m],[]);
% Q3 = BeamReducedAssembly.constrain_tensor(Q3r);
% Q4 =BeamReducedAssembly.constrain_tensor(Q4r);
Q3=full(Q3r);
Q4=full(Q4r);
Q2=tensors_ROM.Q2;
Q3=tensors_ROM.Q3;
Q4=tensors_ROM.Q4;

%% % parametric formulation for defects (only to test) 

% FORMULATION = 'N1'; % N1/N1t/N0
% VOLUME = 1;         % integration over defected (1) or nominal volume (0)
% U = V(:,1:2);       % defect basis
% tensors_DpROM = DpROM_reduced_tensors(FORMULATION, VOLUME, ...
%     BeamAssembly, elements, V, U);
% 
% % evaluate the defected tensors at xi
% xi = rand(size(U,2),1);
% [Q2, Q3, Q4] = DefectedTensors(tensors_DpROM, xi);

%%

durationTensor=toc;
ASEMduration=[ASEMduration durationTensor];
%% this also works (moving from full to reduce tensor)
%  VV=sptensor(Vr);
%  VVt=sptensor(Vr');
% % 
%  Q3nr=BeamAssembly.tensor('tensor_Q3',[myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs],[]);
%  Q4nr=BeamAssembly.tensor('tensor_Q4',[myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs],[]);
%     Q3 = BeamAssembly.constrain_tensor(Q3nr);
%     Q4 =BeamAssembly.constrain_tensor(Q4nr);
%  Q3  = ttt(sptensor(ttt(ttt(VVt,Q3,2,1),VV,3,1)),VV ,2,1);
%  Q4   = ttt(sptensor(ttt(ttt(ttt(VVt,Q4,2,1),VV,4,1),VV ,3,1)),VV ,2,1);
%testing moving from full to reduce
%%

Q3t = Q3 + permute(Q3,[1 3 2]); 
Q4t = Q4 + permute(Q4,[1 3 2 4]) + permute(Q4,[1 4 2 3]);

% q0 = zeros(t,1);
% qd0 = zeros(t,1);
% qdd0 = zeros(t,1);
q0=V'*u0;
qd0=V'*v0;
qdd0=V'*a0;
TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Modal linear Residual evaluation function handle
Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,BeamReducedAssembly,F_ext);

% time integration
tic
TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
durationL=toc
Lduration=[Lduration durationL];
Ltime=[Ltime TI_lin_red.Solution.time'];
LROM=[LROM TI_lin_red.Solution.u(dof,:)'];

% Reduced solution Noninear
% For demonstration purposes, we simply reduce the nonlinear system using
% out-of-plane bending modes. This is expected to produce bad results when 
% in-plane stretching is involved in the response.

%%

TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
% Modal nonlinear Residual evaluation function handle
 Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,BeamReducedAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
 %Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor_diffApp(q,qd,qdd,t,BeamAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t,Vr,Mrrr,Crrr);

% time integration
tic
TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_newmark_red);
TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
durationNL=toc

NLduration=[NLduration durationNL];

NLtime=[NLtime TI_NL_newmark_red.Solution.time'];
NLROM=[NLROM TI_NL_newmark_red.Solution.u(dof,:)'];

%% TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
% 
% 
% TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
% % Modal nonlinear Residual evaluation function handle
% Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,BeamReducedAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
% 
% % time integration
% tic
% TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
% TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
% durationNL=toc
% 
% NLduration=[NLduration durationNL];
% 
% NLtime=[NLtime TI_NL_alpha_red.Solution.time'];
% NLROM=[NLROM TI_NL_alpha_red.Solution.u(dof,:)'];
%%




count=count+1;
VbMode=[VbMode m];
MDeri=[MDeri n];
    end
end
%%
nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

figure;

plot(TI_lin.Solution.time, TI_lin.Solution.u(dof,:),'DisplayName', ['Full linear (Newmark) ' 'time= ' num2str(FullLDuration) ])
hold on
plot(TI_NL.Solution.time, TI_NL.Solution.u(dof,:),'DisplayName', ['Full nonlinear (Newmark) ' 'time= ' num2str(FullNLduration)])
hold on
for i=1:count
Vmode=VbMode(i);
Mderi=MDeri(i);


plot(Ltime(:,i), LROM(:,i),'DisplayName', ['Reduced linear (Newmark) ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(Lduration(i)) 's'])

xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')

plot(NLtime(:,i), NLROM(:,i),'DisplayName', ['Reduced nonlinear (Newmark) T' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(NLduration(i)) 's'])
xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')

end
