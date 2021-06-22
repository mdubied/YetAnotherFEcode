% EXAMPLE: beam meshed with 2D element
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
%% old use of VMS

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
Cc=BeamAssembly.constrain_matrix(C);
n_VMs = 3;
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

%normalization
for ii = 1:size(VMs,2)
    VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
end

% PLOT
for mod=1:5
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    %subplot(4,3,mod)
    figure('units','normalized','position',[.2 .1 .6 .8])
    
    PlotMesh(nodes, elementPlot, 0);
    v1 = reshape(VMs(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor',Ly*1.1);
    title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
    
end
%% MDs
NMDs =n_VMs*(n_VMs+1)/2;
%[MDs,MDs_names,time_md] = BeamAssembly.MDs_compute( n_VMs,NMDs,u0);
[MDs,MDs_names]= modal_derivatives(BeamAssembly,elements,VMs);  
%using%julia
%   figure('units','normalized','position',[.2 .1 .6 .8])
% PNxx=ceil(size(MDs,2)/2); %Plot number
% PNyy=size(MDs,2)-PNxx;
% if PNxx==1
%     PNxx=2;
%     PNyy=1;
% elseif PNxx==2
%     PNxx=3;
%     PNyy=1;
% end
for ii = 1:size(MDs,2)
    MDs(:,ii) = MDs(:,ii)/max(sqrt(sum(MDs(:,ii).^2,2)));
    %MDss(:,ii)=MDss(:,ii)/max(sqrt(sum(MDss(:,ii).^2,2)));
end

for mod=1%:m*(m+1)/2%size(MDs,2)
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    figure('units','normalized','position',[.2 .1 .6 .8])
    %  subplot(2,3,mod)
    PlotMesh(nodes, elementPlot, 0);
    d1 = reshape(MDs(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes, elementPlot, d1, 'factor', Ly*0.3 );
    %title(['MD:' MDs_names{mod}])
    title(['\theta_{' num2str(MDs_names(mod,1)) num2str(MDs_names(mod,2)) '}'])
    %d1s = reshape(MDss(:,mod), 2, []).';
    %PlotFieldonDeformedMesh(nodes, elementPlot, d1s, 'factor', Ly*.5 );
end
% Q3=BeamAssembly.tensor('tensor_Q3',[myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs],[2 3]);
% Q4=BeamAssembly.tensor('tensor_Q4',[myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs myMesh.nDOFs],[2 3 4]);

%deri1=BeamAssembly.stiffness_derivative(u0,VMs);
%% Force

% Define external force:
% % Body force
% Pressure = 1e6;
% F = Pressure*BeamAssembly.uniform_body_force();

% Nodal force
F = zeros(myMesh.nDOFs,1);
%F = zeros(length(myMesh.EBC.unconstrainedDOFs),1);
nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(2)) = 10e6;
%% - Static Analysis

% u_lin = BeamAssembly.solve_system(K, F);
% ULIN = reshape(u_lin,2,[]).';	% Linear response
% u = static_equilibrium(BeamAssembly, F, 'display', 'iter-detailed');
% UNL = reshape(u,2,[]).';        % Nonlinear response
%
% fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ...
%     '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))
%
% % PLOT
% figure('units','normalized','position',[.2 .1 .6 .8])
% scale = 5;
% PlotMesh(nodes, elementPlot, 0);
% PlotFieldonDeformedMesh(nodes,elementPlot,UNL,'factor',scale,'color','k');
% colormap jet
% title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])
%
% figure('units','normalized','position',[.2 .1 .6 .8])
% scale = 5;
% PlotMesh(nodes, elementPlot, 0);
% PlotFieldonDeformedMesh(nodes,elementPlot,ULIN,'factor',scale,'color','k');
% colormap jet
% title(['LINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])
%% Dynamic response Initialiazing
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
tmax = 100*T;
%% Linear Dynamic Response
%tmax=0.002;
tic
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin.Solution.u = BeamAssembly.unconstrain_vector(TI_lin.Solution.q);
FullLDuration=toc
% Animate solution on Mesh (very slow)
%AnimateFieldonDeformedMesh(myMesh.nodes,myMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:2,'filename','lineardisp')

%% NL Dynamic Response
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
f0c_ROM={};

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

count=0;
VbMode=[];
MDeri=[];
for m=3;% m=[2 3 7 10]
    for n=6 % n=[0 m*(m+1)/2]
        
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
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamAssembly.mass_matrix)*V(:,ii));
        end
        BeamReducedAssembly  = ReducedAssembly(myMesh,V);
        
        
        BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
        BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix(0,0,u0);
        BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
        
        f0_ROM_i= sort(sqrt(eig(BeamReducedAssembly.DATA.M\BeamReducedAssembly...
            .DATA.K))/2/pi) ;
        f0c_ROM=[f0c_ROM f0_ROM_i];
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
[max_size, max_index] = max(cellfun('size', f0c_ROM, 1));
f0_ROM=zeros(max_size,size(f0c_ROM,2));
for ii=1:max_index
    f0c=f0c_ROM{ii};
    s=length(f0c);
    f0_ROM(1:s,ii)=f0c;
end
%% HyperReduced Solution
LROMH=[];
LtimeH=[];
LdurationH=[];
NLROMH=[];
NLtimeH=[];
NLdurationH=[]; 
f0c_ROMH={};

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

countH=0;
VbModeH=[];
MDeriH=[];
for m=3;% m=[2 3 7 10]
    for n=6 % n=[0 m*(m+1)/2]
        
        t=m+n;
      
%%
        
        if n==0
            V = [VMs(:,1:m)];
        else
            V = [VMs(:,1:m) MDs(:,1:n)];
        end
        %V=[ MDs(:,1:6)];
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamAssembly.mass_matrix)*V(:,ii));
        end
        BeamReducedAssembly  = ReducedAssembly(myMesh,V);
        
        
        BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
        BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix(0,0,u0);
        BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
        
        %% algorithm 3
        V_H=V(:,1:m)
        Lin_sol=TI_lin.Solution.u;
        Lin_sol_snap=[];
        for ii=1:size(Lin_sol,2)
            if rem(ii,10)==0
                Lin_sol_snap=[Lin_sol_snap Lin_sol(:,ii)];
            end
        end
        eta=V_H'*Lin_sol_snap;
        
        Theta = QM_Theta_from_SMDs(V(:,m+1:m+n), MDs_names(1:n,:));
        %QM uplifting
        uu = zeros(size(V_H,1), size(eta, 2));
        for tt = 1 : size(eta, 2)
            uu(:,tt) = einsum('Iij,iJ,jK->IJK', Theta, eta(:,tt), eta(:,tt));
        end
        u_lin_ECSW = V_H * eta + 1/2*uu;
     
        %Construct Gb
        qq=(V'*V)^-1*V'*u_lin_ECSW;
        tic
        [G,b]=BeamReducedAssembly.constructGb(qq);
        GBconstructTime=toc;
        %%fNNLS
        tic
        [x_fnnls,w_fnnls]=fnnls(G'*G,G'*b,norm(b)*0.01);
        fnnlsTime=toc
        tic;[x_nnls,w_nnls,info]=nnls(full(G),full(b),struct('Accy',1,'Tol',1/((norm(b)*5))));
        nnlTime=toc
        options = optimset('TolX',1/(norm(b)*0.01));
        tic;x_lsq=lsqnonneg(G,b,options);
        lsqTime=toc
        nnz(x_fnnls)
        nnz(x_nnls)
        nnz(x_lsq)
         x_sNNLS=sNNLS(G,b,0.0035);
         nnz(x_sNNLS)
        BeamReducedAssembly.DATA.elementWeights=x_sNNLS;
        %%
        f0_ROM_i= sort(sqrt(eig(BeamReducedAssembly.DATA.M\BeamReducedAssembly...
            .DATA.K))/2/pi) ;
        f0c_ROMH=[f0c_ROMH f0_ROM_i];
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
        LdurationH=[LdurationH durationL];
        LtimeH=[LtimeH TI_lin_red.Solution.time'];
        LROMH=[LROMH TI_lin_red.Solution.u(dof,:)'];
        
        %% Reduced solution Noninear
        % For demonstration purposes, we simply reduce the nonlinear system using
        % out-of-plane bending modes. This is expected to produce bad results when
        % in-plane stretching is involved in the response.
        
        
        %TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
        TI_NL_alpha_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hyper(q,qd,qdd,t,BeamReducedAssembly,F_ext);
        
        % time integration
        tic
        TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
        durationNL=toc
        
        NLdurationH=[NLdurationH durationNL];
        
        NLtimeH=[NLtimeH TI_NL_alpha_red.Solution.time'];
        NLROMH=[NLROMH TI_NL_alpha_red.Solution.u(dof,:)'];
        
        countH=countH+1;
        VbModeH=[VbModeH m];
        MDeriH=[MDeriH n];
    end
end
[max_size, max_index] = max(cellfun('size', f0c_ROMH, 1));
f0_ROM=zeros(max_size,size(f0c_ROMH,2));
for ii=1:max_index
    f0c=f0c_ROMH{ii};
    s=length(f0c);
    f0_ROM(1:s,ii)=f0c;
end
%% Reduced solution  Tensor Approach (normal ROM) (name of variable ***T)
LROMT=[];
LtimeT=[];
LdurationT=[];
NLROMT=[];
NLtimeT=[];
NLdurationT=[];
ASEMdurationT=[];
f0c_ROMT={};
nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

countT=0;
VbModeT=[];
MDeriT=[];

for m=3 %m=[2 3 7 10]
    for n=6% n=[0 m*(m+1)/2]
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
        
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamAssembly.mass_matrix)*V(:,ii));
        end
        
        %V=BeamAssembly.unconstrain_vector(V);
        
        BeamReducedAssembly  = ReducedAssembly(myMesh,V);
        
        
        % Vr=V(BeamReducedAssembly.Mesh.EBC.unconstrainedDOFs(:)',:);
        % Vr=V;
        
        %
        BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
        BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix(0,0,u0);
        BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
        
        
        
        %% --------use of T in yetAnotherFEcode-------------------------------------
        % Q2=BeamReducedAssembly.DATA.K;
        %  % Q2=Krrr;
        % Q3r=BeamReducedAssembly.tensor('tensor_Q3',[n+m n+m n+m],[]);
        % Q4r=BeamReducedAssembly.tensor('tensor_Q4',[n+m n+m n+m n+m],[]);
        % % Q3 = BeamReducedAssembly.constrain_tensor(Q3r);
        % % Q4 =BeamReducedAssembly.constrain_tensor(Q4r);
        % Q3=full(Q3r);
        % Q4=full(Q4r);
        %-------------------------------------------------------------------------
        %% using Julia
        tic
        
        tensors_ROM = reduced_tensors_ROM(BeamAssembly, elements, V);
        Q2=tensors_ROM.Q2;
        Q3=tensors_ROM.Q3;
        Q4=tensors_ROM.Q4;
        durationTensor=toc;
        ASEMdurationT=[ASEMdurationT durationTensor];
        %compute Natural Freq
        f0_ROM_i= sort(sqrt(eig(BeamReducedAssembly...
            .DATA.M\Q2))/2/pi) ;
        
        f0c_ROMT=[f0c_ROMT f0_ROM_i];
        
        
        %% % parametric formulation for defects (only to test)---------------------
        
        % FORMULATION = 'N1'; % N1/N1t/N0
        % VOLUME = 1;         % integration over defected (1) or nominal volume (0)
        % U = V(:,1:2);       % defect basis
        % tensors_DpROM = reduced_tensors_DpROM(FORMULATION, VOLUME, ...
        %     BeamAssembly, elements, V, U);
        %
        % % evaluate the defected tensors at xi
        % xi = rand(size(U,2),1);
        % [Q2, Q3, Q4] = DefectedTensors(tensors_DpROM, xi);
        %% this also works (moving from full to reduce tensor) using yet antoher FE
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
        %% L sim
        
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
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear_tensor(q,qd,qdd,t,BeamReducedAssembly,F_ext,Q2);
        
        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        LdurationT=[LdurationT durationL];
        LtimeT=[LtimeT TI_lin_red.Solution.time'];
        LROMT=[LROMT TI_lin_red.Solution.u(dof,:)'];
        
        
        
        %% NL
        
        TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,BeamReducedAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
        %Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor_diffApp(q,qd,qdd,t,BeamAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t,Vr,Mrrr,Crrr);
        
        % time integration
        tic
        TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_newmark_red);
        TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
        durationNL=toc
        
        NLdurationT=[NLdurationT durationNL];
        
        NLtimeT=[NLtimeT TI_NL_newmark_red.Solution.time'];
        NLROMT=[NLROMT TI_NL_newmark_red.Solution.u(dof,:)'];
        
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
        % NLdurationT=[NLdurationT durationNL];
        %
        % NLtimeT=[NLtimeT TI_NL_alpha_red.Solution.time'];
        % NLROMT=[NLROMT TI_NL_alpha_red.Solution.u(dof,:)'];
        %%
        countT=countT+1;
        VbModeT=[VbModeT m];
        MDeriT=[MDeriT n];
        
        
    end
end

[max_sizeT, max_indexT] = max(cellfun('size', f0c_ROMT, 1));
f0_ROMT=zeros(max_sizeT,size(f0c_ROMT,2));

for ii=1:max_indexT
    f0c=f0c_ROMT{ii};
    s=length(f0c);
    f0_ROMT(1:s,ii)=f0c;
end
%% HFM-d init (Defected mesh) first VM
U=[VMs(:,1)];      %defect Basis
xi=[0.05];
U1 = reshape(U(:,1), 2, []).'*xi(1);
%U2=reshape(U(:,2), 2, []).'*xi(2);
Ut=U1%+U2;
nodes_defected=nodes+ Ut;%U1*xi;
% Mesh
MeshDefected=Mesh(nodes_defected);
MeshDefected.create_elements_table(elements,myElementConstructor);
MeshDefected.set_essential_boundary_condition([nset{1} nset{3}],1:2,0);
% Assembly
BeamDefectedAssembly = Assembly(MeshDefected);
Md = BeamDefectedAssembly.mass_matrix();
[Kd,~] = BeamDefectedAssembly.tangent_stiffness_and_force(u0);
% store matrices
BeamDefectedAssembly.DATA.K = Kd;
BeamDefectedAssembly.DATA.M = Md;
BeamDefectedAssembly.DATA.C= 0*Md+0*Kd;

% Vm and Md computation for defected Mesh
n_Vmd=7;
% Vibration Modes (VM): defected------------------------------------------
Kdc = BeamDefectedAssembly.constrain_matrix(Kd);
Mdc = BeamDefectedAssembly.constrain_matrix(Md);
[VMd,om] = eigs(Kdc, Mdc, n_Vmd, 'SM');
[f0d,ind] = sort(sqrt(diag(om))/2/pi);
VMd = VMd(:,ind);
for ii = 1:n_Vmd
    VMd(:,ii) = VMd(:,ii)/max(sqrt(sum(VMd(:,ii).^2,2)));
end
VMd = BeamAssembly.unconstrain_vector(VMd);

% Modal Derivative Defected------------------------------------------------
[MDd, MDd_names] = modal_derivatives(BeamDefectedAssembly, elements, VMd);
for ii = 1:size(MDd,2)
    MDd(:,ii) = MDd(:,ii)/max(sqrt(sum(MDd(:,ii).^2,2)));
end
% ploting VMd------------------------------------------------------
for mod=1
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    %subplot(4,3,mod)
    figure('units','normalized','position',[.2 .1 .6 .8])
    
    PlotMesh(nodes_defected, elementPlot, 0);
    v1 = reshape(VMd(:,mod), 2, []).';
    factor = 2*max(nodes_defected(:,2));
    PlotFieldonDeformedMesh(nodes_defected, elementPlot, v1, 'factor', factor);
    title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0d(mod),3) ' Hz'])
    grid on; 
end

% Plot MDd
for mod=1%:size(MDd,2)
    
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    figure('units','normalized','position',[.2 .1 .6 .8])
    
    PlotMesh(nodes_defected, elementPlot, 0);
    d1 = reshape(MDd(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes_defected, elementPlot, d1, 'factor', Ly*0.3 );
    %title(['MD:' MDs_names{mod}])
    title(['\theta_{' num2str(MDd_names(mod,1)) num2str(MDd_names(mod,2)) '}'])
    axis on; grid on; box on
end
%% NL full order model with defects (using tangenstiff_defected) 
q0 = BeamAssembly.constrain_vector(u0);
qd0 = BeamAssembly.constrain_vector(v0);
qdd0 = BeamAssembly.constrain_vector(a0);
TI_NL_def = ImplicitNewmark('timestep',h,'alpha',0.005);
BeamAssembly.DATA.Ud=U*xi;
% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear_defected(q,qd,qdd,t,BeamAssembly,F_ext);

% Nonlinear Time Integration

tic
TI_NL_def.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL_def.Solution.u = BeamAssembly.unconstrain_vector(TI_NL_def.Solution.q);
FullNLdefduration=toc
%% HFOM-d Linear 
% Initial condition: equilibrium
u00 = zeros(BeamDefectedAssembly.Mesh.nDOFs, 1);
v00 = zeros(BeamDefectedAssembly.Mesh.nDOFs, 1);
a00 = zeros(BeamDefectedAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0))

q0 = BeamDefectedAssembly.constrain_vector(u00);
qd0 = BeamDefectedAssembly.constrain_vector(v00);
qdd0 = BeamDefectedAssembly.constrain_vector(a00);

% Instantiate object for linear time integration
TI_lin_HFMd = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,BeamDefectedAssembly,F_ext);

tic
TI_lin_HFMd.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin_HFMd.Solution.ud = BeamDefectedAssembly.unconstrain_vector(TI_lin_HFMd.Solution.q);
FullLDurationHFMd=toc
%% HFOM-d NL
% Instantiate object for nonlinear time integration
TI_NL_HFMd = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamDefectedAssembly,F_ext);

% Nonlinear Time Integration

tic
TI_NL_HFMd.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL_HFMd.Solution.ud = BeamDefectedAssembly.unconstrain_vector(TI_NL_HFMd.Solution.q);
FullNLdurationHFMd=toc
%% ROM for defected case ROM-d using tensors
LROMdd=[];
Ltimedd=[];
Ldurationdd=[];
NLROMdd=[];
NLtimedd=[];
NLdurationdd=[];
ASEMdurationdd=[];
f0c_ROMdd={};

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode ); %Check if its the same node
dof=dof(2);

countdd=0;
VbModedd=[];
MDeridd=[];
for m=3 %m=[2 3 7 10]
    for  n=6 % n=[0 m*(m+1)/2]
        t=m+n;
        
        if n==0
            V = [VMd(:,1:m)];
        else
            V = [VMd(:,1:m) MDd(:,1:n)];
        end
        
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamDefectedAssembly.mass_matrix)*V(:,ii));
        end
        
        
        
        BeamReducedDefectedAssembly  = ReducedAssembly(MeshDefected,V);
        
        
        BeamReducedDefectedAssembly.DATA.M = V'*Md*V;
        BeamReducedDefectedAssembly.DATA.C = V'*BeamDefectedAssembly.DATA.C*V;
        %BeamReducedDefectedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
        
        
        
        %-------------------------------------------------------------------------
        %% using Julia
        tic
        
        tensors_ROMd = reduced_tensors_ROM(BeamReducedDefectedAssembly, elements, V);
        Q2=tensors_ROMd.Q2;
        Q3=tensors_ROMd.Q3;
        Q4=tensors_ROMd.Q4;
        Q3t=tensors_ROMd.Q3t;
        Q4t=tensors_ROMd.Q4t;
        durationddensor=toc;
        ASEMdurationdd=[ASEMdurationdd durationddensor];
        %compute Natural Freq
        f0_ROMd_i= sort(sqrt(eig(BeamReducedDefectedAssembly...
            .DATA.M\Q2))/2/pi) ;
        
        f0c_ROMdd=[f0c_ROMdd f0_ROMd_i];
        
        BeamReducedDefectedAssembly.DATA.K=Q2; %used in residuals
        
        %% L sim
        
        
        % q0 = zeros(t,1);
        % qd0 = zeros(t,1);
        % qdd0 = zeros(t,1);
        q0=V'*u0;
        qd0=V'*v0;
        qdd0=V'*a0;
        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
        
        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear_tensor(q,qd,qdd,t, ...
            BeamReducedDefectedAssembly,F_ext,Q2);
        
        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        Ldurationdd=[Ldurationdd durationL];
        Ltimedd=[Ltimedd TI_lin_red.Solution.time'];
        LROMdd=[LROMdd TI_lin_red.Solution.u(dof,:)'];
        
        
        
        %% NL
        
        TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(...
            q,qd,qdd,t,BeamReducedDefectedAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
        %Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_
        %tensor_diffApp(q,qd,qdd,t,BeamAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t,Vr,Mrrr,Crrr);
        
        % time integration
        tic
        TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_newmark_red);
        TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
        durationNL=toc
        
        NLdurationdd=[NLdurationdd durationNL];
        
        NLtimedd=[NLtimedd TI_NL_newmark_red.Solution.time'];
        NLROMdd=[NLROMdd TI_NL_newmark_red.Solution.u(dof,:)'];
        
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
        % NLdurationdd=[NLdurationdd durationNL];
        %
        % NLtimedd=[NLtimedd TI_NL_alpha_red.Solution.time'];
        % NLROMdd=[NLROMdd TI_NL_alpha_red.Solution.u(dof,:)'];
        %%
        countdd=countdd+1;
        VbModedd=[VbModedd m];
        MDeridd=[MDeridd n];
        
        
    end
end

[max_sizeT, max_indexT] = max(cellfun('size', f0c_ROMdd, 1));
f0_ROMdd=zeros(max_sizeT,size(f0c_ROMdd,2));

for ii=1:max_indexT
    f0c=f0c_ROMdd{ii};
    s=length(f0c);
    f0_ROMdd(1:s,ii)=f0c;
end
%% DPROM init and Defect Sensitivities (choose the defect, 1st VM)
FORMULATION='N1' % N1/N1t/N0
VOLUME=1 ;        % defected volume =1 . nominal volume=0

% U = V(:,1:2);       % defect basis
% xi = rand(size(U,2),1)*0;
ndef=size(U,2); %number of defects
[DS, DSnames] = defect_sensitivities(BeamAssembly, elements, VMs, U, ...
    FORMULATION); % for DpROM

%%  DpROM VM (name of variable **DpROM) 1st VM
LDPROM=[];
LtimeDP=[];
LdurationDP=[];
NLDPROM=[];
NLtimeDP=[];
NLdurationDP=[];
ASEMdurationDP=[];
f0c_DPROM={};
nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

countDP=0;
VbModeDp=[];
MDeriDp=[];
DSenti=[];
for m=3 %=[2 3 7 10]
    for  n=6 %n=[0 m*(m+1)/2]
        for k=[ndef*m] % always include DS or so
            
            if n==0 && k~=0
                break
            end
            
            t=m+n+k;
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
                V = [VMs(:,1:m) MDs(:,1:n) DS(:,1:k)];
            end
            
            %mass normalization
            for ii = 1 : size(V, 2)
                V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamAssembly.mass_matrix)*V(:,ii));
            end
            
            %V=BeamAssembly.unconstrain_vector(V);
            
            BeamDefectedReducedAssembly  = ReducedAssembly(myMesh,V);
            
            
            % Vr=V(BeamDefectedReducedAssembly.Mesh.EBC.unconstrainedDOFs(:)',:);
            % Vr=V;
            
            %
            BeamDefectedReducedAssembly.DATA.M = BeamDefectedReducedAssembly.mass_matrix();
            BeamDefectedReducedAssembly.DATA.C = BeamDefectedReducedAssembly.damping_matrix(0,0,u0);
            BeamDefectedReducedAssembly.DATA.K =  BeamDefectedReducedAssembly.stiffness_matrix(u0);
            
            
            
            %% --------use of T in yetAnotherFEcode-------------------------------------
            % Q2=BeamDefectedReducedAssembly.DATA.K;
            %  % Q2=Krrr;
            % Q3r=BeamDefectedReducedAssembly.tensor('tensor_Q3',[n+m n+m n+m],[]);
            % Q4r=BeamDefectedReducedAssembly.tensor('tensor_Q4',[n+m n+m n+m n+m],[]);
            % % Q3 = BeamDefectedReducedAssembly.constrain_tensor(Q3r);
            % % Q4 =BeamDefectedReducedAssembly.constrain_tensor(Q4r);
            % Q3=full(Q3r);
            % Q4=full(Q4r);
            %-------------------------------------------------------------------------
            %% using Julia
            tic
            
            tensors_DpROM = reduced_tensors_DpROM(BeamAssembly, elements, V,U ...
                ,FORMULATION,VOLUME);
            % evaluate the defected tensors at xi
            [Q2, Q3, Q4, Q3t, Q4t] = DefectedTensors(tensors_DpROM, xi);
            
            durationDPensor=toc;
            ASEMdurationDP=[ASEMdurationDP durationDPensor];
            %compute Natural Freq
            f0_ROM_i= sort(sqrt(eig(BeamDefectedReducedAssembly...
                .DATA.M\Q2))/2/pi) ;
            
            f0c_DPROM=[f0c_DPROM f0_ROM_i];
            
            
            %% % parametric formulation for defects (only to test)---------------------
            
            % FORMULATION = 'N1'; % N1/N1t/N0
            % VOLUME = 1;         % integration over defected (1) or nominal volume (0)
            % U = V(:,1:2);       % defect basis
            % tensors_DpROM = reduced_tensors_DpROM(FORMULATION, VOLUME, ...
            %     BeamAssembly, elements, V, U);
            %
            % % evaluate the defected tensors at xi
            % xi = rand(size(U,2),1);
            % [Q2, Q3, Q4] = DefectedTensors(tensors_DpROM, xi);
            %% this also works (moving from full to reduce tensor) using yet antoher FE
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
            %% L sim
            
            % Q3t = Q3 + permute(Q3,[1 3 2]);
            % Q4t = Q4 + permute(Q4,[1 3 2 4]) + permute(Q4,[1 4 2 3]);
            
            % q0 = zeros(t,1);
            % qd0 = zeros(t,1);
            % qdd0 = zeros(t,1);
            q0=V'*u0;
            qd0=V'*v0;
            qdd0=V'*a0;
            TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
            
            % Modal linear Residual evaluation function handle
            Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear_tensor(q,qd,qdd,t,BeamDefectedReducedAssembly,F_ext,Q2);
            
            % time integration
            tic
            TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
            TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
            durationL=toc
            LdurationDP=[LdurationDP durationL];
            LtimeDP=[LtimeDP TI_lin_red.Solution.time'];
            LDPROM=[LDPROM TI_lin_red.Solution.u(dof,:)'];
            
            
            
            %% NL
            
            TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
            % Modal nonlinear Residual evaluation function handle
            Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,BeamDefectedReducedAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
            %Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor_diffApp(q,qd,qdd,t,BeamAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t,Vr,Mrrr,Crrr);
            
            % time integration
            tic
            TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_newmark_red);
            TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
            durationNL=toc
            
            NLdurationDP=[NLdurationDP durationNL];
            
            NLtimeDP=[NLtimeDP TI_NL_newmark_red.Solution.time'];
            NLDPROM=[NLDPROM TI_NL_newmark_red.Solution.u(dof,:)'];
            
            %% TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
            %
            %
            % TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
            % % Modal nonlinear Residual evaluation function handle
            % Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,BeamDefectedReducedAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
            %
            % % time integration
            % tic
            % TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
            % TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
            % durationNL=toc
            %
            % NLdurationDP=[NLdurationDP durationNL];
            %
            % NLtimeDP=[NLtimeDP TI_NL_alpha_red.Solution.time'];
            % NLDPROM=[NLDPROM TI_NL_alpha_red.Solution.u(dof,:)'];
            %%
            countDP=countDP+1;
            VbModeDp=[VbModeDp m];
            MDeriDp=[MDeriDp n];
            DSenti=[DSenti k];
            
        end
    end
end
[max_sizeT, max_indexT] = max(cellfun('size', f0c_DPROM, 1));
f0_DpROM=zeros(max_sizeT,size(f0c_DPROM,2));

for ii=1:max_indexT
    f0c=f0c_DPROM{ii};
    s=length(f0c);
    f0_DpROM(1:s,ii)=f0c;
end
%% HyperReduced Solution defected
LROMHd=[];
LtimeHd=[];
LdurationHd=[];
NLROMHd=[];
NLtimeHd=[];
NLdurationHd=[];
f0c_ROMHd={};

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

countHd=0;
VbModeHd=[];
MDeriHd=[];
DSentiHd=[];
for m=3;% m=[2 3 7 10]
    for n=m*(m+1)/2
     for k=m
        
        t=m+n+k;
      
%%
        
        if n==0
            V = [VMs(:,1:m)];
        else
            V = [VMs(:,1:m) MDs(:,1:n) DS(:,1:k)];
        end
        %V=[ MDs(:,1:6)];
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamAssembly.mass_matrix)*V(:,ii));
        end
        BeamReducedAssemblyDefected  = ReducedAssembly(myMesh,V);
        
        
        BeamReducedAssemblyDefected.DATA.M = BeamReducedAssemblyDefected.mass_matrix();
        BeamReducedAssemblyDefected.DATA.C = BeamReducedAssemblyDefected.damping_matrix(0,0,u0);
        BeamReducedAssemblyDefected.DATA.K =  BeamReducedAssemblyDefected.stiffness_matrix(u0);
        BeamReducedAssemblyDefected.DATA.Ud=U*xi;

        %% algorithm 3d
        V_H=V(:,1:m)
        Lin_sol=TI_lin.Solution.u;
        Lin_sol_snap=[];
        for ii=1:size(Lin_sol,2)
            if rem(ii,10)==0
                Lin_sol_snap=[Lin_sol_snap Lin_sol(:,ii)];
            end
        end
        eta=V_H'*Lin_sol_snap;
        
        Theta = QM_Theta_from_SMDs(V(:,m+1:m+n), MDs_names(1:n,:));
        Xi = DS_Xi_QM(V(:,m+n+1:m+n+k), DSnames);
        %QM uplifting
        uu = zeros(size(V_H,1), size(eta, 2));
        for tt = 1 : size(eta, 2)
            uu(:,tt) = einsum('Iij,iJ,jK->IJK', Theta, eta(:,tt), eta(:,tt));
        end
        %
                uud = zeros(size(V_H,1), size(eta, 2));
        for tt = 1 : size(eta, 2)
            for ii= 1:size(xi,2)
                uud(:,tt) =uud(:,tt)+ einsum('Iij,iJ,jK->IJK', Xi, eta(:,tt), xi(:,ii));
            end
        end
        
        u_lin_ECSW = V_H * eta + 1/2*uu+uud;
     
        %Construct Gb
        qq=(V'*V)^-1*V'*u_lin_ECSW;
        tic
        [G,b]=BeamReducedAssemblyDefected.constructGb(qq);
        GBconstructTime2=toc;
        %fNNLS
        tic
        [x_fnnls,w_fnnls]=fnnls(G'*G,G'*b,(norm(b)*0.01));
        fnnlsTime=toc
        tic;[x_nnls,w_nnls,info]=nnls(full(G),full(b),struct('Accy',1,'Tol',1/((norm(b)*5))));
        nnlTime=toc
%         options = optimset('TolX',(norm(b)*0.01));
%         tic;x_lsq=lsqnonneg(G,b,options);
        lsqTime=toc
        nnz(x_fnnls)
        nnz(x_nnls)
%         nnz(x_lsq)
 x_sNNLS=sNNLS(G,b,0.0035);
         nnz(x_sNNLS)
        BeamReducedAssemblyDefected.DATA.elementWeights=x_sNNLS;
        %%
        f0_ROM_i= sort(sqrt(eig(BeamReducedAssemblyDefected.DATA.M\BeamReducedAssemblyDefected...
            .DATA.K))/2/pi) ;
        f0c_ROMHd=[f0c_ROMHd f0_ROM_i];
        %
        q0 = zeros(t,1);
        qd0 = zeros(t,1);
        qdd0 = zeros(t,1);
        
        
        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
        
        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,BeamReducedAssemblyDefected,F_ext);
        
        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        LdurationHd=[LdurationHd durationL];
        LtimeHd=[LtimeHd TI_lin_red.Solution.time'];
        LROMHd=[LROMHd TI_lin_red.Solution.u(dof,:)'];
        
        %% Reduced solution Noninear
        % For demonstration purposes, we simply reduce the nonlinear system using
        % out-of-plane bending modes. This is expected to produce bad results when
        % in-plane stretching is involved in the response.
        
        
        %TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
        TI_NL_alpha_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hyper_defected(q,qd,qdd,t,BeamReducedAssemblyDefected,F_ext);
        
        % time integration
        tic
        TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
        durationNL=toc
        
        NLdurationHd=[NLdurationHd durationNL];
        
        NLtimeHd=[NLtimeHd TI_NL_alpha_red.Solution.time'];
        NLROMHd=[NLROMHd TI_NL_alpha_red.Solution.u(dof,:)'];
        
        countHd=countHd+1;
        VbModeHd=[VbModeHd m];
        MDeriHd=[MDeriHd n];
          DSentiHd=[DSentiHd k];
    end
    end
end
[max_size, max_index] = max(cellfun('size', f0c_ROMHd, 1));
f0_ROM=zeros(max_size,size(f0c_ROMHd,2));
for ii=1:max_index
    f0c=f0c_ROMHd{ii};
    s=length(f0c);
    f0_ROM(1:s,ii)=f0c;
end
%% Plotting
nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

figure;
hold on

%plot full simulation
plot(TI_lin.Solution.time, TI_lin.Solution.u(dof,:),'DisplayName', ['Full linear (Newmark) ' 'time= ' num2str(FullLDuration) ])
hold on
plot(TI_NL.Solution.time, TI_NL.Solution.u(dof,:),'DisplayName', ['Full nonlinear (Newmark) ' 'time= ' num2str(FullNLduration)])
hold on
%plot full simulation for defected case
plot(TI_lin_HFMd.Solution.time, TI_lin_HFMd.Solution.ud(dof,:),'DisplayName', ['Full linear Defected(Newmark) ' 'time= ' num2str(FullLDurationHFMd) ])
hold on
plot(TI_NL_HFMd.Solution.time, TI_NL_HFMd.Solution.ud(dof,:),'DisplayName', ['Full nonlinear Defected(Newmark) ' 'time= ' num2str(FullNLdurationHFMd)])
hold on
%plot FOM with tangenstiff def
plot(TI_NL_def.Solution.time, TI_NL_def.Solution.u(dof,:),'DisplayName', ['Full nonlinear defKF ' 'time= ' num2str(FullNLdefduration)])
hold on
%plot normal ROM
for i=1:count
    Vmode=VbMode(i);
    Mderi=MDeri(i);
    plot(Ltime(:,i), LROM(:,i),'DisplayName', ['Reduced linear (Newmark) ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(Lduration(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtime(:,i), NLROM(:,i),'DisplayName', ['Reduced nonlinear (Newmark) ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(NLduration(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot ROM-Hyper
for i=1:countH
    VmodeH=VbModeH(i);
    MderiH=MDeriH(i);
    plot(LtimeH(:,i), LROMH(:,i),'DisplayName', ['Hyperreduced linear (Newmark) ' ' VM=' num2str(VmodeH) ' MD=' num2str(MderiH) ' time=' num2str(LdurationH(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeH(:,i), NLROMH(:,i),'DisplayName', ['Hyperreduced nonlinear (Newmark) ' ' VM=' num2str(VmodeH) ' MD=' num2str(MderiH) ' time=' num2str(NLdurationH(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot Tensor normal ROM
for i=1:countT
    VmodeT=VbModeT(i);
    MderiT=MDeriT(i);
    ASEMdurationT_i=ASEMdurationT(i);
    NLdurationT_i=NLdurationT(i);
    totalduration_i=ASEMdurationT_i+NLdurationT_i;  % time with assem time of tensors
    plot(LtimeT(:,i), LROMT(:,i),'DisplayName', ['Reduced linear Tensor(Newmark)  ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(LdurationT(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeT(:,i), NLROMT(:,i),'DisplayName', ['Reduced nonlinear (Newmark) TensorApproach ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(totalduration_i) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot Tensor  DpROM
for i=1:countDP
    VmodeDP=VbModeDp(i);
    MderiDP=MDeriDp(i);
    DSsenti=DSenti(i)
    ASEMdurationT_i=ASEMdurationDP(i);
    NLdurationT_i=NLdurationDP(i);
    totalduration_i=ASEMdurationT_i+NLdurationT_i;  % time with assem time of tensors
    plot(LtimeDP(:,i), LDPROM(:,i),'DisplayName', ['DPROM L (Newmark)  ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(LdurationDP(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeDP(:,i), NLDPROM(:,i),'DisplayName', ['DPROM NL (Newmark)  ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(totalduration_i) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot Tensor defected ROM
for i=1:countdd
    Vmoded=VbModedd(i);
    Mderid=MDeridd(i);
    ASEMdurationd_i=ASEMdurationdd(i);
    NLdurationd_i=NLdurationdd(i);
    totalduration_i=ASEMdurationd_i+NLdurationd_i;  % time with assem time of tensors
    plot(Ltimedd(:,i), LROMdd(:,i),'DisplayName', ['Reduced linear defected Mesh (Newmark)  ' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(Ldurationdd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimedd(:,i), NLROMdd(:,i),'DisplayName', ['Reduced nonlinear defected Mesh(Newmark) TensorApproach ' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(totalduration_i) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot defected Hyper
for i=1:countHd
    VmodeHd=VbModeHd(i);
    MderiHd=MDeriHd(i);
    DSsentiHd=DSentiHd(i)
      plot(LtimeHd(:,i), LROMHd(:,i),'DisplayName', ['hyperDPROM L (Newmark)  ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' DS=' num2str(DSsentiHd) ' time=' num2str(LdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeHd(:,i), NLROMHd(:,i),'DisplayName', ['HyperDPROM NL (Newmark)  ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' DS=' num2str(MderiHd) ' time=' num2str(NLdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot Tensor normal ROM