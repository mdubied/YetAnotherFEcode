% FreqDIV
clear;
close all;
clc

whichModel = 'ANSYS'; % or "ABAQUS" or %CUSTOM


%% PREPARE MODEL

% DATA ____________________________________________________________________
E       = 170e9;     % Young's modulus [Pa]
rho     = 2329;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio
thickness = 1e-5;     % [m] beam's out-of-plane thickness
ngauss=2;
% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial,ngauss);

% MESH_____________________________________________________________________

switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'Job-BeamQuad';
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
    case 'ANSYS'
        % Alternatively, one can write an input file in ANSYS and read it as:
        
        meshfile='QGmesh.dat';
        meshtype='QUAD8';
        model=Ansys2Matlab_mesh(meshfile,meshtype);
        model.nodes=model.nodes(:,2:3);          % 2:3 because 2D
        nodes=model.nodes;
        model.elements=model.elements(:,2:end);
        elements=model.elements;
end

myMesh = Mesh(nodes);
myMesh.create_elements_table(model.elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
switch upper( whichModel )
    case 'CUSTOM'
        myMesh.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
    case 'ABAQUS'
        myMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
    case 'ANSYS'
        fixed_nodes = model.fixedNodes(:)';
        myMesh.set_essential_boundary_condition((fixed_nodes),1:2,0);
end

% plot
elementPlot =elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
hold on
plot(nodes(fixed_nodes,1),nodes(fixed_nodes,2),'r*')
%plot3(nodes(fixed_nodes,1),nodes(fixed_nodes,2),nodes(fixed_nodes,3),'r*')
hold off
%% ASSEMBLY ________________________________________________________________
FreqDivAssembly = Assembly(myMesh);
M = FreqDivAssembly.mass_matrix();
nNodes = size(nodes,1);
u0_1 = zeros( myMesh.nDOFs, 1);
[K,~] = FreqDivAssembly.tangent_stiffness_and_force(u0_1);

C = 761*M +1.96e-10*K;

alfa=761;
beta=1.96e-10;
%  beta=0;
% C=0*K+0*M;

%C2=FreqDivAssembly.damping_matrix(alfa,beta,u0_1);
% store matrices
FreqDivAssembly.DATA.K = K;
FreqDivAssembly.DATA.M = M
FreqDivAssembly.DATA.C= C;

%% old use of VMs
%  Kc = FreqDivAssembly.constrain_matrix(K);
%  Mc = FreqDivAssembly.constrain_matrix(M);
%
% Eigenvalue problem_______________________________________________________
% n_VMs = 8; % first n_VMs modes with lowest frequency calculated

% [V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
% [f0,ind] = sort(sqrt(diag(om))/2/pi);
% V0 = V0(:,ind);
% for ii = 1:n_VMs
%     V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
% end
% V0 = FreqDivAssembly.unconstrain_vector(V0);
%or this below


% [VMs,f0,time_vm] = VMs_computes(Kc,Mc,n_VMs,1);
% for ii = 1:n_VMs
%     VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
% end
% VMs = FreqDivAssembly.unconstrain_vector(VMs);
%% VMs
Kc = FreqDivAssembly.constrain_matrix(K);
Mc = FreqDivAssembly.constrain_matrix(M);

n_VMs = 15;
[VMs,f0,time_vm]=FreqDivAssembly.VMs_compute(n_VMs,1);
figure('units','normalized','position',[.2 .1 .6 .8])
hold on
%normalization dispalcement
for ii = 1:size(VMs,2)
    VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
end
% for ii = 1:size(A,2)
%     VMsA(:,ii) = A(:,ii)/max(sqrt(sum(A(:,ii).^2,2)));
% end
% PNx=ceil(n_VMs/2); %Plot number
% PNy=n_VMs-PNx;
% if PNx==1
%     PNx=2;
%     PNy=1;
% elseif PNx==2
%     PNx=3;
%     PNy=1;
% end

% PLOT
for mod=1:5
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    figure('units','normalized','position',[.2 .1 .6 .8])
    
    %subplot(4,4,mod)
    
    PlotMesh(nodes, elementPlot, 0);
    v1 = reshape(VMs(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes, elementPlot, -v1, 'factor', 1e-5);
    title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
    
end

%% MDs
NMDs = 120;
% [MDs,MDs_names,time_md] = FreqDivAssembly.MDs_compute( n_VMs,NMDs,u0_1);
%using julia
[MDs,MDs_names]= modal_derivatives(FreqDivAssembly,elements,VMs);

%[MDss,MDss_names]= modal_derivatives(FreqDivAssembly,elements,VMs);
%using julia
figure('units','normalized','position',[.2 .1 .6 .8])
hold on
% PNxx=ceil(size(MDs,2)/2); %Plot number
% PNyy=size(MDs,2)-PNxx;
% if PNxx==1
%     PNxx=2;
%     PNyy=1;
% elseif PNxx==2
%     PNxx=3;
%     PNyy=1;
% end
%normalization
for ii = 1:size(MDs,2)
    MDs(:,ii) = MDs(:,ii)/max(sqrt(sum(MDs(:,ii).^2,2)));
    %MDss(:,ii)=MDss(:,ii)/max(sqrt(sum(MDss(:,ii).^2,2)));
end
%plot
for mod=[1] %:size(MDs,2)
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    %figure('units','normalized','position',[.2 .1 .6 .8])
    if mod<5
        subplot(4,4,mod)
    elseif mod>5 && mod<30
        subplot(4,4,mod-10)
    elseif mod>29 && mod<32
        subplot(4,4,mod-19)
    elseif mod>32
        subplot(4,4,mod-27)
    end
    
    %PlotMesh(nodes, elementPlot, 0);
    d1 = reshape(MDs(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes, elementPlot, d1, 'factor', 1e-5);
    title(['\theta_{' num2str(MDs_names(mod,1)) num2str(MDs_names(mod,2)) '}'])
    
end
close all
%% Force

% Define external force:
% % Body force
% Pressure = 1e6;
% F = Pressure*FreqDivAssembly.uniform_body_force();

% Nodal force
F = zeros(myMesh.nDOFs,1);
%nf = find_node(9.25e-07,0.000183,[],nodes); % node where to put the force
nf=2564;
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(1)) = 135e-6 % -135e-6;
%% Static Analysis
u_lin = FreqDivAssembly.solve_system(K, F);
ULIN = reshape(u_lin,2,[]).';	% Linear response
u = static_equilibrium(FreqDivAssembly, F, 'display', 'iter-detailed');
UNL = reshape(u,2,[]).';        % Nonlinear response

fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ...
    '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 1500;
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,UNL,'factor',scale,'color','k');
colormap jet
title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])

%% Dynamic response Initlialization
% forcing frequency of the average of first two natural frequencies
omega_ext = 2*pi*868000;
%omega_ext = pi*2*1.02e5;
T =  2*pi/omega_ext; % time period of forcing

% load amplification factor
amplification_factor = 1;

% forcing function

F_ext = @(t) amplification_factor * F * sin(omega_ext * t);

% Initial condition: equilibrium
u0 = zeros(FreqDivAssembly.Mesh.nDOFs, 1);
v0 = zeros(FreqDivAssembly.Mesh.nDOFs, 1);
a0 = zeros(FreqDivAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0))

q0 = FreqDivAssembly.constrain_vector(u0);
qd0 = FreqDivAssembly.constrain_vector(v0);
qdd0 = FreqDivAssembly.constrain_vector(a0);

% time step for integration
h = T/150;

% Precompute data for Assembly object
% FreqDivAssembly.DATA.M = M;
% FreqDivAssembly.DATA.K = K;
% FreqDivAssembly.DATA.C = C; %rayleigh

tmax = 700*T;
%% DOF used in Ploting and Saving
%nfsense1 = find_node(9.25e-07,91.5e-6,[],nodes); % node where to see the results
nfsense1=2754;
dof1 = get_index(nfsense1, myMesh.nDOFPerNode );
dof1=dof1(1);
%nfsense2 = find_node(148.85e-6,182.07e-6,[],nodes); % node where to see the results
nfsense2=2144;
dof2 = get_index(nfsense2, myMesh.nDOFPerNode );
dof2=dof2(2);
%nfsense3 = find_node(16.775e-6,22.85e-6,[],nodes); % node where to see the results
nfsense3=1483;
dof3 = get_index(nfsense3, myMesh.nDOFPerNode );
dof3=dof3(1);
%% Linear Dynamic Response
% Instantiate object for linear time integration

TI_lin = ImplicitNewmark('timestep',T/50,'alpha',0.005,'linear',true); %h is different fast

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,FreqDivAssembly,F_ext);

tic
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution

TI_lin.Solution.u = FreqDivAssembly.unconstrain_vector(TI_lin.Solution.q);

% lin_SOLUTION=TI_lin.Solution;

FullLDuration=toc

TI_lin_time= TI_lin.Solution.time';
TI_lin_dof1= TI_lin.Solution.u(dof1,:)';
TI_lin_dof2= TI_lin.Solution.u(dof2,:)';
TI_lin_dof3= TI_lin.Solution.u(dof3,:)';

save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\full_solution\TI_lin', ...
    'FullLDuration','TI_lin_time','TI_lin_dof1','TI_lin_dof2','TI_lin_dof3','-v7.3');

% Animate solution on Mesh (very slow)
%AnimateFieldonDeformedMesh(myMesh.nodes,myMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:2,'filename','lineardisp')
%% Non Linear Dyanmic Response
tmax = T;
% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,FreqDivAssembly,F_ext);

% Nonlinear Time Integration

tic
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = FreqDivAssembly.unconstrain_vector(TI_NL.Solution.q);
FullNLduration=toc
% NL_SOLUTION=TI_NL.Solution;



TI_NL_time= TI_NL.Solution.time';
TI_NL_dof1= TI_NL.Solution.u(dof1,:)';
TI_NL_dof2= TI_NL.Solution.u(dof2,:)';
TI_NL_dof3= TI_NL.Solution.u(dof3,:)';

save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\full_solution\TI_NL', ...
    'FullNLduration','TI_NL_time','TI_NL_dof1','TI_NL_dof2','TI_NL_dof3','-v7.3');
%  save('TI_NL_FREQ_test','-struct','NL_SOLUTION','-v7.3');

% Generalized alpha scheme
% linear
% TI_lin_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7, 'linear',true);
% TI_lin_alpha.Integrate(q0,qd0,qdd0,tmax,residual_lin);
% TI_lin_alpha.Solution.u = FreqDivAssembly.unconstrain_vector(TI_lin_alpha.Solution.q);
%
% % nonlinear
% TI_NL_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
% TI_NL_alpha.Integrate(q0,qd0,qdd0,tmax,residual);
% TI_NL_alpha.Solution.u = FreqDivAssembly.unconstrain_vector(TI_NL_alpha.Solution.q);
%% Maixmum Modal interaction  MMI for normal case

%Linear Solution is Needed to Get the Best MDs
%Sorts the MDs index based on highest to lowest modal interaction (time
% depedent)
W=[];
W_names=[];
for ii = 1 : size(VMs, 2)
    V_MMI(:,ii) = VMs(:,ii) / (VMs(:,ii)'*(FreqDivAssembly.DATA.M)*VMs(:,ii));
end
%
eta=V_MMI'*TI_lin.Solution.u;
for ii=1:size(eta,1)
    for jj=ii:size(eta,1)
        W_names=[ii ;jj];
        MMI= max(eta(ii,:).*eta(jj,:));
        W=[W [W_names ; MMI]];
    end
end

figure()

hold on

scatterbar3(W(1,:),W(2,:),normalize(W(3,:),'range'),1)
scatter3(W(1,:),W(2,:),normalize(W(3,:),'range'),'*')
grid on ;box on;
xticks(unique(W(1,:)))
yticks(unique(W(2,:)))

sortedW=sortrows(W',3,'descend')';

sortedMDs=[];
for ii=1:size(sortedW,2)
    I_index=sortedW(1,ii);
    J_index=sortedW(2,ii);
    MD_ii=find(and(MDs_names(:,1)== I_index,MDs_names(:,2)== J_index));
    if isempty(MD_ii)
        MD_ii=find(and(MDs_names(:,1)== J_index,MDs_names(:,2)== I_index));
    end
    sortedMDs=[sortedMDs MD_ii];
end
save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\full_solution\sortedMDs', ...
    'sortedMDs','-v7.3');
%% Reduced solution Linear/NL
% tmax = 10*T;
% LROM1=[];
% LROM2=[];
% LROM3=[];
%
% Ltime=[];
% Lduration=[];
%
% NLROM1=[];
% NLROM2=[];
% NLROM3=[];
%
% NLtime=[];
% NLduration=[];
% f0c_ROM={};
%
% count=0;
% VbMode=[];
% MDeri=[];
% VMss=orth(VMs);
% MDss=orth(MDs);
% for m=[13]
%     for n=[45]% m*(m+1)/2]
%
%         t=m+n;
%
%         if n==0
%             V = [VMs(:,1:m)];
%         else
%             V = [VMs(:,1:m) MDs(:,1:n)];
%         end
%         %mass normalization
%         for ii = 1 : size(V, 2)
%             V(:,ii) = V(:,ii) / (V(:,ii)'*(FreqDivAssembly.DATA.M)*V(:,ii));
%         end
%
%         FreqDivReducedAssembly  = ReducedAssembly(myMesh,V);
%
%
%         %FreqDivReducedAssembly.DATA.M = FreqDivReducedAssembly.mass_matrix();
%         FreqDivReducedAssembly.DATA.M =V'*FreqDivAssembly.DATA.M*V;
%         %FreqDivReducedAssembly.DATA.K =  FreqDivReducedAssembly.stiffness_matrix(u0);
%         FreqDivReducedAssembly.DATA.K =V'*FreqDivAssembly.DATA.K*V;
%
%         FreqDivReducedAssembly.DATA.C = alfa*FreqDivReducedAssembly.DATA.M+beta*FreqDivReducedAssembly.DATA.K;
%         %calculate NF for reduced case
%         f0_ROM_i= sort(sqrt(eig(FreqDivReducedAssembly.DATA.M\FreqDivReducedAssembly...
%             .DATA.K))/2/pi) ;
%         f0c_ROM=[f0c_ROM f0_ROM_i];
%         %
%         %         q0 = zeros(t,1);
%         %         qd0 = zeros(t,1);
%         %         qdd0 = zeros(t,1);
%
%         q0=V'*u0;
%         qd0=V'*v0;
%         qdd0=V'*a0;
%
%         TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
%
%         % Modal linear Residual evaluation function handle
%         Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,FreqDivReducedAssembly,F_ext);
%
%         % time integration
%         tic
%         TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
%         TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
%         durationL=toc
%
%         Lduration=[Lduration durationL];
%         Ltime=[Ltime TI_lin_red.Solution.time'];
%         LROM1=[LROM1 TI_lin_red.Solution.u(dof1,:)'];
%         LROM2=[LROM2 TI_lin_red.Solution.u(dof2,:)'];
%         LROM3=[LROM3 TI_lin_red.Solution.u(dof3,:)'];
%
%         %% Reduced solution Noninear
%
%         TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
%
%         % Modal nonlinear Residual evaluation function handle
%         Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,FreqDivReducedAssembly,F_ext);
%
%         % time integration
%         tic
%         TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
%         TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
%         durationNL=toc
%
%         NLduration=[NLduration durationNL];
%
%         NLtime=[NLtime TI_NL_newmark_red.Solution.time'];
%         NLROM1=[NLROM1 TI_NL_newmark_red.Solution.u(dof1,:)'];
%         NLROM2=[NLROM2 TI_NL_newmark_red.Solution.u(dof2,:)'];
%         NLROM3=[NLROM3 TI_NL_newmark_red.Solution.u(dof3,:)'];
%         %%
%         count=count+1;
%         VbMode=[VbMode m];
%         MDeri=[MDeri n];
%     end
% end
% [max_size, max_index] = max(cellfun('size', f0c_ROM, 1));
% f0_ROM=zeros(max_size,size(f0c_ROM,2));
% for ii=1:max_index
%     f0c=f0c_ROM{ii};
%     s=length(f0c);
%     f0_ROM(1:s,ii)=f0c;
% end
%% HyperReduced solution Linear/NL
tmax = 400*T;
LROM1H=[];
LROM2H=[];
LROM3H=[];

LtimeH=[];
LdurationH=[];
NLROM1H=[];
NLROM2H=[];
NLROM3H=[];

NLtimeH=[];
NLdurationH=[];
f0c_ROMH={};

countH=0;
VbModeH=[];
MDeriH=[];
VMssH=orth(VMs);
MDssH=orth(MDs);
for m=[15]
    for n=[45]% m*(m+1)/2]

        t=m+n;

        if n==0
            V = [VMs(:,1:m)];
        else
            V = [VMs(:,1:m)  MDs(:,[sortedMDs(1:n)]) ]; %MDs(:,1:n)]; %
        end
        %mass normalization
        V=orth(V);
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(FreqDivAssembly.DATA.M)*V(:,ii));
        end

        FreqDivReducedAssembly  = ReducedAssembly(myMesh,V);


        %FreqDivReducedAssembly.DATA.M = FreqDivReducedAssembly.mass_matrix();
        FreqDivReducedAssembly.DATA.M =V'*FreqDivAssembly.DATA.M*V;
        %FreqDivReducedAssembly.DATA.K =  FreqDivReducedAssembly.stiffness_matrix(u0);
        FreqDivReducedAssembly.DATA.K =V'*FreqDivAssembly.DATA.K*V;

        FreqDivReducedAssembly.DATA.C = alfa*FreqDivReducedAssembly.DATA.M+beta*FreqDivReducedAssembly.DATA.K;
        %% algorithm 3
         V_H=V(:,1:m);
        Lin_sol=TI_lin.Solution.u;
        Lin_sol_snap=[];
        for ii=1:size(Lin_sol,2)
            if rem(ii,350)==0
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
        [G,b]=FreqDivReducedAssembly.constructGb(qq);
        GBconstructTime=toc;
        
         %%fNNLS
%         tic
%         [x_fnnls,w_fnnls]=fnnls(G'*G,G'*b,norm(b)*0.01);
%         fnnlsTime=toc
%         nnz(x_fnnls)
        %nnls
%         tic;[x_nnls,w_nnls,info]=nnls(full(G),full(b),struct('Accy',1,'Tol',1/((norm(b)*10))));
%         nnlTime=toc
%         nnz(x_nnls)
%         %lsq
%         options = optimset('TolX',1/(norm(b)*0.01));
%         tic;x_lsq=lsqnonneg(G,b,options);
%         lsqTime=toc
%         nnz(x_lsq)
%         %mostafa snlls
%          x_sNNLS=sNNLS(G,b,0.0035);
%          nnz(x_sNNLS)
        %jain_snnls
        [E_jain, xi_jain]=snnls_j(full(G),full(b),25e-9)
        nnz(xi_jain)
        FreqDivReducedAssembly.DATA.elementWeights=xi_jain;
        %%
        %calculate NF for reduced case
        f0_ROM_i= sort(sqrt(eig(FreqDivReducedAssembly.DATA.M\FreqDivReducedAssembly...
            .DATA.K))/2/pi) ;
        f0c_ROMH=[f0c_ROMH f0_ROM_i];
        %
        q0 = zeros(t,1);
        qd0 = zeros(t,1);
        qdd0 = zeros(t,1);
        
        %         q0=V'*u0;
        %         qd0=V'*v0;
        %         qdd0=V'*a0;

        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,FreqDivReducedAssembly,F_ext);

        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc

        LdurationH=[LdurationH durationL];
        LtimeH=[LtimeH TI_lin_red.Solution.time'];
        LROM1H=[LROM1H TI_lin_red.Solution.u(dof1,:)'];
        LROM2H=[LROM2H TI_lin_red.Solution.u(dof2,:)'];
        LROM3H=[LROM3H TI_lin_red.Solution.u(dof3,:)'];

        %% Reduced solution Noninear

        TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);

        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hyper(q,qd,qdd,t,FreqDivReducedAssembly,F_ext);

        % time integration
        tic
        TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
        durationNL=toc

        NLdurationH=[NLdurationH durationNL];

        NLtimeH=[NLtimeH TI_NL_newmark_red.Solution.time'];
        NLROM1H=[NLROM1H TI_NL_newmark_red.Solution.u(dof1,:)'];
        NLROM2H=[NLROM2H TI_NL_newmark_red.Solution.u(dof2,:)'];
        NLROM3H=[NLROM3H TI_NL_newmark_red.Solution.u(dof3,:)'];
        %%
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
%% Redu   ced solution  Tensor Approach (normal ROM) (name of variable ***T)
tmax =700*T;
% h=T/200;
LROMT1=[];
LROMT2=[];
LROMT3=[];

LtimeT=[];
LdurationT=[];

NLROMT1=[];
NLROMT2=[];
NLROMT3=[];
NLtimeT=[];

NLdurationT=[];
ASEMdurationT=[];
f0c_ROMT={};

countT=0;
VbModeT=[];
MDeriT=[];
%   VMss=orth(VMs);
%  MDss=orth(MDs);
for m=[15]
    for  n=[55] %m*(m+1)/2]
        t=m+n;
        if n==0
            V = [VMs(:,1:m)];
        else
            %V = [VMs(:,1:m) MDs(:,1:n)];
            V = [VMs(:,1:m) MDs(:,[sortedMDs(1:n)])];
        end
        
        V=orth(V);
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(FreqDivAssembly.DATA.M)*V(:,ii));
        end
        
        
        FreqDivReducedAssemblyTensor  = ReducedAssembly(myMesh,V);
        
        
        %         FreqDivReducedAssembly.DATA.M = FreqDivReducedAssembly.mass_matrix();
        %         FreqDivReducedAssembly.DATA.C = FreqDivReducedAssembly.damping_matrix(0,0,u0);
        %         FreqDivReducedAssembly.DATA.K =  FreqDivReducedAssembly.stiffness_matrix(u0);
        
        FreqDivReducedAssemblyTensor.DATA.M =V'*FreqDivAssembly.DATA.M*V;
        FreqDivReducedAssemblyTensor.DATA.K =V'*FreqDivAssembly.DATA.K*V;
        
        FreqDivReducedAssemblyTensor.DATA.C = alfa*FreqDivReducedAssemblyTensor.DATA.M+beta*FreqDivReducedAssemblyTensor.DATA.K;
        
        %% using Julia
        if n==56
            FreqDivReducedAssemblyTensor.DATA.C = alfa*FreqDivReducedAssemblyTensor.DATA.M;
            C = 761*M +0*K;
            FreqDivAssembly.DATA.C= C;
        end
        
        tic
        
        tensors_ROM = reduced_tensors_ROM(FreqDivAssembly, elements, V);
        Q2=tensors_ROM.Q2;
        Q3=tensors_ROM.Q3;
        Q4=tensors_ROM.Q4;
        Q3t = Q3 + permute(Q3,[1 3 2]);
        Q4t = Q4 + permute(Q4,[1 3 2 4]) + permute(Q4,[1 4 2 3]);
        durationTensor=toc;
        ASEMdurationT=[ASEMdurationT durationTensor];
        
        %compute Natural Freq
        f0_ROM_i= sort(sqrt(eig(FreqDivReducedAssemblyTensor...
            .DATA.M\Q2))/2/pi) ;
        
        f0c_ROMT=[f0c_ROMT f0_ROM_i];
        
    save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROM\V15M55ROMTensors', ...
    'Q2','Q3','Q4','Q3t','Q4t','ASEMdurationT','-v7.3');
        %% L sim
        
        
        
        % q0 = zeros(t,1);
        % qd0 = zeros(t,1);
        % qdd0 = zeros(t,1);
        q0=V'*u0;
        qd0=V'*v0;
        qdd0=V'*a0;
        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
        
        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear_tensor(q,qd,qdd,t,FreqDivReducedAssemblyTensor,F_ext,Q2);
        
        % time integration
        tic
        
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
%         TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        
        LdurationT=[LdurationT durationL];
        LtimeT=[LtimeT TI_lin_red.Solution.time'];
        LROMT1=[LROMT1 (V(dof1,:)*TI_lin_red.Solution.q)'];
        LROMT2=[LROMT2 (V(dof2,:)*TI_lin_red.Solution.q)'];
        LROMT3=[LROMT3 (V(dof3,:)*TI_lin_red.Solution.q)'];
        
        
        save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROM\LROM', ...
    'LdurationT','LtimeT','LROMT1','LROMT2','LROMT3','-v7.3');
        %% NL
        
        TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,FreqDivReducedAssemblyTensor,F_ext,Q2,Q3,Q4,Q3t,Q4t);
        
        % time integration
        tic
        TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_newmark_red);
        %         TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
        durationNL=toc
        
        NLdurationT=[NLdurationT durationNL];
        
        NLtimeT=[NLtimeT TI_NL_newmark_red.Solution.time'];
        %         NLROMT1=[NLROMT1 TI_NL_newmark_red.Solution.u(dof1,:)'];
        %         NLROMT2=[NLROMT2 TI_NL_newmark_red.Solution.u(dof2,:)'];
        %         NLROMT3=[NLROMT3 TI_NL_newmark_red.Solution.u(dof3,:)'];
        
        NLROMT1=[NLROMT1 (V(dof1,:)*TI_NL_newmark_red.Solution.q)'];
        NLROMT2=[NLROMT2 (V(dof2,:)*TI_NL_newmark_red.Solution.q)'];
        NLROMT3=[NLROMT3 (V(dof3,:)*TI_NL_newmark_red.Solution.q)'];
        
        save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROM\NLROM', ...
    'NLdurationT','NLtimeT','NLROMT1','NLROMT2','NLROMT3','-v7.3');
        %%
        countT=countT+1;
        VbModeT=[VbModeT m];
        MDeriT=[MDeriT n];
        
        save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROM\ROMcount', ...
    'countT','VbModeT','MDeriT','-v7.3');
    end
end

[max_sizeT, max_indexT] = max(cellfun('size', f0c_ROMT, 1));
f0_ROMT=zeros(max_sizeT,size(f0c_ROMT,2));

for ii=1:max_indexT
    f0c=f0c_ROMT{ii};
    s=length(f0c);
    f0_ROMT(1:s,ii)=f0c;
end
 save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROM\f0_ROMT', ...
    'f0_ROMT','-v7.3');
1;
%% HFM-d init (Defected mesh) first VM
U=[VMs(:,2)];%[VMs(:,1) VMs(:,2) VMs(:,4)]     %defect Basis
xi=[2e-6];   %xi=[2e-6 2e-6 2e-6];
U1 = reshape(U(:,1), 2, []).'*xi(1);  %U1r = reshape(U(:,1), 2, []).'*xi; U=U1+U2+..
Ut=U1;
nodes_defected=nodes+ Ut;     %+U1r
% Mesh
MeshDefected=Mesh(nodes_defected);
MeshDefected.create_elements_table(elements,myElementConstructor);
%MeshDefected.set_essential_boundary_condition([nset{1} nset{3}],1:2,0);
MeshDefected.set_essential_boundary_condition((fixed_nodes),1:2,0)
% Assembly
FreqDivDefectedAssembly = Assembly(MeshDefected);
Md = FreqDivDefectedAssembly.mass_matrix();
[Kd,~] = FreqDivDefectedAssembly.tangent_stiffness_and_force(u0);
% store matrices
FreqDivDefectedAssembly.DATA.K = Kd;
FreqDivDefectedAssembly.DATA.M = Md;
FreqDivDefectedAssembly.DATA.C= alfa*Md+beta*Kd;
%plot the defected Mesh

elementPlot =elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes_defected, elementPlot, 0);
hold on
plot(nodes_defected(fixed_nodes,1),nodes_defected(fixed_nodes,2),'r*')


% Vm and Md computation for defected Mesh
n_Vmd=15;
% Vibration Modes (VM): defected------------------------------------------
Kdc = FreqDivDefectedAssembly.constrain_matrix(Kd);
Mdc = FreqDivDefectedAssembly.constrain_matrix(Md);
[VMd,om] = eigs(Kdc, Mdc, n_Vmd, 'SM');
[f0d,ind] = sort(sqrt(diag(om))/2/pi);
VMd = VMd(:,ind);
for ii = 1:n_Vmd
    VMd(:,ii) = VMd(:,ii)/max(sqrt(sum(VMd(:,ii).^2,2)));
end
VMd = FreqDivDefectedAssembly.unconstrain_vector(VMd);

% Modal Derivative Defected------------------------------------------------
[MDd, MDd_names] = modal_derivatives(FreqDivDefectedAssembly, elements, VMd);
for ii = 1:size(MDd,2)
    MDd(:,ii) = MDd(:,ii)/max(sqrt(sum(MDd(:,ii).^2,2)));
end
% ploting VMd------------------------------------------------------
for mod=3
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    %subplot(4,3,mod)
    figure('units','normalized','position',[.2 .1 .6 .8])
    
    PlotMesh(nodes_defected, elementPlot, 0);
    v1 = reshape(VMd(:,mod), 2, []).';
    factor = 0.1*max(nodes_defected(:,2));
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
    PlotFieldonDeformedMesh(nodes_defected, elementPlot, d1, 'factor', 1e-5 );
    %title(['MD:' MDs_names{mod}])
    title(['\theta_{' num2str(MDd_names(mod,1)) num2str(MDd_names(mod,2)) '}'])
    grid on; box on
end
%% HFOM-d Linear
tmax = 700*T;
% Initial condition: equilibrium
u00 = zeros(FreqDivDefectedAssembly.Mesh.nDOFs, 1);
v00 = zeros(FreqDivDefectedAssembly.Mesh.nDOFs, 1);
a00 = zeros(FreqDivDefectedAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0))

q0 = FreqDivDefectedAssembly.constrain_vector(u00);
qd0 = FreqDivDefectedAssembly.constrain_vector(v00);
qdd0 = FreqDivDefectedAssembly.constrain_vector(a00);

% Instantiate object for linear time integration
TI_lin_HFMd = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,FreqDivDefectedAssembly,F_ext);

tic
TI_lin_HFMd.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin_HFMd.Solution.ud = FreqDivDefectedAssembly.unconstrain_vector(TI_lin_HFMd.Solution.q);
FullLDurationHFMd=toc

TI_lin_HFMd_time= TI_lin_HFMd.Solution.time';
TI_lin_HFMd_dof1= TI_lin_HFMd.Solution.ud(dof1,:)';
TI_lin_HFMd_dof2= TI_lin_HFMd.Solution.ud(dof2,:)';
TI_lin_HFMd_dof3= TI_lin_HFMd.Solution.ud(dof3,:)';

save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\HFMd\TI_lin_d', ...
    'FullLDurationHFMd','TI_lin_HFMd_time','TI_lin_HFMd_dof1','TI_lin_HFMd_dof2','TI_lin_HFMd_dof3','-v7.3');
%% HFOM-d NL
tmax = 1*T;
% Instantiate object for nonlinear time integration
TI_NL_HFMd = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,FreqDivDefectedAssembly,F_ext);

% Nonlinear Time Integration

tic
TI_NL_HFMd.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL_HFMd.Solution.ud = FreqDivDefectedAssembly.unconstrain_vector(TI_NL_HFMd.Solution.q);
FullNLdurationHFMd=toc

TI_NL_HFMd_time= TI_NL_HFMd.Solution.time';
TI_NL_HFMd_dof1= TI_NL_HFMd.Solution.ud(dof1,:)';
TI_NL_HFMd_dof2= TI_NL_HFMd.Solution.ud(dof2,:)';
TI_NL_HFMd_dof3= TI_NL_HFMd.Solution.ud(dof3,:)';

save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\HFMd\TI_NL_d', ...
    'FullNLdurationHFMd','TI_NL_HFMd_time','TI_NL_HFMd_dof1','TI_NL_HFMd_dof2','TI_NL_HFMd_dof3','-v7.3');
%% Maixmum Modal interaction  MMI for defected case

%Linear Solution is Needed to Get the Best MDs
%Sorts the MDs index based on highest to lowest modal interaction (time
% depedent)
Wd=[];
Wd_names=[];
eta_d=VMd'*TI_lin_HFMd.Solution.ud;

for ii=1:size(eta_d,1)
    for jj=ii:size(eta_d,1)
        Wd_names=[ii ;jj];
        MMI= max(eta_d(ii,:).*eta_d(jj,:));
        Wd=[Wd [Wd_names ; MMI]];
    end
end

figure()

hold on

scatterbar3(Wd(1,:),Wd(2,:),Wd(3,:),1)
scatter3(Wd(1,:),Wd(2,:),Wd(3,:),'*')
grid on ;box on;
xticks(unique(Wd(1,:)))
yticks(unique(Wd(2,:)))

sortedWd=sortrows(Wd',3,'descend')';

sortedMDd=[];
for ii=1:size(sortedWd,2)
    I_index=sortedWd(1,ii);
    J_index=sortedWd(2,ii);
    MDd_ii=find(and(MDd_names(:,1)== I_index,MDd_names(:,2)== J_index));
    if isempty(MDd_ii)
        MDd_ii=find(and(MDd_names(:,1)== J_index,MDd_names(:,2)== I_index));
    end
    sortedMDd=[sortedMDd MDd_ii];
end
save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\HFMd\sortedMDd', ...
    'sortedMDd','-v7.3');
%% ROM for defected Mesh ROM-d using tensors
tmax = 700*T;
LROMdd1=[];
LROMdd2=[];
LROMdd3=[];

Ltimedd=[];
Ldurationdd=[];

NLROMdd1=[];
NLROMdd2=[];
NLROMdd3=[];
NLtimedd=[];

NLdurationdd=[];
ASEMdurationdd=[];
f0c_ROMdd={};

countdd=0;
VbModedd=[];
MDeridd=[];
%   VMss=orth(VMs);
%  MDss=orth(MDs);
for m=[15]
    for  n=[55] %m*(m+1)/2]
        t=m+n;
        if n==0
            V = [VMd(:,1:m)];
        else
            V = [VMd(:,1:m) MDd(:,[sortedMDd(1:n)])];
        end
        V=orth(V);
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(FreqDivDefectedAssembly.DATA.M)*V(:,ii));
        end
        
        
        FreqDivDefectedReducedAssembly  = ReducedAssembly(MeshDefected,V);
        
        
        %         FreqDivReducedAssembly.DATA.M = FreqDivReducedAssembly.mass_matrix();
        %         FreqDivReducedAssembly.DATA.C = FreqDivReducedAssembly.damping_matrix(0,0,u0);
        %         FreqDivReducedAssembly.DATA.K =  FreqDivReducedAssembly.stiffness_matrix(u0);
        
        FreqDivDefectedReducedAssembly.DATA.M =V'*FreqDivDefectedAssembly.DATA.M*V;
        FreqDivDefectedReducedAssembly.DATA.K =V'*FreqDivDefectedAssembly.DATA.K*V;
        
        FreqDivDefectedReducedAssembly.DATA.C = alfa*FreqDivDefectedReducedAssembly.DATA.M+beta*FreqDivDefectedReducedAssembly.DATA.K;
        
        %% using Julia
        tic
        
        tensors_ROM = reduced_tensors_ROM(FreqDivDefectedAssembly, elements, V);
        Q2=tensors_ROM.Q2;
        Q3=tensors_ROM.Q3;
        Q4=tensors_ROM.Q4;
         
        Q3t = Q3 + permute(Q3,[1 3 2]);
        Q4t = Q4 + permute(Q4,[1 3 2 4]) + permute(Q4,[1 4 2 3]);
        durationTensor=toc;
        ASEMdurationdd=[ASEMdurationdd durationTensor];
        
        %compute Natural Freq
        f0_ROM_i= sort(sqrt(eig(FreqDivDefectedReducedAssembly...
            .DATA.M\Q2))/2/pi) ;
        
        f0c_ROMdd=[f0c_ROMdd f0_ROM_i];
        
          save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROMd\V15M55ROMdTensors', ...
    'Q2','Q3','Q4','Q3t','Q4t','ASEMdurationdd','-v7.3');
        %% L sim
       
        
        % q0 = zeros(t,1);
        % qd0 = zeros(t,1);
        % qdd0 = zeros(t,1);
        q0=V'*u0;
        qd0=V'*v0;
        qdd0=V'*a0;
        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
        
        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear_tensor(q,qd,qdd,t,FreqDivDefectedReducedAssembly,F_ext,Q2);
        
        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
%         TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        Ldurationdd=[Ldurationdd durationL];
        Ltimedd=[Ltimedd TI_lin_red.Solution.time'];

        LROMdd1=[LROMdd1 (V(dof1,:)*TI_lin_red.Solution.q)'];
        LROMdd2=[LROMdd2 (V(dof2,:)*TI_lin_red.Solution.q)'];
        LROMdd3=[LROMdd3 (V(dof3,:)*TI_lin_red.Solution.q)'];
        
        
         save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROMd\LROMd', ...
    'Ldurationdd','Ltimedd','LROMdd1','LROMdd2','LROMdd3','-v7.3');
        %% NL
        
        TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,FreqDivDefectedReducedAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
        
        % time integration
        tic
        TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_newmark_red);
        %TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
        durationNL=toc
        
        NLdurationdd=[NLdurationdd durationNL];
        
        NLtimedd=[NLtimedd TI_NL_newmark_red.Solution.time'];
        %         NLROMdd1=[NLROMdd1 TI_NL_newmark_red.Solution.u(dof1,:)'];
        %         NLROMdd2=[NLROMdd2 TI_NL_newmark_red.Solution.u(dof2,:)'];
        %         NLROMdd3=[NLROMdd3 TI_NL_newmark_red.Solution.u(dof3,:)'];
        
        
        NLROMdd1=[NLROMdd1 (V(dof1,:)*TI_NL_newmark_red.Solution.q)'];
        
        NLROMdd2=[NLROMdd2 (V(dof2,:)*TI_NL_newmark_red.Solution.q)'];
        
        NLROMdd3=[NLROMdd3 (V(dof3,:)*TI_NL_newmark_red.Solution.q)'];
        
        save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROMd\NLROMd', ...
    'NLdurationdd','NLtimedd','NLROMdd1','NLROMdd2','NLROMdd3','-v7.3');
        %%
        countdd=countdd+1;
        VbModedd=[VbModedd m];
        MDeridd=[MDeridd n];
        
         save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROMd\ROMcountd', ...
    'countdd','VbModedd','MDeridd','-v7.3');
    end
end

[max_sizeT, max_indexT] = max(cellfun('size', f0c_ROMdd, 1));
f0_ROMdd=zeros(max_sizeT,size(f0c_ROMdd,2));

for ii=1:max_indexT
    f0c=f0c_ROMdd{ii};
    s=length(f0c);
    f0_ROMdd(1:s,ii)=f0c;
end
 save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\ROMd\f0_ROMdd', ...
    'f0_ROMdd','-v7.3');
%% DPROM init and Defect Sensitivities (choose the defect, 1st VM)
FORMULATION='N1' % N1/N1t/N0
VOLUME=1 ;        % defected volume =1 . nominal volume=0

% U = V(:,1:2);       % defect basis
% xi = rand(size(U,2),1)*0;
ndef=size(U,2); %number of defects
[DS, DSnames] = defect_sensitivities(FreqDivAssembly, elements, VMs, U, ...
    FORMULATION); % for DpROM
% for ii = 1:size(DS,2)
%     DS(:,ii) = DS(:,ii)/max(sqrt(sum(DS(:,ii).^2,2)));
% end
%%  DpROM VM (name of variable **DpROM) 1st VM
tmax = 700*T;
LDPROM1=[];
LDPROM2=[];
LDPROM3=[];


LtimeDP=[];
LdurationDP=[];
NLDPROM1=[];
NLDPROM2=[];
NLDPROM3=[];



NLtimeDP=[];
NLdurationDP=[];
ASEMdurationDP=[];
f0c_DPROM={};


countDP=0;
VbModeDp=[];
MDeriDp=[];
DSenti=[];
for m=[15]
    for  n=[55]
        for k=[m]%ndef*m] % always include DS or so
            
            if n==0 && k~=0
                break
            end
            
            t=m+n+k;
            
            if n==0
                V = [VMs(:,1:m)];
            else
                %                 if k<7
                V = [VMs(:,1:m) MDs(:,[sortedMDs(1:n)]) DS(:,1:k)];
                %                 elseif k>7
                %                     V = [VMs(:,1:m) MDs(:,1:n) DS(:,1:7) DS(:,9:k)];
                %                 end
            end
            V=orth(V);
            %mass normalization
            for ii = 1 : size(V, 2)
                V(:,ii) = V(:,ii) / (V(:,ii)'*(FreqDivAssembly.mass_matrix)*V(:,ii));
            end
            
            %             V=orth(V);
            
            FreqDivDefectedDPROMAssembly  = ReducedAssembly(myMesh,V);
            
            
            
            
            FreqDivDefectedDPROMAssembly.DATA.M =V'*FreqDivAssembly.DATA.M*V;
            FreqDivDefectedDPROMAssembly.DATA.K =V'*FreqDivAssembly.DATA.K*V;
            
            FreqDivDefectedDPROMAssembly.DATA.C = alfa*FreqDivDefectedDPROMAssembly.DATA.M+beta*FreqDivDefectedDPROMAssembly.DATA.K;
            
            
            
            %% using Julia
            tic
            
            tensors_DpROM = reduced_tensors_DpROM(FreqDivAssembly, elements, V,U ...
                ,FORMULATION,VOLUME);
            % evaluate the defected tensors at xi
            [Q2, Q3, Q4, Q3t, Q4t] = DefectedTensors(tensors_DpROM, xi);
            
            durationDPensor=toc;
            ASEMdurationDP=[ASEMdurationDP durationDPensor];
            %compute Natural Freq
            f0_ROM_i= sort(sqrt(eig(FreqDivDefectedDPROMAssembly...
                .DATA.M\Q2))/2/pi) ;
            
            f0c_DPROM=[f0c_DPROM f0_ROM_i];
            
             save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\DpROM\V15M55ROMDpROMTensors', ...
    'Q2','Q3','Q4','Q3t','Q4t','ASEMdurationDP','-v7.3');
            %% L
            
            % q0 = zeros(t,1);
            % qd0 = zeros(t,1);
            % qdd0 = zeros(t,1);
            q0=V'*u0;
            qd0=V'*v0;
            qdd0=V'*a0;
            TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
            
            % Modal linear Residual evaluation function handle
            Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear_tensor(q,qd,qdd,t,FreqDivDefectedDPROMAssembly,F_ext,Q2);
            
            % time integration
            tic
            TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
%             TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
            durationL=toc
            LdurationDP=[LdurationDP durationL];
            LtimeDP=[LtimeDP TI_lin_red.Solution.time'];
            
            LDPROM1=[LDPROM1  (V(dof1,:)*TI_lin_red.Solution.q)'];
            LDPROM2=[LDPROM2 (V(dof2,:)*TI_lin_red.Solution.q)'];
            LDPROM3=[LDPROM3 (V(dof3,:)*TI_lin_red.Solution.q)'];
            
            
         save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\DpROM\DpLROM', ...
    'LdurationDP','LtimeDP','LDPROM1','LDPROM2','LDPROM3','-v7.3');
            
            
            %% NL
          
            TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);
            % Modal nonlinear Residual evaluation function handle
            Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor(q,qd,qdd,t,FreqDivDefectedDPROMAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t);
            %Residual_NL_newmark_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensor_diffApp(q,qd,qdd,t,BeamAssembly,F_ext,Q2,Q3,Q4,Q3t,Q4t,Vr,Mrrr,Crrr);
            
            % time integration
            tic
            
            TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_newmark_red);
            %TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
            durationNL=toc
            
            NLdurationDP=[NLdurationDP durationNL];
            
            NLtimeDP=[NLtimeDP TI_NL_newmark_red.Solution.time'];
            %             NLDPROM1=[NLDPROM1 TI_NL_newmark_red.Solution.u(dof1,:)'];
            %             NLDPROM2=[NLDPROM2 TI_NL_newmark_red.Solution.u(dof2,:)'];
            %             NLDPROM3=[NLDPROM3 TI_NL_newmark_red.Solution.u(dof3,:)'];
            
            
            NLDPROM1=[NLDPROM1 (V(dof1,:)*TI_NL_newmark_red.Solution.q)'];
            
            NLDPROM2=[NLDPROM2 (V(dof2,:)*TI_NL_newmark_red.Solution.q)'];
            
            NLDPROM3=[NLDPROM3 (V(dof3,:)*TI_NL_newmark_red.Solution.q)'];
            
              save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\DpROM\DpNLROM', ...
    'NLdurationDP','NLtimeDP','NLDPROM1','NLDPROM2','NLDPROM3','-v7.3');
            %%
            countDP=countDP+1;
            VbModeDp=[VbModeDp m];
            MDeriDp=[MDeriDp n];
            DSenti=[DSenti k];
             save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\DpROM\DpROMcount', ...
    'countDP','VbModeDp','MDeriDp','DSenti','-v7.3');
            
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
save('C:\Users\mosta\OneDrive\Documents\GitHub\saved_var\DpROM\f0_DpROM', ...
    'f0_DpROM','-v7.3');
%% HyperReduced solution Linear/NL defected
tmax = 700*T;
LROM1Hd=[]; 
LROM2Hd=[];
LROM3Hd=[];

LtimeHd=[]; 
LdurationHd=[]; 

NLROM1Hd=[];
NLROM2Hd=[];
NLROM3Hd=[];

NLtimeHd=[];
NLdurationHd=[];
f0c_ROMHd={};

countHd=0;
VbModeHd=[];
MDeriHd=[];
DSentiHd=[];
VMssHd=orth(VMs);
MDssHd=orth(MDs);
for m=[15]
    for n=[55]% m*(m+1)/2]
        for k=m

        t=m+n+k;

        if n==0
            V = [VMs(:,1:m)];
        else
            V = [VMs(:,1:m) MDs(:,[sortedMDs(1:n)]) DS(:,1:k)];
            
        end
        %mass normalization
         V=orth(V);
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(FreqDivAssembly.DATA.M)*V(:,ii));
        end

        FreqDivReducedAssemblyDefected  = ReducedAssembly(MeshDefected,V);


%         FreqDivReducedAssemblyDefected.DATA.M = FreqDivReducedAssemblyDefected.mass_matrix();
%         FreqDivReducedAssemblyDefected.DATA.K =  FreqDivReducedAssemblyDefected.stiffness_matrix(u0);
%         FreqDivReducedAssemblyDefected.DATA.C = FreqDivReducedAssemblyDefected.damping_matrix(0,0,u0);
%         
        FreqDivReducedAssemblyDefected.DATA.M =V'*FreqDivDefectedAssembly.DATA.M*V;
        FreqDivReducedAssemblyDefected.DATA.K =V'*FreqDivDefectedAssembly.DATA.K*V;
        FreqDivReducedAssemblyDefected.DATA.C = alfa*FreqDivReducedAssemblyDefected.DATA.M+beta*FreqDivReducedAssemblyDefected.DATA.K;
        
        FreqDivReducedAssemblyDefected.DATA.Ud=U*xi;
     %% algorithm 3d
        V_H=V(:,1:m);
        Lin_sol=TI_lin.Solution.u;
        Lin_sol_snap=[];
        for ii=1:size(Lin_sol,2)
            if rem(ii,350)==0
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
        [G,b]=FreqDivReducedAssemblyDefected.constructGb(qq);
        GBconstructTime2=toc;
%         %fNNLS
%          tic
%          [x_fnnls,w_fnnls]=fnnls(G'*G,G'*b,(norm(b)*1))%10e-9));
%          fnnlsTime=toc
%          nnz(x_fnnls)

%           tic
%           [x_fnnls,w_fnnls]=fnnls(G'*G,G'*b)%0.001),0.01);
%           fnnlsTime=toc
%           nnz(x_fnnls)
% 
%          tic;[x_nnls,w_nnls,info]=nnls(full(G),full(b),struct('Accy',1)),'Tol',1/((norm(b)*10e-150))));
%          nnlTime=toc
%          nnz(x_nnls)
% %         options = optimset('TolX',(norm(b)*0.01));
% %         tic;x_lsq=lsqnonneg(G,b,options);
%         lsqTime=toc

% %         nnz(x_lsq)
% %         x_sNNLS=sNNLS(G,b,0.0035);
% %          nnz(x_sNNLS)
%         BeamReducedAssemblyDefected.DATA.elementWeights=x_nnls;
        [E_jain_d, xi_jain_d]=snnls_j(full(G),full(b),25e-9) %25e-9
        nnz(xi_jain_d)
        FreqDivReducedAssemblyDefected.DATA.elementWeights=round(xi_jain_d); %xi_jain_d; %ones(636,1)%xi_jain_d;
        %%
        %calculate NF for reduced case
        f0_ROM_i= sort(sqrt(eig(FreqDivReducedAssemblyDefected.DATA.M\FreqDivReducedAssemblyDefected...
            .DATA.K))/2/pi) ;
        f0c_ROMHd=[f0c_ROMHd f0_ROM_i];
        %
        q0 = zeros(t,1);
        qd0 = zeros(t,1);
        qdd0 = zeros(t,1);
        
        %         q0=V'*u0;
        %         qd0=V'*v0;
        %         qdd0=V'*a0;

        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,FreqDivReducedAssemblyDefected,F_ext);

        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc

        LdurationHd=[LdurationHd durationL];
        LtimeHd=[LtimeHd TI_lin_red.Solution.time'];
        LROM1Hd=[LROM1Hd TI_lin_red.Solution.u(dof1,:)'];
        LROM2Hd=[LROM2Hd TI_lin_red.Solution.u(dof2,:)'];
        LROM3Hd=[LROM3Hd TI_lin_red.Solution.u(dof3,:)'];

        %% Reduced solution Noninear

        TI_NL_newmark_red = ImplicitNewmark('timestep',h,'alpha',0.005);

        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hyper(q,qd,qdd,t,FreqDivReducedAssemblyDefected,F_ext);

        % time integration
        
        TI_NL_newmark_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_newmark_red.Solution.u = V * TI_NL_newmark_red.Solution.q;
        durationNL=toc

        NLdurationHd=[NLdurationHd durationNL];

        NLtimeHd=[NLtimeHd TI_NL_newmark_red.Solution.time'];
        NLROM1Hd=[NLROM1Hd TI_NL_newmark_red.Solution.u(dof1,:)'];
        NLROM2Hd=[NLROM2Hd TI_NL_newmark_red.Solution.u(dof2,:)'];
        NLROM3Hd=[NLROM3Hd TI_NL_newmark_red.Solution.u(dof3,:)'];
        %%
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

%% Comparison A
% Linear
%dof = 66; % random degree of freedom at which time response is compared
%
% nfsense1 = find_node(9.25e-07,91.5e-6,[],nodes); % node where to see the results
% dof1 = get_index(nfsense1, myMesh.nDOFPerNode );
% dof1=dof1(1);
% nfsense2 = find_node(148.85e-6,182.07e-6,[],nodes); % node where to see the results
% dof2 = get_index(nfsense2, myMesh.nDOFPerNode );
% dof2=dof2(1);
% nfsense3 = find_node(16.775e-6,22.85e-6,[],nodes); % node where to see the results
% dof3 = get_index(nfsense3, myMesh.nDOFPerNode );
% dof3=dof3(1);


figure;
hold on
%plot full simulation
plot(TI_lin_time, TI_lin_dof1,'DisplayName',[ 'Full linear (Newmark)' 'time=' num2str(FullLDuration) 's'])
plot(TI_lin_time, TI_NL_dof1,'DisplayName', ['Full Non linear ' ' time=' num2str(FullNLduration) 's'])
xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')

%plot full simulation for defected case
plot(TI_lin_HFMd_time, TI_lin_HFMd_dof1,'DisplayName', ['Full linear Defected ' 'time= ' num2str(FullLDurationHFMd) ])
hold on
plot(TI_NL_HFMd_time, TI_NL_HFMd_dof1,'DisplayName', ['Full nonlinear Defected ' 'time= ' num2str(FullNLdurationHFMd)])
hold on

%plot normal ROM
% for i=1:count
%     Vmode=VbMode(i);
%     Mderi=MDeri(i);
%     plot(Ltime(:,i), LROM1(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(Lduration(i)) 's'])
%     xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
%     plot(NLtime(:,i), NLROM1(:,i),'DisplayName', ['Reduced nonlinear ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(NLduration(i)) 's'])
%     xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
% end
%plot hyperreduced ROM
for i=1:countH
    VmodeH=VbModeH(i);
    MderiH=MDeriH(i);
    plot(LtimeH(:,i), LROM1H(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(VmodeH) ' MD=' num2str(MderiH) ' time=' num2str(LdurationH(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeH(:,i), NLROM1H(:,i),'DisplayName', ['Reduced nonlinear ' ' VM=' num2str(VmodeH) ' MD=' num2str(MderiH) ' time=' num2str(NLdurationH(i)) 's' ' Gb=' GBconstructTime])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot Tensor normal ROM
for i=1:countT
    VmodeT=VbModeT(i);
    MderiT=MDeriT(i);
    ASEMdurationT_i=ASEMdurationT(i);
    NLdurationT_i=NLdurationT(i);
    totalduration_i=ASEMdurationT_i+NLdurationT_i;  % time with assem time of tensors
          plot(LtimeT(:,i), LROMT1(:,i),'DisplayName', ['Reduced linear TensorApproach  ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(LdurationT(i)) 's' ])
          xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeT(:,i), NLROMT1(:,i),'DisplayName', ['Reduced nonlinear TensorApproach ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationT_i) 's' ])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end

%plot Tensor defected ROM (meshdefected)
for i=1:countdd
    Vmoded=VbModedd(i);
    Mderid=MDeridd(i);
    ASEMdurationd_i=ASEMdurationdd(i);
    NLdurationd_i=NLdurationdd(i);
    totalduration_i=ASEMdurationd_i+NLdurationd_i;  % time with assem time of tensors
    plot(Ltimedd(:,i), LROMdd1(:,i),'DisplayName', ['Reduced linear defected Mesh TensorApproach' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(Ldurationdd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimedd(:,i), NLROMdd1(:,i),'DisplayName', ['Reduced nonlinear defected Mesh TensorApproach ' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationd_i) 's' ])
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
    plot(LtimeDP(:,i), LDPROM1(:,i),'DisplayName', ['DPROM L  ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(LdurationDP(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeDP(:,i), NLDPROM1(:,i),'DisplayName', ['DPROM NL ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationT_i) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot hyperreduced defected ROM
for i=1:countHd
    VmodeHd=VbModeHd(i);
    MderiHd=MDeriHd(i);
    DSentiHd=DSentiHd(i)
    plot(LtimeHd(:,i), LROM1Hd(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' time=' num2str(LdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeHd(:,1), NLROM1Hd(:,i),'DisplayName', ['Hyper Reduced  Defected nonlinear ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' DS=' num2str(DSentiHd)  ' time=' num2str(NLdurationHd(i)) 's' ' Gb=' GBconstructTime2])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%% Comparison B



figure;
hold on
%plot full simulation
plot(TI_lin_time, TI_lin_dof2,'DisplayName',[ 'Full linear (Newmark)' 'time=' num2str(FullLDuration) 's'])
plot(TI_NL_time, TI_NL_dof2,'DisplayName', ['Full Non linear ' ' time=' num2str(FullNLduration) 's'])
xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')

%plot full simulation for defected case
plot(TI_lin_HFMd_time, TI_lin_HFMd_dof2,'DisplayName', ['Full linear Defected ' 'time= ' num2str(FullLDurationHFMd) ])
hold on
plot(TI_NL_HFMd_time, TI_NL_HFMd_dof2,'DisplayName', ['Full nonlinear Defected ' 'time= ' num2str(FullNLdurationHFMd)])
hold on

%plot normal ROM
% for i=1:count
%     Vmode=VbMode(i);
%     Mderi=MDeri(i);
%     plot(Ltime(:,i), LROM2(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(Lduration(i)) 's'])
%     xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
%     plot(NLtime(:,i), NLROM2(:,i),'DisplayName', ['Reduced nonlinear ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(NLduration(i)) 's'])
%     xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
% end
%plot hyperreduced ROM
for i=1:countH
    VmodeH=VbModeH(i);
    MderiH=MDeriH(i);
    plot(LtimeH(:,i), LROM2H(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(VmodeH) ' MD=' num2str(MderiH) ' time=' num2str(LdurationH(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeH(:,i), NLROM2H(:,i),'DisplayName', ['Reduced nonlinear ' ' VM=' num2str(VmodeH) ' MD=' num2str(MderiH) ' time=' num2str(NLdurationH(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%plot Tensor normal ROM
for i=1:countT
    VmodeT=VbModeT(i);
    MderiT=MDeriT(i);
    ASEMdurationT_i=ASEMdurationT(i);
    NLdurationT_i=NLdurationT(i);
    totalduration_i=ASEMdurationT_i+NLdurationT_i;  % time with assem time of tensors
    plot(LtimeT(:,i), LROMT2(:,i),'DisplayName', ['Reduced linear TensorApproach  ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(LdurationT(i)) 's' ])
   xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeT(:,i), NLROMT2(:,i),'DisplayName', ['Reduced nonlinear TensorApproach ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationT_i) 's' ])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end

%plot Tensor defected ROM (meshdefected)
for i=1:countdd
    Vmoded=VbModedd(i);
    Mderid=MDeridd(i);
    ASEMdurationd_i=ASEMdurationdd(i);
    NLdurationd_i=NLdurationdd(i);
    totalduration_i=ASEMdurationd_i+NLdurationd_i;  % time with assem time of tensors
    plot(Ltimedd(:,i), LROMdd2(:,i),'DisplayName', ['Reduced linear defected Mesh TensorApproach' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(Ldurationdd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimedd(:,i), NLROMdd2(:,i),'DisplayName', ['Reduced nonlinear defected Mesh TensorApproach ' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationd_i) 's' ])
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
    plot(LtimeDP(:,i), LDPROM2(:,i),'DisplayName', ['DPROM L  ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(LdurationDP(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeDP(:,i), NLDPROM2(:,i),'DisplayName', ['DPROM NL ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationT_i) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
for i=1:countHd
    VmodeHd=VbModeHd(i);
    MderiHd=MDeriHd(i);
    DSentiHd=DSentiHd(i)
    plot(LtimeHd(:,i), LROM2Hd(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' time=' num2str(LdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeHd(:,1), NLROM2Hd(:,i),'DisplayName', ['Reduced nonlinear ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' time=' num2str(NLdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%% Comparison C

figure;
hold on
%plot full simulation
plot(TI_lin_time, TI_lin_dof3,'DisplayName',[ 'Full linear (Newmark)' 'time=' num2str(FullLDuration) 's'])
plot(TI_NL_time, TI_NL_dof3,'DisplayName', ['Full Non linear ' ' time=' num2str(FullNLduration) 's'])
xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')

%plot full simulation for defected case
plot(TI_lin_HFMd_time, TI_lin_HFMd_dof3,'DisplayName', ['Full linear Defected ' 'time= ' num2str(FullLDurationHFMd) ])
hold on
plot(TI_NL_HFMd_time, TI_NL_HFMd_dof3,'DisplayName', ['Full nonlinear Defected ' 'time= ' num2str(FullNLdurationHFMd)])
hold on

%plot normal ROM
% for i=1:count
%     Vmode=VbMode(i);
%     Mderi=MDeri(i);
%     plot(Ltime(:,i), LROM3(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(Lduration(i)) 's'])
%     xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
%     plot(NLtime(:,i), NLROM3(:,i),'DisplayName', ['Reduced nonlinear ' ' VM=' num2str(Vmode) ' MD=' num2str(Mderi) ' time=' num2str(NLduration(i)) 's'])
%     xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
% end
%plot Tensor normal ROM
for i=1:countT
    VmodeT=VbModeT(i);
    MderiT=MDeriT(i);
    ASEMdurationT_i=ASEMdurationT(i);
    NLdurationT_i=NLdurationT(i);
    totalduration_i=ASEMdurationT_i+NLdurationT_i;  % time with assem time of tensors
    plot(LtimeT(:,i), LROMT3(:,i),'DisplayName', ['Reduced linear TensorApproach  ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(LdurationT(i)) 's' ])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeT(:,i), NLROMT3(:,i),'DisplayName', ['Reduced nonlinear TensorApproach ' ' VM=' num2str(VmodeT) ' MD=' num2str(MderiT) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationT_i) 's' ])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end

%plot Tensor defected ROM (meshdefected)
for i=1:countdd
    Vmoded=VbModedd(i);
    Mderid=MDeridd(i);
    ASEMdurationd_i=ASEMdurationdd(i);
    NLdurationd_i=NLdurationdd(i);
    totalduration_i=ASEMdurationd_i+NLdurationd_i;  % time with assem time of tensors
    plot(Ltimedd(:,i), LROMdd3(:,i),'DisplayName', ['Reduced linear defected Mesh TensorApproach' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(Ldurationdd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimedd(:,i), NLROMdd3(:,i),'DisplayName', ['Reduced nonlinear defected Mesh TensorApproach ' ' VM=' num2str(Vmoded) ' MD=' num2str(Mderid) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationd_i) 's' ])
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
    plot(LtimeDP(:,i), LDPROM3(:,i),'DisplayName', ['DPROM L  ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(LdurationDP(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeDP(:,i), NLDPROM3(:,i),'DisplayName', ['DPROM NL ' ' VM=' num2str(VmodeDP) ' MD=' num2str(MderiDP) ' DS=' num2str(DSsenti) ' time=' num2str(totalduration_i) 's' ' Assembly=' num2str(ASEMdurationT_i) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
for i=1:countHd
    VmodeHd=VbModeHd(i);
    MderiHd=MDeriHd(i);
    DSentiHd=DSentiHd(i)
    plot(LtimeHd(:,i), LROM3Hd(:,i),'DisplayName', ['Reduced linear ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' time=' num2str(LdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeHd(:,1), NLROM3Hd(:,i),'DisplayName', ['Reduced nonlinear ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' time=' num2str(NLdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
%% FFT
% Ts = h;
TIME_FFT=[NLtimeT(:,1)];
Ts=TIME_FFT(2,1);
X_FFT1=[NLROMT1(:,1) NLROMT2(:,1) NLROMT3(:,1)];
for i=1:3
    
    t=TIME_FFT(1:end);
    x=X_FFT1(1:end,i);
    
    figure(1000)
    hold on
    plot(t,x)
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    y = fft(x);
    fs = 1/Ts;
    f = (0:length(y)-1)*fs/length(y);
    legend('A','B','C')
    
    figure(2000)
    hold on
    plot(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
    figure(3000)
    hold on
    loglog(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    [y1,f1]=fft_n([t,x],fs);
    
    figure(4000)
    hold on
    loglog(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    figure(5000)
    hold on
    plot(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
end
%% FFT ANSYS
Ts = [2.30410000000000e-08];
TIME_FFT=[QGANSYS135868(:,1)];
X_FFT1=[QGANSYS135868(:,2) QGANSYS135868(:,3) QGANSYS135868(:,4)];
for i=1:3
    
    % t=TIME_FFT(1:56419);
    % x=X_FFT1(1:56419,i);
    t=TIME_FFT(150000:end-1);
    x=X_FFT1(150000:end-1,i);
    figure(1000)
    hold on
    plot(t,x)
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    y = fft(x);
    fs = 1/Ts;
    f = (0:length(y)-1)*fs/length(y);
    legend('A-ansys','B-ansys','C-ansys')
    
    figure(2000)
    hold on
    plot(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A-ansys','B-ansys','C-ansys')
    
    figure(3000)
    hold on
    loglog(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    legend('A-ansys','B-ansys','C-ansys')
    
    [y1,f1]=fft_n([t,x],fs);
    
    figure(4000)
    hold on
    loglog(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A-ansys','B-ansys','C-ansys')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    figure(10000)
    hold on
    plot(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A-ansys','B-ansys','C-ansys')
    
end

%% FFT Defected
% Ts = h;
TIME_FFT=[NLtimedd(:,1)];
Ts=TIME_FFT(2,1);

X_FFT1=[NLROMdd1(:,1) NLROMdd2(:,1) NLROMdd3(:,1)];
for i=1:3
    
    t=TIME_FFT(1:end);
    x=X_FFT1(1:end,i);
    
    figure(1000)
    hold on
    plot(t,x)
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    y = fft(x);
    fs = 1/Ts;
    f = (0:length(y)-1)*fs/length(y);
    legend('A','B','C')
    
    figure(2000)
    hold on
    plot(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
    figure(3000)
    hold on
    loglog(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    [y1,f1]=fft_n([t,x],fs);
    
    figure(4000)
    hold on
    loglog(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    figure(5000)
    hold on
    plot(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
end

%% FFT Defected DPROM
% Ts = h;
TIME_FFT=[NLtimeDP(:,1)];
Ts=TIME_FFT(2,1);

X_FFT1=[NLDPROM1(:,1) NLDPROM2(:,1) NLDPROM3(:,1)];
for i=1:3
    
    t=TIME_FFT(1:end);
    x=X_FFT1(1:end,i);
    
    figure(1000)
    hold on
    plot(t,x)
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    y = fft(x);
    fs = 1/Ts;
    f = (0:length(y)-1)*fs/length(y);
    legend('A','B','C')
    
    figure(2000)
    hold on
    plot(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
    figure(3000)
    hold on
    loglog(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    [y1,f1]=fft_n([t,x],fs);
    
    figure(4000)
    hold on
    loglog(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    figure(5000)
    hold on
    plot(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
end
%% FFT Defected DPROM
% Ts = h;
TIME_FFT=[NLtimeHd(:,1)];
Ts=TIME_FFT(2,1);

X_FFT1=[NLROM1Hd(:,1) NLROM2Hd(:,1) NLROM3Hd(:,1)];
for i=1:3
    
    t=TIME_FFT(1:end);
    x=X_FFT1(1:end,i);
    
    figure(1000)
    hold on
    plot(t,x)
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    y = fft(x);
    fs = 1/Ts;
    f = (0:length(y)-1)*fs/length(y);
    legend('A','B','C')
    
    figure(2000)
    hold on
    plot(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
    figure(3000)
    hold on
    loglog(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    [y1,f1]=fft_n([t,x],fs);
    
    figure(4000)
    hold on
    loglog(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    figure(5000)
    hold on
    plot(f1,abs(y1(:,2)))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
    legend('A','B','C')
    
end
figure(1000)
legend('A-ROM','B-ROM','C-ROM','A-ansys','B-ansys','C-ansys','A-HFMd','B-HFMd','C-HFMd','A-ROMd','B-ROMd','C-ROMd')
figure(2000)
legend('A-ROM','B-ROM','C-ROM','A-ansys','B-ansys','C-ansys','A-HFMd','B-HFMd','C-HFMd','A-ROMd','B-ROMd','C-ROMd')
figure(3000)
legend('A-ROM','B-ROM','C-ROM','A-ansys','B-ansys','C-ansys','A-HFMd','B-HFMd','C-HFMd','A-ROMd','B-ROMd','C-ROMd')
figure(4000)
legend('A-ROM','B-ROM','C-ROM','A-ansys','B-ansys','C-ansys','A-HFMd','B-ROMd','C-HFMd','A-ROMd','B-ROMd','C-ROMd')
figure(5000)
legend('A-ROM','B-ROM','C-ROM','A-ansys','B-ansys','C-ansys','A-HFMd','B-HFMd','C-HFMd','A-ROMd','B-ROMd','C-ROMd')
