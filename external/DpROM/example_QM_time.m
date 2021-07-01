% EXAMPLE: beam meshed with 2D element
clear; 
close all; 
clc
format short g

FORMULATION = 'N1'; % N1/N1t/N0


%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2800;     % density [kg/m^3]
nu      = 0.30;     % Poisson's ratio 
thickness = .2;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = false;	% set "true" for plane_stress
% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 2;
Ly = .050;
h = 1;
nx = 20*h;
ny = 1*h;
[nodes, elements, nset] = mesh_2Drectangle(Lx, Ly, nx, ny);

% nominal mesh
MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);
MeshNominal.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)

% defected mesh
    % arch defect
    xi = 1;                                     % defect amplitude
    yd = Ly * sin(pi/Lx * nodes(:,1));          % y-displacement 
    nodes_defected = nodes + [yd*0 yd]*xi;   	% new nodes
    arc_defect = zeros(numel(nodes),1);         
    arc_defect(2:2:end) = yd;                   % vectorized defect-field
    U = arc_defect;                             % defect basis
MeshDefected = Mesh(nodes_defected);
MeshDefected.create_elements_table(elements,myElementConstructor);
MeshDefected.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
fprintf(' Arc defect: %.1f * beam thickness \n\n', xi)

% ASSEMBLY ________________________________________________________________
% nominal
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;

% defected
DefectedAssembly = Assembly(MeshDefected);
Md = DefectedAssembly.mass_matrix();
[Kd,~] = DefectedAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    DefectedAssembly.DATA.K = Kd;
    DefectedAssembly.DATA.M = Md;

% External force __________________________________________________________
node_dofs = MeshNominal.get_DOF_from_location([Lx/2, Ly/2]);
forced_dof = node_dofs(2);
Fext = zeros(nNodes*2, 1);
Fext( forced_dof ) = 1;
Fextc = NominalAssembly.constrain_vector( Fext );

% Let us also define the index of the forced dof in the constrained vector:
forced_dof_c = NominalAssembly.free2constrained_index( forced_dof );    


%% Eigenmodes                                                       

% Eigenvalue problem_______________________________________________________
n_VMs = 4; % first n_VMs modes with lowest frequency calculated

% Vibration Modes (VM): nominal
Knc = DefectedAssembly.constrain_matrix(Kn);
Mnc = DefectedAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(Knc, Mnc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssembly.unconstrain_vector(VMn);

% Vibration Modes (VM): defected
Kdc = DefectedAssembly.constrain_matrix(Kd);
Mdc = DefectedAssembly.constrain_matrix(Md);
[VMd,om] = eigs(Kdc, Mdc, n_VMs, 'SM');
[f0d,ind] = sort(sqrt(diag(om))/2/pi);
VMd = VMd(:,ind);
for ii = 1:n_VMs
    VMd(:,ii) = VMd(:,ii)/max(sqrt(sum(VMd(:,ii).^2,2)));
end
VMd = DefectedAssembly.unconstrain_vector(VMd);

% % PLOT (defected)
% mod = 1;
% elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
% figure
% PlotMesh(nodes_defected, elementPlot, 0);
% v1 = reshape(VMd(:,mod), 2, []).';
% S = 2*max(nodes_defected(:,2));
% PlotFieldonDeformedMesh(nodes_defected, elementPlot, v1, 'factor', S);
% title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0d(mod),3) ' Hz'])
% axis on; grid on; box on

% Damping _________________________________________________________________
alfa = 3.1;
beta = 0;%6.3*1e-6;
D = alfa*Mn + beta*Kn;
NominalAssembly.DATA.D = D;
DefectedAssembly.DATA.C = D;
Dc = NominalAssembly.constrain_matrix(D);


%% Modal Derivatives & Defect Sensitivities                         

[MDn, MDnames] = modal_derivatives(NominalAssembly, elements, VMn); % nominal
MDd = modal_derivatives(DefectedAssembly, elements, VMd);           % defected

% [DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
%     FORMULATION); % for DpROM

% % PLOT a MD
% nwho = 1;
% elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
% figure
% v1 = reshape(MDn(:, nwho), 2, []).' / max(abs(v1(:)));
% S = Ly/ny;
% PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S, 'component', 'U');
% title(['\theta_{' num2str(MDnames(nwho,1)) num2str(MDnames(nwho,2)) '}'])
% axis on; grid on; box on
% drawnow


%% Quadratic Manifold                                               

Phi = VMd;
for ii = 1 : size(VMd, 2)
    Phi(:,ii) = VMd(:,ii) / sqrt(VMd(:,ii)'*Md*VMd(:,ii));
end
Th = modal_derivatives(DefectedAssembly, elements, Phi);

Theta = QM_Theta_from_SMDs(Th, MDnames);

stiff_tensors_QM = reduced_tensors_QMROM(DefectedAssembly, elements, Phi, Theta);
mass_tensors_QM  = QM_mass_tensors(DefectedAssembly, Phi, Theta);
damp_tensors_QM  = QM_damp_tensors(alfa, beta, DefectedAssembly, Phi, Theta);

QMtensors.K = stiff_tensors_QM;
QMtensors.C = damp_tensors_QM;
QMtensors.M = mass_tensors_QM;


%% (Dp)ROM                                                          

% V  = [VMn MDn DS]; 	% reduced order basis (DpROM)
Vn = [VMn MDn];   	% reduced order basis (ROM-n)
Vd = [VMd MDd];   	% reduced order basis (ROM-d)

% mass-normalize
for ii = 1 : size(Vn, 2)
    Vn(:,ii) = Vn(:,ii) / sqrt(Vn(:,ii)'*Mn*Vn(:,ii));
end
% for ii = 1 : size(V, 2)
%     V(:,ii) = V(:,ii) / sqrt(V(:,ii)'*Mn*V(:,ii));
% end
for ii = 1 : size(Vd, 2)
    Vd(:,ii) = Vd(:,ii) / sqrt(Vd(:,ii)'*Md*Vd(:,ii));
end

% Mr  = V'*Mn*V;      % reduced mass matrix (DpROM)
Mnr = Vn'*Md*Vn; 	% reduced mass matrix (ROM-n)
Mdr = Vd'*Md*Vd; 	% reduced mass matrix (ROM-d)

% standard reduced order model (defects in the mesh)
tensors_ROM = reduced_tensors_ROM(DefectedAssembly, elements, Vd);

% % parametric formulation for defects
% VOLUME = 1;         % integration over defected (1) or nominal volume (0)
% tensors_DpROM = reduced_tensors_DpROM(NominalAssembly, elements, Vn, U, ...
%     FORMULATION, VOLUME);
% 
% % evaluate the defected tensors at xi
% [Q2, Q3, Q4, Q3t, Q4t] = DefectedTensors(tensors_DpROM, xi);
% 
% % check eigenfrequencies
% f0_ROMd = sort(sqrt(eig(Mdr\tensors_ROM.Q2))/2/pi);
% f0_DpROM = sort(sqrt(eig(Mnr\Q2))/2/pi);
% id = 1 : n_VMs;
% f0_ROMd = f0_ROMd(id);
% f0_DpROM = f0_DpROM(id);
% disp(table(f0n, f0d, f0_ROMd, f0_DpROM))


%% transient analysis                                               

% SETTINGS
% excitation type
LOADING = 'impulse';

% load amplification factor
amplification_factor = 15000;

switch LOADING
    case 'harmonic'
        % forcing frequency of the average of first two natural frequencies
        omega_ext = mean(2*pi*f0d(1)); 
        T =  2*pi/omega_ext; % time period of forcing
        
        % forcing function
        F_ext = @(t) amplification_factor * Fext * sin(omega_ext * t);
        
        % time step for integration
        h = T/100;        
    case 'impulse'
        % shortest period
        T = 1/max(f0d);
        
        % time step for integration
        h = T/100;
        
        % pulse duration
        figure
        for ii = 1%[1 5 10]
            DeltaT = ii*h;
            F_ext = @(t) Fext * (t<=DeltaT);
            
            t = 0 : h : 0.05;
            y = F_ext(t);
            y = y(forced_dof,:);
            [freq, spectrummag, spectrumangle] = fft_yFs(y, 1/h);
            
            subplot 211
            plot(t, y,'linewidth',1,'DisplayName',[num2str(ii) ' cazzesimi']); hold on
            subplot 212
            plot(freq, spectrummag,'linewidth',1); hold on;
        end
        subplot 211; xlabel('time'); title(['impulse (T=' num2str(t(end)) ')']); ylabel('A')
        xlim([0 DeltaT*2]); ylim([-.1 1.1]); grid on; legend
        subplot 212; xlabel('freq'); ylabel('A'); grid on
        plot(max(f0d)*[1 1], ylim, 'k--')
        
        F_ext = @(t) amplification_factor / max(spectrummag) * Fext * (t<=DeltaT);
end


%% transient analysis (Linear & NL-LM)                              
% integtation time limit
tmax = 5 * 1/f0d(1);


% % % FULL ORDER MODEL INTEGRATION ________________________________________
% % Initial condition: equilibrium
% u0 = zeros(DefectedAssembly.Mesh.nDOFs, 1);
% v0 = zeros(DefectedAssembly.Mesh.nDOFs, 1);
% a0 = zeros(DefectedAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0)) 
% 
% q0 = DefectedAssembly.constrain_vector(u0);
% qd0 = DefectedAssembly.constrain_vector(v0);
% qdd0 = DefectedAssembly.constrain_vector(a0);
% 
% % linear
% residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,DefectedAssembly,F_ext);
% TI_linF = GeneralizedAlpha('timestep',h,'rho_inf',0.7, 'linear',true);
% TI_linF.Integrate(q0,qd0,qdd0,tmax,residual_lin);
% TI_linF.Solution.u = DefectedAssembly.unconstrain_vector(TI_linF.Solution.q);
% 
% % nonlinear
% residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,DefectedAssembly,F_ext);
% TI_NLF = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
% TI_NLF.Integrate(q0,qd0,qdd0,tmax,residual);
% TI_NLF.Solution.u = DefectedAssembly.unconstrain_vector(TI_NLF.Solution.q);


% REDUCED ORDER MODEL INTEGRATION _________________________________________
m = size(Vd,2);
myReducedAssembly = ReducedAssembly(MeshDefected, Vd);

myReducedAssembly.DATA.M = Vd'*Md*Vd;
myReducedAssembly.DATA.C = Vd'*D*Vd;
myReducedAssembly.DATA.K = Vd'*Kd*Vd;
myReducedAssembly.DATA.tensors = tensors_ROM;

q0 = zeros(m,1);
qd0 = zeros(m,1);
qdd0 = zeros(m,1);


% LINEAR ******************************************************************
% Instantiate object for time integration
TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
% Modal linear Residual evaluation function handle
Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,...
    myReducedAssembly,F_ext);
% time integration
TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
TI_lin_red.Solution.u = Vd * TI_lin_red.Solution.q;


% NONLINEAR ***************************************************************
% Instantiate object for time integration
TI_NL_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
% TI_NL_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',false);
% Modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensors(q,...
    qd,qdd,t,myReducedAssembly,F_ext);
% time integration
TI_NL_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_red.Solution.u = Vd * TI_NL_red.Solution.q;


%% transient analysis (NL-QM)                                       

myReducedAssembly.DATA.tensors = QMtensors;
myReducedAssembly.DATA.Phi = Phi;
myReducedAssembly.DATA.Theta = Theta;

m = size(Phi,2);
q0 = zeros(m,1);
qd0 = zeros(m,1);
qdd0 = zeros(m,1);

% Instantiate object for time integration
TI_QM_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
% TI_NL_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',false);
% Modal nonlinear Residual evaluation function handle
Residual_QM_red = @(q,qd,qdd,t)residual_reduced_nonlinear_tensors_QM(q,...
    qd,qdd,t,myReducedAssembly,F_ext);
% time integration
TI_QM_red.Integrate(q0,qd0,qdd0,tmax,Residual_QM_red);

q = TI_QM_red.Solution.q;
uu = zeros(size(Phi,1), size(q, 2));
for tt = 1 : size(q, 2)
    uu(:,tt) = einsum('Iij,iJ,jK->IJK', Theta, q(:,tt), q(:,tt));
end
TI_QM_red.Solution.u = Phi * q + 1/2*uu;


%% Plot all responses                                               

dof = forced_dof;

D = Ly; % characteristic dimension (thickness)

figure;
% Linear
plot(TI_lin_red.Solution.time, TI_lin_red.Solution.u(dof,:)/D, ...
    'linewidth',2,'DisplayName', ['Linear (' num2str(size(TI_lin_red.Solution.q,1)) ' dofs, ' num2str(TI_lin_red.Solution.soltime) ' s)'])
hold on
% Nonlinear LM
plot(TI_NL_red.Solution.time, TI_NL_red.Solution.u(dof,:)/D, ...
    'linewidth',2,'DisplayName', ['LM (' num2str(size(TI_NL_red.Solution.q,1)) ' dofs, ' num2str(TI_NL_red.Solution.soltime) ' s)'])
% Nonlinear QM
plot(TI_QM_red.Solution.time, TI_QM_red.Solution.u(dof,:)/D, '--', ...
    'linewidth',2,'DisplayName', ['QM (' num2str(size(TI_QM_red.Solution.q,1)) ' dofs, ' num2str(TI_QM_red.Solution.soltime) ' s)'])
xlabel('time'); ylabel('u/D'); 
grid on; axis tight; legend('show')
% Nonlinear FULL
plot(TI_NLF.Solution.time, TI_NLF.Solution.u(dof,:)/D, 'k-', ...
    'linewidth',1,'DisplayName', ['FOM (' num2str(size(TI_NLF.Solution.q,1)) ' dofs, ' num2str(TI_NLF.Solution.soltime) ' s)'])















%% auxiliary functions                                              

function Theta = QM_Theta_from_SMDs(SMDs, names)

    n = size(SMDs,1);
    m = max(names(:));
    I = names(:,1);
    J = names(:,2);

    Theta = zeros(n,m,m);
    for k = 1:size(SMDs,2)
        Theta(:,I(k),J(k)) = SMDs(:,k);
        Theta(:,J(k),I(k)) = SMDs(:,k);
    end

end

function [r, drdqdd,drdqd,drdq, c0] = residual_reduced_nonlinear_tensors( q, qd, qdd, t, Assembly, Fext)
    ten = Assembly.DATA.tensors;
    M_V = Assembly.DATA.M;
    C_V = Assembly.DATA.C;
    V = Assembly.V;
    [K_V,F_V] = tensors_KF(ten.Q2, ten.Q3, ten.Q4, ten.Q3t, ten.Q4t, q);
    % Residual is computed according to the formula above:
    F_inertial = M_V * qdd;
    F_damping = C_V * qd;
    F_ext_V =  V.'*Fext(t);
    r = F_inertial + F_damping + F_V - F_ext_V ;
    drdqdd = M_V;
    drdqd = C_V;
    drdq = K_V;
    % We use the following measure to comapre the norm of the residual $\mathbf{r}$
    % 
    % $$\texttt{c0} = \|\mathbf{M_V}\ddot{\mathbf{q}}\| + \|\mathbf{C_V}\dot{\mathbf{q}}\| 
    % + \|\mathbf{F_V}(\mathbf{q})\| + \|\mathbf{V}^{\top}\mathbf{F}_{ext}(t)\|$$
    c0 = norm(F_inertial) + norm(F_damping) + norm(F_V) + norm(F_ext_V);
end

function [r, drdqdd,drdqd,drdq, c0] = residual_reduced_nonlinear_tensors_QM( q, qd, qdd, t, Assembly, Fext)
    
    Phi = Assembly.DATA.Phi;
    Theta = Assembly.DATA.Theta;
    qmt = Assembly.DATA.tensors;
    [K,D,M,F_V,F_damping,F_inertial,F_ext_V] = tensors_KF_QM(qmt,Phi,Theta,Fext(t),q,qd,qdd);
    
    r = F_V + F_damping + F_inertial - F_ext_V;
    drdqdd = M;
    drdqd  = D;
    drdq   = K;
    % We use the following measure to comapre the norm of the residual $\mathbf{r}$
    % 
    % $$\texttt{c0} = \|\mathbf{M_V}\ddot{\mathbf{q}}\| + \|\mathbf{C_V}\dot{\mathbf{q}}\| 
    % + \|\mathbf{F_V}(\mathbf{q})\| + \|\mathbf{V}^{\top}\mathbf{F}_{ext}(t)\|$$
    c0 = norm(F_inertial) + norm(F_damping) + norm(F_V) + norm(F_ext_V);
end
