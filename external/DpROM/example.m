% EXAMPLE: beam meshed with 2D element
clear; 
close all; 
clc
format short g

FORMULATION = 'N1'; % N1/N1t/N0


%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .2;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = false;	% set "true" for plane_stress
% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 2;
Ly = .05;
nx = 40;
ny = 2;
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
n_VMs = 3; % first n_VMs modes with lowest frequency calculated

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
VMd = NominalAssembly.unconstrain_vector(VMd);

% PLOT (defected)
mod = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure
PlotMesh(nodes_defected, elementPlot, 0);
v1 = reshape(VMd(:,mod), 2, []).';
S = 2*max(nodes_defected(:,2));
PlotFieldonDeformedMesh(nodes_defected, elementPlot, v1, 'factor', S);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0d(mod),3) ' Hz'])
axis on; grid on; box on

% Damping _________________________________________________________________
alfa = 3.1;
beta = 6.3*1e-6;
D = alfa*Mn + beta*Kn;
NominalAssembly.DATA.D = D;
Dc = NominalAssembly.constrain_matrix(D);


%% Modal Derivatives & Defect Sensitivities                         

[MDn, MDnames] = modal_derivatives(NominalAssembly, elements, VMn); % nominal
MDd = modal_derivatives(DefectedAssembly, elements, VMd);           % defected

[DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    FORMULATION); % for DpROM

% PLOT a MD
nwho = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure
v1 = reshape(MDn(:, nwho), 2, []).';
S = Ly/2;
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S, 'component', 'U');
title(['\theta_{' num2str(MDnames(1,nwho)) num2str(MDnames(2,nwho)) '}'])
axis on; grid on; box on


%% (Dp)ROM                                                          

Vn = [VMn MDn DS]; 	% reduced order basis (DpROM)
Vd = [VMd MDd];   	% reduced order basis (ROM-d)

% mass-normalize
for ii = 1 : size(Vn, 2)
    Vn(:,ii) = Vn(:,ii) / (Vn(:,ii)'*Mn*Vn(:,ii));
end
for ii = 1 : size(Vd, 2)
    Vd(:,ii) = Vd(:,ii) / (Vd(:,ii)'*Md*Vd(:,ii));
end

Mnr = Vn'*Mn*Vn; 	% reduced mass matrix (DpROM)
Mdr = Vd'*Md*Vd; 	% reduced mass matrix (ROM-d)

% standard reduced order model (defects in the mesh)
tensors_ROM = reduced_tensors_ROM(DefectedAssembly, elements, Vd);

% parametric formulation for defects
VOLUME = 1;         % integration over defected (1) or nominal volume (0)
tensors_DpROM = reduced_tensors_DpROM(NominalAssembly, elements, Vn, U, ...
    FORMULATION, VOLUME);

% evaluate the defected tensors at xi
[Q2, Q3, Q4, Q3t, Q4t] = DefectedTensors(tensors_DpROM, xi);

% check eigenfrequencies
f0_ROMd = sort(sqrt(eig(Mdr\tensors_ROM.Q2))/2/pi);
f0_DpROM = sort(sqrt(eig(Mnr\Q2))/2/pi);
id = 1 : n_VMs;
f0_ROMd = f0_ROMd(id);
f0_DpROM = f0_DpROM(id);
disp(table(f0n, f0d, f0_ROMd, f0_DpROM))


%% FRF - Harmonic Balance (with NLvib)                              

% PREPARE MODEL ___________________________________________________________
% ROMd: reduced matrices
Kr = Vd' * Kd * Vd;     % reduced linear stiffness matrix
Mr = Vd' * Md * Vd;     % reduced mass matrix
Dr = Vd' * D  * Vd;     % reduced damping matrix
Fr = Vd'*Fext;          % reduced external force vector

% Let us defined the ROMd Assembly (although the ReducedAssembly class is
% not strictly necessary, we want to define a different object - remember
% that the synthax obj1=obj2 does NOT copy an object)
ROMd_Assembly = Assembly(MeshDefected, Vd); 
ROMd_Assembly.DATA.K = Kr;
ROMd_Assembly.DATA.M = Mr;
ROMd_Assembly.DATA.D = Dr;
ROMd_System = FE_system(ROMd_Assembly, Fr, 'custom');

% the function to compute the nonlinear forces and the jacobian are passed
% to NLvib through a global function handle:
global fnl_CUSTOM
fnl_CUSTOM = @(myAssembly, q) tensors_KF_NLvib(tensors_ROM.Q3,...
    tensors_ROM.Q4, tensors_ROM.Q3t, tensors_ROM.Q4t, q);

% ANALYSIS PARAMETERS _____________________________________________________
imod = 1;               % eigenfreq to study
omi = 2*pi*f0d(imod); 	% linear eigenfrequency
n = size(Vd, 2);        % number of DOFs of the reduced system
H = 5;                  % harmonic order
N = 3*H+1;              % number of time samples per period
Om_s = omi * 0.90;   	% start frequency
Om_e = omi * 1.1;    	% end frequency
ds =  1;                % Path continuation step size
exc_lev = [4000];

% COMPUTE FRs _____________________________________________________________
fprintf('\n\n FRF from %.2f to %.2f rad/s \n\n', Om_s, Om_e)
r2 = cell(length(exc_lev),1);
for iex = 1 : length(exc_lev)
    % Set excitation level
    ROMd_System.Fex1 = Fr * exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*Mr + 1i*Om_s*Dr + Kr) \ ROMd_System.Fex1;
    y0 = zeros( (2*H+1)*n , 1);
    y0( n + (1:2*n) ) = [real(Q1);-imag(Q1)];
    
    % stuff for scaling
    qscl = max(abs((-omi^2*Mr + 1i*omi*Dr + Kr) \ ROMd_System.Fex1));
    dscale = [y0*0+qscl; omi];
    Sopt = struct('Dscale', dscale, 'dynamicDscale', 1);
    
    % Solve and continue w.r.t. Om	
    [X, Solinfo, Sol] = solve_and_continue(y0, ...
        @(X) HB_residual(X, ROMd_System, H, N, 'FRF'), ...
        Om_s, Om_e, ds, Sopt);
    
    % Interpret solver output
    r2{iex} = nlvib_decode(X, Solinfo, Sol, 'FRF', 'HB', n, H);
    
    results.FRF.HB{iex} = r2{iex};
end


%% PLOT FRs                                                         

r2 = results.FRF.HB;

figure
h = 1;
for iex = 1 : length(exc_lev)
    % 1st harmonic amplitude of the forced dof (use force_dof_c!)
    A = r2{iex}.Qre(:, :, h);
    B = r2{iex}.Qim(:, :, h);
    
    c = Vd * (A + 1i*B);  % project back reduced solution to full space
    c = c(forced_dof, :);
    
    W = r2{iex}.omega;
    a1 = squeeze(sqrt( A.^2 + B.^2 ));
    plot(W, abs(c) / Ly, '.-', 'linewidth', 1); hold on
end
grid on
axis tight
xlabel('\omega [rad/s]')
ylabel('|Q_1| / L_y [-]')
title('FRF with HB (beam mid-span)')


% LINEAR RESPONSE
% compute the linear FRF for comparison
nw = 201;
w_linear = linspace(Om_s, Om_e, nw);
for iex = 1 : length(exc_lev)
    fr_linear = zeros(nNodes*2, nw);
    for ii = 1:nw
        w = w_linear(ii);
        frl = (-w^2*Mr + 1i*w*Dr + Kr) \ Fr * exc_lev(iex);
        fr_linear(:, ii) = Vd * frl;
    end
    plot(w_linear, abs(fr_linear(forced_dof, :))/Ly, 'k--')
end
drawnow


