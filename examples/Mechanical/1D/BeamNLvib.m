% 1D-Beam Example with NLvib
% Example based on:
% Sombroek et al. (2018). Numerical computation of nonlinear normal modes 
% in a modal derivative subspace. Computers and Structures, 195, 34–46. 
% https://doi.org/10.1016/j.compstruc.2017.08.016
%
% Created: 21 April 2021
% Jacopo Marconi, Politecnico di Milano
clear
clc
close

LOADDATA = true;

imod = 1; % mode analyze

%% parameters
% geometry
len = 1;        	% length
height = 1e-2;    	% height in the bending direction
thickness = 1e-2;	% thickness in the third dimension

% mesh
nElements = 10;
dx = len/nElements;

% material properties
E       = 210e9;    % Young's modulus
rho     = 7800;     % density
nu      = 0.3;      % nu

if LOADDATA == 1
    % load example data
    load('BeamNLvib.mat')
    fprintf('Data loaded \n')
end

%% Structural Model
% Material
myBeamMaterial  = KirchoffMaterial();
set(myBeamMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

% Element
myElementConstructor = @()BeamElement(thickness, height, myBeamMaterial); % same element all across the domain

% Mesh
x = (0:dx:len).';
nNodes = size(x,1);
Nodes = [x, zeros(nNodes,1)];

Elements = [1:nNodes-1;2:nNodes].';
BeamMesh = Mesh(Nodes);
BeamMesh.create_elements_table(Elements,myElementConstructor);

nNodes = BeamMesh.nNodes;
BeamMesh.set_essential_boundary_condition([1 nNodes],1:3,0);

ndofs = length( BeamMesh.EBC.unconstrainedDOFs );

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(BeamMesh);
M = BeamAssembly.mass_matrix();
u0 = zeros( BeamMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;

%% Eigenvalue problem
n_VMs = ndofs; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[Phi,om2] = eigs(Kc, Mc, n_VMs, 'SM');
[om, ind] = sort(sqrt(diag(om2)));
f0 = om/2/pi;
Phi = Phi(:,ind);
for ii = 1:n_VMs
    Phi(:,ii) = Phi(:,ii)/max(sqrt(Phi(1:3:end,ii).^2+Phi(2:3:end,ii).^2));
end
Phi = BeamAssembly.unconstrain_vector(Phi);

u = Phi(1:3:end, imod);
v = Phi(2:3:end, imod);
x = Nodes(:, 1);
y = Nodes(:, 2);

figure
scale = len / 3;
plot(x, y, 'k--o'); hold on
plot(x+u*scale, y+v*scale, 'b-o')
grid on; axis equal tight; 
xlabel('x [m]'); ylabel('y [m]'); title(['mode ' num2str(imod)])
drawnow

% Damping _________________________________________________________________
Qfactor = 100;
csi = 1./(2*Qfactor);   % dimensionless damping
om0 = 2*pi*f0(1);       % first eigenfrequency
alfa = 2 * om0 * csi;
D = alfa*M;             % Mass-proportinal damping: D = a*M
BeamAssembly.DATA.D = D;
Dc = BeamAssembly.constrain_matrix(D);

%% NLvib: HB+NMA
% Example adapted from "09_beam_cubicSpring_NM" from NLvib

BeamSystem = FE_system( BeamAssembly );

PHI_lin = BeamAssembly.constrain_vector(Phi);

analysis = 'NMA';

% Analysis parameters _____________________________________________________
H = 9;              % harmonic order (size of the problem is (2H+1)*ndof,
                    % 	    i.e. A0 and (Ai,Bi) for i=1...H for each dof)
N = 3*H+1;       	% number of time samples per period
log10a_s = -5;      % start vibration level (log10 of modal mass)
log10a_e = -2;      % end vibration level (log10 of modal mass)
inorm = ndofs-1;   	% coordinate for phase normalization

% INITIALIZATION __________________________________________________________
% Initial guess vector x0 = [Psi; om; del], where del is the modal
% damping ratio, estimate from underlying linear system
omi = om(imod);                 % linear eigenfrequency
phi = PHI_lin(:,imod);          % linear eigenmode
Psi = zeros((2*H+1)*ndofs, 1);  % initialize ALL the harmonics to zero
Psi(ndofs+(1:ndofs)) = phi;     % initialize the first harmonic with the 
                                % linear eigenmode. Psi(1:ndofs) are the
                                % 0th harmonic terms (static).
x0 = [Psi; omi; 0];          	% initial guess
q_scl = max(abs(Psi));        	% used for scaling (scl)

% Solve and continue w.r.t. Om ____________________________________________
ds = .02;                       % arclength parameter
Sopt=struct('dynamicDscale',1);	% set scaling and other options here
f_scl = mean(abs(Kc*phi));      % scaling of dynamic force equilibrium
fun_postprocess = {};           % (optional) feval(fun_postprocess, X),
                                %            results are stored in Sol
if LOADDATA == 0
    % solve
    [X_HB,Solinfo,Sol] = solve_and_continue(x0,...
        @(X) HB_residual(X, BeamSystem, H, N, analysis, inorm, f_scl),...
        log10a_s, log10a_e, ds, [Sopt, fun_postprocess]);
end

% Interpret solver output _________________________________________________
Psi_HB = X_HB(1:end-3,:); 	% Normalized solution. Each column contains 
                          	% 2H+1 coefficients for the ndofs, i.e.:
                          	% Psi_HB(        1 :   ndofs, j) --> A0
                          	% Psi_HB(  ndofs+1 : 2*ndofs, j) --> A1 (Re)
                           	% Psi_HB(2*ndofs+1 : 3*ndofs, j) --> B1 (Im)
                         	% ...
                            % Psi_HB((2*H-1)*ndofs+1 : 2*H*ndofs, j) --> AH
                            % Psi_HB(2*H*ndofs+1 : (2*H+1)*ndofs, j) --> BH
                            % for the j-th frequency
                            
om_HB  = X_HB(end-2,:);   	% circular frequency
del_HB = X_HB(end-1,:);  	% modal damping ratio
log10a_HB = X_HB(end,:);   	% log10 of the reference a(imod)
a_HB = 10.^log10a_HB;
Q_HB = Psi_HB .*repmat(a_HB,... % de-normalize solution
    size(Psi_HB,1), 1);

% POSTPROCESSING __________________________________________________________
% Determine total energy in the system from the displacement and velocity
% at t=0
energy = zeros(size(a_HB));

k = repmat(1:H, ndofs, 1); % harmonic numbers

for ww = 1 : length( om_HB )
    Q_HBi = Q_HB(:, ww);              	% ((2H+1)*ndofs) vector
    Qi = reshape(Q_HBi, ndofs, 2*H+1); 	% (ndofs) * (2H+1) matrix
    
    Q0 = Qi(:,1);           % constant part
    Qre = Qi(:,2:2:end);    % real part
    Qim = Qi(:,3:2:end);    % imaginary part
    Qc = Qre + 1i*Qim;
    
    % *******************************
    % adapted from EXAMPLE 9 of NLvib
    % *******************************
    % q(t) = real(sum_k(Qc*exp(i*k*W*t))), q0 = q(t=0)
    % u(t) = dq/dt = real(sum_k(Qc*(i*k*W)*exp(i*k*W*t))), u0 = u(t=0)
    q0 = Q0 + sum(Qre, 2);
    u0 = sum(Qim .*(k*om_HB(ww)),2);
    
    % the total energy is constant on the orbit (therefore we compute it
    % for t=0):
    energy(ww) = 1/2*u0'*Mc*u0 + 1/2*q0'*Kc*q0;
    % ... we should add +1/3*K3*q0^3 +1/4*K4*q0^4, but we don't have K3,K4
end
fprintf(['\nBE WARE: the energy computed here is not correct, only \n', ...
    'convenient (we don''t have third and fourth order tensors K3 \n', ...
    'and K4 of the full model to compute 1/3*K3*q^3 and 1/4*K4*q^4).\n\n'])

%% PLOT
figure('units','normalized','position',[.3 .1 .4 .8])
subplot 211
semilogy(om_HB, a_HB, 'k.-');       % Amplitude vs frequency
    xlabel('\omega [rad/s]'); ylabel('a_{HB}');
    grid on; box on; axis tight
    hold on
    plot(xlim, 10^log10a_s*[1 1],'r--')
    plot(xlim, 10^log10a_e*[1 1],'r--')
    ylim(10.^[log10a_s-.5 log10a_e+.5])
    title(['mode ' num2str(imod)])
subplot 212
semilogx(energy, om_HB, 'b.-')      % Energy vs frequency (linear)
    grid on; axis tight
    xlabel('energy [J]'); ylabel('\omega [rad/s]');
    ylim([320 500])
    xlim([1e-4 1e1])
    title('Energy-Frequency plot')

%% NLvib: HB+FR
% to do

%% NLvib: shooting+NMA
% to do

%% NLvib: shooting+FR
% to do






