% 1D-Beam Example with NLvib
%
% Created: 21 April 2021
% Jacopo Marconi, Politecnico di Milano
clear
clc
close

imod = 1; % mode analyze

%% parameters
% geometry
len = .7;        	% length
height = .014;    	% height in the bending direction
thickness = .014;  	% thickness in the third dimension

% mesh
nElements = 10;
dx = len/nElements;

% material properties
E       = 205e9;    % Young's modulus
rho     = 7800;     % density
nu      = 0.3;      % nu

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
scale = 0.3;
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

%% NLvib: HB+NMA
% Example adapted from "09_beam_cubicSpring_NM" from NLvib

BeamSystem = FE_system( BeamAssembly );

PHI_lin = BeamAssembly.constrain_vector(Phi);

analysis = 'NMA';

% Analysis parameters _____________________________________________________
H = 5;              % harmonic order (size of the problem is (2H+1)*ndof,
                    % 	    i.e. A0 and (Ai,Bi) for i=1...H for each dof)
N = 3*H+1;       	% number of time samples per period
log10a_s = -5;      % start vibration level (log10 of modal mass)
log10a_e = -3;      % end vibration level (log10 of modal mass)
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
[X_HB,Solinfo,Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X, BeamSystem, H, N, analysis, inorm, f_scl),...
    log10a_s, log10a_e, ds, [Sopt, fun_postprocess]);

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
for ww = 1 : length( om_HB )
    Q_HBi = Q_HB(:, ww);              	% ((2H+1)*ndofs) vector
    Qi = reshape(Q_HBi, ndofs, 2*H+1); 	% (ndofs) * (2H+1) matrix
    
    % displacement
    q0 = Qi(:,1)+sum(Qi(:,2:2:end),2);	% ndofs vector (sum of the real 
                                       	% coefficients A)    
	% velocity
    u0 = sum(Qi(:,3:2:end), 2)*om_HB(ww); % ndofs vector (sum of the imag.
                                        % coefficients B, times omega)
	% "linear" energy
    energy(ww) = 1/2*u0'*Mc*u0 + 1/2*q0'*Kc*q0; % +1/3*K3*q0^3 + 1/4*K4*q0^4
end

%% PLOT
figure('units','normalized','position',[.3 .1 .4 .8])
subplot 211
semilogy(om_HB/(2*pi), (a_HB), 'k.-');     % Amplitude vs frequency
xlabel('f [Hz]'); ylabel('a_{HB}');
grid on; box on; axis tight
hold on
plot(xlim, 10^log10a_s*[1 1],'r--')
plot(xlim, 10^log10a_e*[1 1],'r--')
ylim(10.^[log10a_s-.5 log10a_e+.5])
title(['mode ' num2str(imod)])

subplot 212
semilogx(energy, om_HB/(2*pi), 'k.-')       % Energy vs frequency
grid on; axis tight
xlabel('energy [J]'); ylabel('f [Hz]');
omlims = [om_HB(1) om_HB(end)];
dom = diff( omlims )/10;
omlims = [om_HB(1)-dom om_HB(end)+dom];
ylim( omlims/2/pi )

%% NLvib: HB+FR
% to do

%% NLvib: shooting+NMA
% to do

%% NLvib: shooting+FR
% to do






