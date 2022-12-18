% build_PROM
%
% Synthax:
% [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM] = build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA)
%
% Description: Builds a PROM based on the nominal mesh and shape variation
% basis.
%
% INPUTS: 
% (1) MeshNominal:  nominal mesh converted from Abaqus              
% (2) nodes:        nodes and their coordinates
% (3) elements:     elements and corresponding nodes
% (4) U:            shape variation basis
% (5) FORMULATION:  order of the Neumann approximation (N0/N1/N1t)
% (6) VOLUME:       integration over defected (1) or nominal volume (0)
% (7) USEJULIA:     use of JULIA (1) for the computation of internal forces
%                   tensors
%
% OUTPUTS:
% (1) V:                    ROB    
% (2) PROM_Assembly:        PROM assembly
% (3) tensors_PROM:         (reduced) tensors for the internal forces 
% (4) tensors_hydro_PROM:   (reduced) tensors for the hydrdynamic forces
%     
%
% Additional notes:
%
% Last modified: 17/12/2022, Mathieu Dubied, ETH Zürich

function [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM] = build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA)
    
    % ASSEMBLY ____________________________________________________________
    NominalAssembly = Assembly(MeshNominal);
    Mn = NominalAssembly.mass_matrix();
    nNodes = size(nodes,1);
    u0 = zeros( MeshNominal.nDOFs, 1);
    [Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;

    % DAMPING _____________________________________________________________ 
    alfa = 3.1;
    beta = 6.3*1e-6;
    Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
    NominalAssembly.DATA.D = Dn;
    Dc = NominalAssembly.constrain_matrix(Dn);

    % ROB _________________________________________________________________
    
    % vibration modes
    n_VMs = 2;
    Kc = NominalAssembly.constrain_matrix(Kn);
    Mc = NominalAssembly.constrain_matrix(Mn);
    [VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
    [f0n,ind] = sort(sqrt(diag(om))/2/pi);
    VMn = VMn(:,ind);
    for ii = 1:n_VMs
        VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
    end
    VMn = NominalAssembly.unconstrain_vector(VMn);

    % modal derivatives
    [MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);

    % shape variation/defect sensitivities (PROM)
    [DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    FORMULATION);
    
    % ROB formulation
    V  = [VMn MDn DS];
    V  = orth(V);

    % reduced assembly
    PROM_Assembly = ReducedAssembly(MeshNominal, V);
    PROM_Assembly.DATA.M = PROM_Assembly.mass_matrix();  % reduced mass matrix 
    PROM_Assembly.DATA.C = V.'*Dn*V;    % reduced damping matrix, using C as needed by the residual function of Newmark integration
    PROM_Assembly.DATA.K = V.'*Kn*V;    % reduced stiffness matrix 

    % INTERNAL FORCES TENSORS _____________________________________________
    tensors_PROM = reduced_tensors_DpROM(NominalAssembly, elements, ...
    V, U, FORMULATION, VOLUME, USEJULIA);
    %[Q2, Q3, Q4, Q3t, Q4t, M] = DefectedTensors(tensors_PROM, xi_k);
    
    % HYDRODYNAMIC FORCES TENSORS _________________________________________
    [~,~,skinElements, skinElementFaces] = getSkin2D(elements);
    vwater = [1;0.1];
    rho = 1;
    fourthOrder = 0;
    tensors_hydro_PROM = reduced_tensors_hydro_PROM(NominalAssembly, elements, V, U, fourthOrder, skinElements, skinElementFaces, vwater, rho);
end