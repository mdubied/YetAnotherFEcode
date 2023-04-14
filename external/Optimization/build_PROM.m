% build_PROM
%
% Synthax:
% [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM,tensors_actu_top_PROM, tensors_actu_bottom_PROM] 
%   = build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,,FOURTHORDER,ACTUATION)
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
% (8) FOURTHORDER:  computing the 4th order tensors (1) or not (0)
% (9) ACTUATION:    model with actuation (1) or without (0)
%
% OUTPUTS:
% (1) V:                        ROB    
% (2) PROM_Assembly:            PROM assembly
% (3) tensors_PROM:             (reduced) tensors for the internal forces 
% (4) tensors_hydro_PROM:       (reduced) tensors for the hydrdynamic forces
% (5) tensors_actu_top_PROM:    (reduced) tensors for actuation forces (only
%                               computed if ACTUATION is set to 1. return 0
%                               otherwise). Top muscle
% (5) tensors_actu_bottom_PROM: (reduced) tensors for actuation forces (only
%                               computed if ACTUATION is set to 1. return 0
%                               otherwise). Bottom muscle
%   
%
% Last modified: 26/03/2023, Mathieu Dubied, ETH Zurich

function [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM,tensors_actu_top_PROM,tensors_actu_bottom_PROM] = ...
    build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION)
    
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
    %VMn = NominalAssembly.unconstrain_vector(VMn);
%     MDn = NominalAssembly.unconstrain_vector(MDn);
%     DS = NominalAssembly.unconstrain_vector(DS)
    V  = [VMn MDn DS];
    V  = orth(V);

    %with rigid body mode
    mSingle = [1 0]; % horizontal displacement
    m1 = repmat(mSingle,1,nNodes)';
    V  = [m1 VMn MDn DS];
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
    % used in B
    vwater = [1;0.1];
    rho = 997*0.01;
    % used in C
    vwater = [0.05;0.00001];%[1;0.01];
     rho = 997*0.005;
    tensors_hydro_PROM = reduced_tensors_hydro_PROM(NominalAssembly, elements, V, U, FOURTHORDER, skinElements, skinElementFaces, vwater, rho);
    
    % ACTUATION FORCES TENSORS ____________________________________________
    if ACTUATION
        Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of the nominal fish
        Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of the nominal fish

        nel = size(elements,1);
        actuationDirection = [1;0;0];%[1;0]-->[1;0;0] (Voigt notation)
        
        % top muscle
        topMuscle = zeros(nel,1);
        for el=1:nel
            elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
            elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
            if elementCenterY>0.00 &&  elementCenterX > Lx*0.25
                topMuscle(el) = 1;
            end    
        end
        tensors_actu_top_PROM = reduced_tensors_actuation_PROM(NominalAssembly, V, U, topMuscle, actuationDirection);
        
        % bottom muscle
        bottomMuscle = zeros(nel,1);
        for el=1:nel
            elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
            elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
            if elementCenterY<0.00 &&  elementCenterX > Lx*0.25
                bottomMuscle(el) = 1;
            end    
        end
        tensors_actu_bottom_PROM = reduced_tensors_actuation_PROM(NominalAssembly, V, U, bottomMuscle, actuationDirection);
    else
        tensors_actu_top_PROM = 0;
        tensors_actu_bottom_PROM = 0;
    end




end