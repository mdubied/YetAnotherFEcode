% build_ROM_non_optimisation
%
% Synthax:
% [V,ROM_Assembly,tensors_ROM,tensors_hydro_ROM] = build_ROM_non_optimisation(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA)
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
% (2) ROM_Assembly:         ROM assembly
% (3) tensors_ROM:          (reduced) tensors for the internal forces 
% (4) tensors_hydro_ROM:    (reduced) tensors for the hydrdynamic forces
%     
%
% Additional notes: -
%
% Last modified: 09/10/2023, Mathieu Dubied, ETH Zürich

function [V,ROM_Assembly,tensors_ROM,tailProperties,actuTop,actuBottom] = ...
    build_ROM_non_optimisation(MeshNominal,nodes,elements,USEJULIA,ACTUATION)
    
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
    n_VMs = 3;
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

    % % shape variation/defect sensitivities (PROM)
    % [DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    % FORMULATION);
    
    % ROB formulation
    % V  = [VMn MDn DS];
    mSingle = [-1 0];    % horizontal displacement, rigid body mode
    m1 = repmat(mSingle,1,nNodes)';
    mSingle = [0 1];    % vertical displacement, rigid body mode
    m2 = repmat(mSingle,1,nNodes)';
    V  = [m1 VMn MDn];
    % V = [VMn MDn]
    V  = orth(V);

    % reduced assembly
    ROM_Assembly = ReducedAssembly(MeshNominal, V);
    ROM_Assembly.DATA.M = ROM_Assembly.mass_matrix();  % reduced mass matrix 
    ROM_Assembly.DATA.C = V.'*Dn*V;    % reduced damping matrix, using C as needed by the residual function of Newmark integration
    ROM_Assembly.DATA.K = V.'*Kn*V;    % reduced stiffness matrix 

    % INTERNAL FORCES TENSORS _____________________________________________
    tensors_ROM = reduced_tensors_ROM(NominalAssembly, elements, V, USEJULIA);  
    
    % HYDRODYNAMIC FORCES TENSORS _________________________________________
    % [~,~,skinElements, skinElementFaces] = getSkin2D(elements);
    % vwater = [1;0.1];
    % rho = 1;
    % tensors_hydro_PROM = reduced_tensors_hydro_PROM(NominalAssembly, elements, V, U, FOURTHORDER, skinElements, skinElementFaces, vwater, rho);
    
    % find spine and tail elements
    [spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TRI3(elements,nodes);
    [tailNode, tailElement, tailElementWeights] = find_tail_TRI3(elements,nodes,spineElements,nodeIdxPosInElements);

    % get normalisation factors
    normalisationFactors = compute_normalisation_factors(nodes, elements, spineElements, nodeIdxPosInElements);
    wTail = normalisationFactors(tailElement);

    % get matrices for tail pressure force
    [A,B] = compute_AB_tail_pressure(nodeIdxPosInElements(tailElement,:));
    nodesTailEl = elements(tailElement,:);
    iDOFs = [nodesTailEl(1)*2-1,nodesTailEl(1)*2,...
             nodesTailEl(2)*2-1,nodesTailEl(2)*2,...
             nodesTailEl(3)*2-1,nodesTailEl(3)*2];
    VTail = V(iDOFs,:);

    % group tail quantities in a struct
    tailProperties.A = A;
    tailProperties.B = B;
    tailProperties.w = wTail;
    tailProperties.V = VTail;

    % ACTUATION FORCES MATRICES ___________________________________________
    if ACTUATION
        Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of the nominal fish
        Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of the nominal fish

        nel = size(elements,1);
        actuationDirection = [1;0;0];               %[1;0]-->[1;0;0] (Voigt notation)
        
        % top muscle
        topMuscle = zeros(nel,1);
        for el=1:nel
            elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
            elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
            if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
                topMuscle(el) = 1;
            end    
        end
        actuTop = reduced_tensors_actuation_ROM(NominalAssembly, V, topMuscle, actuationDirection);
        
        % bottom muscle
        bottomMuscle = zeros(nel,1);
        for el=1:nel
            elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
            elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
            if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
                bottomMuscle(el) = 1;
            end    
        end
        actuBottom = reduced_tensors_actuation_ROM(NominalAssembly, V, bottomMuscle, actuationDirection);
    else
        actuTop = 0;
        actuBottom = 0;
    end



end 