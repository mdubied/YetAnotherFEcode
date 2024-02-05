% build_ROM_non_optimisation
%
% Synthax:
% [V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
%    build_ROM_non_optimisation(MeshNominal,nodes,elements,USEJULIA)
%
% Description: Builds a ROM based on the nominal mesh
%
% INPUTS: 
% (1) MeshNominal:  nominal mesh converted from Abaqus              
% (2) nodes:        nodes and their coordinates
% (3) elements:     elements and corresponding nodes
% (4) USEJULIA:     use of JULIA (1) for the computation of internal forces
%                   tensors 
%
% OUTPUTS:
% (1) V:                    ROB    
% (2) ROM_Assembly:         ROM assembly
% (3) tensors_ROM:          (reduced) tensors for the internal forces 
% (4) tailProperties:       properties of the tail pressure force
%                           (matrices, tail elements etc.)
% (5) spineProperties:      properties of the spine change in momentum
%                           (tensor, spine elements etc.)
% (6) draProperties         properties of the form drag forces
% (7) actuTop:              vectors and matrices related to the actuation
%                           muscle at the top
% (8) actuBottom:           vectors and matrices related to the actuation
%                           muscle at the bottom
%     
%
% Additional notes: -
%
% Last modified: 12/11/2023, Mathieu Dubied, ETH Zürich

function [V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_ROM_3D(MeshNominal,nodes,elements,USEJULIA)

    startROMBuilding = tic;
    
    % 2D or 3D? ___________________________________________________________
    fishDim = size(nodes,2);

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

    % % shape variation/defect sensitivities (PROM)
    % [DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    % FORMULATION);
    
    % ROB formulation
    % V  = [VMn MDn DS];
    mSingle = [1 0 0];    % horizontal displacement, rigid body mode
    m1 = repmat(mSingle,1,nNodes)';
    mSingle = [0 1 0];    % vertical displacement, rigid body mode
    m2 = repmat(mSingle,1,nNodes)';
%     mSingle = [0 0 1];
%     m3 = repmat(mSingle,1,nNodes)';
    V  = [m1 m2 VMn MDn];
    V  = orth(V);

    % plot
    % mod = 1;
    % if fishDim == 2
    %     elementPlot = elements(:,1:3);  % plot only corners (otherwise it's a mess)
    %     v1 = reshape(VMn(:,mod), 2, []).';
    % else
    %     elementPlot = elements(:,1:4); 
    %     v1 = reshape(VMn(:,mod), 3, []).';
    % end
    % figure('units','normalized','position',[.2 .1 .6 .8])
    % PlotMesh(nodes, elementPlot, 0);
    % PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
    % title(['\Phi_' num2str(mod)])

    % reduced assembly
    ROM_Assembly = ReducedAssembly(MeshNominal, V);
    ROM_Assembly.DATA.M = ROM_Assembly.mass_matrix();  % reduced mass matrix 
    ROM_Assembly.DATA.C = V.'*Dn*V;    % reduced damping matrix, using C as needed by the residual function of Newmark integration
    ROM_Assembly.DATA.K = V.'*Kn*V;    % reduced stiffness matrix 

    % INTERNAL FORCES TENSORS _____________________________________________
    tensors_ROM = reduced_tensors_ROM(NominalAssembly, elements, V, USEJULIA);  
    
    % HYDRODYNAMIC FORCES _________________________________________________
    
    % find spine and tail elements
    if fishDim == 2
        [spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TRI3(elements,nodes);
        [tailNode, tailElement, ~] = find_tail(elements,nodes,spineElements,nodeIdxPosInElements);
    else
        [spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TET4(elements,nodes);
        [tailNode, tailElement, ~] = find_tail(elements,nodes,spineElements,nodeIdxPosInElements);
    end
    
    % get spine normalisation factors
    normalisationFactors = compute_normalisation_factors(nodes, elements, spineElements, nodeIdxPosInElements);
    wTail = normalisationFactors(tailElement);

    % get dorsal nodes
    [~,matchedDorsalNodesIdx,~,matchedDorsalNodesZPos] = ....
        find_dorsal_nodes(elements, nodes, spineElements, nodeIdxPosInElements);

    % tail pressure force: get matrices
    [A,B] = compute_AB_tail_pressure_TET4(nodeIdxPosInElements(tailElement,:));
    nodesTailEl = elements(tailElement,:);
    iDOFs = [nodesTailEl(1)*3-2,nodesTailEl(1)*3-1,nodesTailEl(1)*3,...
             nodesTailEl(2)*3-2,nodesTailEl(2)*3-1,nodesTailEl(2)*3,...
             nodesTailEl(3)*3-2,nodesTailEl(3)*3-1,nodesTailEl(3)*3,...
             nodesTailEl(4)*3-2,nodesTailEl(4)*3-1,nodesTailEl(4)*3];
    VTail = V(iDOFs,:);
    R = [0 -1 0 0 0 0 0 0 0 0 0 0;
         1 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 -1 0 0 0 0 0 0 0;
         0 0 0 1 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 -1 0 0 0 0;
         0 0 0 0 0 0 1 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 -1 0;
         0 0 0 0 0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0];     % 90 degrees rotation counterclock-wise

    % tail pressure force: group tail quantities in a struct
    tailProperties.A = A;
    tailProperties.B = B;
    tailProperties.w = wTail;
    tailProperties.V = VTail;
    tailProperties.R = R;
    tailProperties.tailNode = tailNode;
    tailProperties.tailElement = tailElement;
    tailProperties.iDOFs = iDOFs;
    tailProperties.zDOFIdx = matchedDorsalNodesIdx(spineElements==tailElement)*3;
    tailProperties.z = matchedDorsalNodesZPos(tailElement);
    
    % spine momentum change tensor (reduced order)
    spineTensors = compute_spine_momentum_tensor_TET4(ROM_Assembly, spineElementWeights,nodeIdxPosInElements,normalisationFactors,matchedDorsalNodesZPos);
    spineProperties.tensors = spineTensors;
    spineProperties.spineNodes = spineNodes;
    spineProperties.spineElements = spineElements;
    spineProperties.nodeIdxPosInElements = nodeIdxPosInElements;
    spineProperties.dorsalNodeIdx = matchedDorsalNodesIdx;
    spineProperties.zPos = matchedDorsalNodesZPos;
     
    % drag force (reduced order)
    [~,~,skinElements, skinElementFaces] = getSkin3D(elements);
    headNode = find_node(0,0,0,nodes);
    headxDOF = 3*headNode-2;
    VHead = V(headxDOF,:);
    rho = 1000;
    kFactor = 1;
    tensors_drag = compute_drag_tensors_ROM(ROM_Assembly, skinElements, skinElementFaces, kFactor*rho,VHead) ;
    dragProperties.tensors = tensors_drag;
    dragProperties.skinElements = skinElements;
    dragProperties.skinElementFaces = skinElementFaces;

    % ACTUATION FORCES ____________________________________________________
    
    Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of the nominal fish
    
    nel = size(elements,1);
    actuationDirection = [1;0;0;0;0;0];               %[1;0;0]-->[1;0;0;0;0;0] (Voigt notation)

    % left muscle
    topMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
        if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
            topMuscle(el) = 1;
        end    
    end
    actuTop = reduced_tensors_actuation_ROM(NominalAssembly, V, topMuscle, actuationDirection);

    % right muscle
    bottomMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
        if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
            bottomMuscle(el) = 1;
        end    
    end
    actuBottom = reduced_tensors_actuation_ROM(NominalAssembly, V, bottomMuscle, actuationDirection);


    fprintf('Time to build ROM: %.2fsec\n',toc(startROMBuilding))

end 