% build_FOM_3D
%
% Synthax:
% [Assembly,tensors_FOM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
%    build_FOM_3D(MeshNominal,nodes,elements,USEJULIA)
%
% Description: Builds a FOM based on the nominal mesh
%
% INPUTS: 
% (1) MeshNominal:  nominal mesh converted from Abaqus              
% (2) nodes:        nodes and their coordinates
% (3) elements:     elements and corresponding nodes

%
% OUTPUTS:   
% (1) ROM_Assembly:         FOM assembly
% (2) tailProperties:       properties of the tail pressure force
%                           (matrices, tail elements etc.)
% (3) spineProperties:      properties of the spine change in momentum
%                           (tensor, spine elements etc.)
% (4) draProperties         properties of the form drag forces
% (5) actuTop:              vectors and matrices related to the actuation
%                           muscle at the top
% (6) actuBottom:           vectors and matrices related to the actuation
%                           muscle at the bottom
%     
%
% Additional notes: -
%
% Last modified: 06/02/2024, Mathieu Dubied, ETH Zürich

function [FOM_Assembly,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_FOM_3D(MeshNominal,nodes,elements)

    startROMBuilding = tic;
    
    % 2D or 3D? ___________________________________________________________
    fishDim = size(nodes,2);

    % ASSEMBLY ____________________________________________________________
    FOM_Assembly = Assembly(MeshNominal);
    Mn = FOM_Assembly.mass_matrix();
    nNodes = size(nodes,1);
    u0 = zeros( MeshNominal.nDOFs, 1);
    [Kn,~] = FOM_Assembly.tangent_stiffness_and_force(u0);
    % store matrices
    FOM_Assembly.DATA.K = Kn;
    FOM_Assembly.DATA.M = Mn;

    % DAMPING _____________________________________________________________ 
    alfa = 3.1;
    beta = 6.3*1e-6;
    Cn = alfa*Mn + beta*Kn; % Rayleigh damping 
    FOM_Assembly.DATA.C = Cn;
    C = FOM_Assembly.constrain_matrix(Cn);

    % INTERNAL FORCES TENSORS _____________________________________________
    % not computed for the FOM. The internal forces are kept in a 
    % non-polynomial form 
    
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
    tailProperties.R = R;
    tailProperties.tailNode = tailNode;
    tailProperties.tailElement = tailElement;
    tailProperties.iDOFs = iDOFs;
    tailProperties.zDOFIdx = matchedDorsalNodesIdx(spineElements==tailElement)*3;
    tailProperties.z = matchedDorsalNodesZPos(tailElement);
    tailProperties.mTilde = 0.25*pi*1000*(tailProperties.z*2)^2;
    
    
%     % spine momentum change tensor (reduced order)
%     spineTensors = compute_spine_momentum_tensor_TET4(Assembly, spineElementWeights,nodeIdxPosInElements,normalisationFactors,matchedDorsalNodesZPos);
%     spineProperties.tensors = spineTensors;
    spineProperties.spineNodes = spineNodes;
    spineProperties.spineElements = spineElements;
    spineProperties.nodeIdxPosInElements = nodeIdxPosInElements;
    spineProperties.dorsalNodeIdx = matchedDorsalNodesIdx;
    spineProperties.zPos = matchedDorsalNodesZPos;
    spineProperties.normalisationFactors = normalisationFactors;
     
    % drag force (reduced order)
    [~,~,skinElements, skinElementFaces] = getSkin3D(elements);
    headNode = find_node(0,0,0,nodes);
    headxDOF = 3*headNode-2;
    vecHeadX = zeros(1,MeshNominal.nDOFs);
    rho = 1000;
    kFactor = 2;
    tensors_drag = compute_drag_tensors_FOM(FOM_Assembly,skinElements,skinElementFaces,kFactor*rho) ;
    dragProperties.tensors = tensors_drag;
    dragProperties.skinElements = skinElements;
    dragProperties.skinElementFaces = skinElementFaces;
    dragProperties.headxDOF = headxDOF;

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
    
    actuTop = compute_actuation_tensors_FOM(FOM_Assembly,topMuscle,actuationDirection);

    % right muscle
    bottomMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
        if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
            bottomMuscle(el) = 1;
        end    
    end
    actuBottom = compute_actuation_tensors_FOM(FOM_Assembly, bottomMuscle, actuationDirection);

    fprintf('Time to build FOM: %.2fsec\n',toc(startROMBuilding))

end 