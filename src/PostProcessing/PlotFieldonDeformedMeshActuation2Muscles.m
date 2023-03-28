function h = PlotFieldonDeformedMeshActuation2Muscles(Nodes,Elements,ActuationElements,ActuationValue,ActuationElements2,ActuationValue2,disp,varargin)
%--------------------------------------------------------------------------
% Purpose:
%         To plot the profile of a component on deformed mesh, currently
%         this function only works for meshes with a constant number of
%         nodes per element. It highlights the 2 muscles with a color
%         scheme proportional to the actuation intensity
% Variable Description:
%           Nodes - The nodal coordinates of the mesh
%           -----> Nodes = [X Y Z]
%           Elements - The nodal connectivity of the elements
%           -----> Elements = [node1 node2......]
%           disp -  Nodal displacements [UX UY UZ]
%           Optional parameters:
%           color - the color of the mesh (black by default)
%           factor - Amplification factor (Change accordingly, trial)
%           cameraPos - camera placement using view() for 3D plot
%           upVec - set vectors that should be oriented vertically
%           component -  The components whose profile to be plotted 
%           -----> components  can be given in the following form: 
%                               1. a column vector in the order of node
%                               numbers 
%                               2. empty array [] in which case it is
%                               treated as zeros
%                               3. 'U' : taken as norm of displacements in
%                               'disp' (default)
%                               4. 'U1', 'U2' or 'U3': disp components:
%                               disp(:,1), disp(:,2) or disp (:,3)
%                               repectively.
%
% Additional notes: only tested for TRI3 elements
%
% Last modified: 22/03/2023, Mathieu Dubied, ETH Zurich
%--------------------------------------------------------------------------


%%
[meshcolor,factor,cameraPos,upVec,c] = parse_inputs(varargin{:});

nnodes = size(Nodes,1);      % number of nodes
dimension = size(Nodes,2) ;  % Dimension of the mesh
elementdim = rank(diff(Nodes(Elements(1,:),:))); % Dimension of elements

nel = size(Elements,1);      % total number of elements   
nnel = size(Elements,2);     % number of nodes per element


if isempty(c)
    c = zeros(nnodes,1);
end

hold on

if dimension == 3   % For 3D plots   
    ux = disp(:,1) ;
    uy = disp(:,2) ;
    uz = disp(:,3) ;
    d = sqrt(ux.^2 + uy.^2 + uz.^2);
    
    if elementdim == 3 % solid in 3D when we simply plot the skin elements 
        originalElements = Elements;
        [skin,~,skinElements,skinElementFaces] = getSkin3D(Elements);      
        skinFaces = skin.';
        nSkinFaces = size(skinFaces,1);      % total number of faces
        nodePerSkinFace = size(skinFaces,2);     % number of nodes per face
    end
        
        
    switch c
        case 'U'
            c = d;
        case 'U1'
            c = ux;
        case 'U2'
            c = uy;
        case 'U3'
            c = uz;
    end
    % Preparing for plot with patch (Elements contains the faces in 2D)
    % X,Y,Z are nnel-by-nel matrices containing x,y,z coordinates for
    % each node in the mesh. E.g., X(i,j) contains the x-coordinate of
    % the i-th node of the j-th element in the mesh, i.e., the node
    % with the index Elements(j,i).

    X = Nodes(skinFaces',1); X = reshape(X, nodePerSkinFace, nSkinFaces);
    Y = Nodes(skinFaces',2); Y = reshape(Y, nodePerSkinFace, nSkinFaces);
    Z = Nodes(skinFaces',3); Z = reshape(Z, nodePerSkinFace, nSkinFaces);
    
    UX = ux(skinFaces',1); UX = reshape(UX, nodePerSkinFace, nSkinFaces);
    UY = uy(skinFaces',1); UY = reshape(UY, nodePerSkinFace, nSkinFaces);
    UZ = uz(skinFaces',1); UZ = reshape(UZ, nodePerSkinFace, nSkinFaces);
    profile = c(skinFaces',1);
    profile = reshape(profile,[nodePerSkinFace length(profile)/nodePerSkinFace]);
    
    % Plotting the profile of a property on the deformed mesh
    defoX = X+factor*UX ;
    defoY = Y+factor*UY ;
    defoZ = Z+factor*UZ ;
    
    view(cameraPos); hold on;
    h{1} = patch(defoX,defoY,defoZ,profile,'EdgeColor',meshcolor,...
        'DisplayName','Deformed Mesh');
    hIdx = 2;
    for idx=1:size(ActuationElements,1)
        % check if an element is part of the skin elements and actuation
        % elements

        if ActuationElements(idx) == 1 && skinElements(idx) == 1
            % find nodes of the skin faces of this element
            faceNodes = faceFromsSkinElement(originalElements(idx,:),skinElementFaces(idx,:));
            % find idx in `skin' (i.e. skin faces) that corresponds to the
            % faceNodes we just found (order can be different)
            [i1, ~]=find(skinFaces==faceNodes(1,1));
            [i2, ~]=find(skinFaces==faceNodes(1,2));
            [i3, ~]=find(skinFaces==faceNodes(1,3));
            i12 = intersect(i1,i2);
            ColumnIdx = intersect(i12,i3);


            h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'white','EdgeColor',meshcolor);
            hIdx = hIdx + 1;

            if ActuationValue >=0
                h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'red','EdgeColor',meshcolor);
            else
                h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'blue','EdgeColor',meshcolor);
            end  
            set(h{hIdx},'FaceAlpha',abs(ActuationValue))
            hIdx = hIdx + 1;
            % possibly a second element in this face being par of the skin
            if faceNodes(2,1) ~= 0
                [i1, ~]=find(skinFaces==faceNodes(2,1));
                [i2, ~]=find(skinFaces==faceNodes(2,2));
                [i3, ~]=find(skinFaces==faceNodes(2,3));
                i12 = intersect(i1,i2);
                ColumnIdx = intersect(i12,i3);

                h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'white','EdgeColor',meshcolor);
                hIdx = hIdx + 1;
    
                if ActuationValue >=0
                    h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'red','EdgeColor',meshcolor);
                else
                    h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'blue','EdgeColor',meshcolor);
                end  
                set(h{hIdx},'FaceAlpha',abs(ActuationValue))
                hIdx = hIdx + 1;
            end

        end
        
    end
    
    camup(upVec) 
    
    


    rotate3d on;

elseif dimension == 2           % For 2D plots
    ux = disp(:,1) ;
    uy = disp(:,2) ;
    d = sqrt(ux.^2 + uy.^2);
    switch c
            case 'U'
                c = d;
            case 'U1'
                c = ux;
            case 'U2'
                c = uy;               
    end
    hold on;
    
    if elementdim == 2 % surface in 2D
        
        X = Nodes(Elements',1); X = reshape(X,[nnel nel]);
        Y = Nodes(Elements',2); Y = reshape(Y,[nnel nel]);
        UX = ux(Elements',1); UX = reshape(UX,[nnel nel]);
        UY = uy(Elements',1); UY = reshape(UY,[nnel nel]);
        profile = c(Elements',1); 
        profile = reshape(profile,[nnel length(profile)/nnel]);
        % Plotting the profile of a property on the deformed mesh
        defoX = X + factor*UX ;
        defoY = Y + factor*UY ;
        
        h{1} = patch(defoX,defoY,profile,'EdgeColor',meshcolor);
        h{2} = plot(defoX(:),defoY(:),'.','Color', meshcolor, 'Markersize',10);
        plotIndex = 2;
        for idx=1:size(ActuationElements,1)
            if ActuationElements(idx) == 1
                h{idx*2+1} = patch(defoX(idx*3-2:idx*3),defoY(idx*3-2:idx*3),'white','EdgeColor',meshcolor);
                if ActuationValue >=0
                    h{idx*2+2} = patch(defoX(idx*3-2:idx*3),defoY(idx*3-2:idx*3),'red','EdgeColor',meshcolor);
                else
                    h{idx*2+2} = patch(defoX(idx*3-2:idx*3),defoY(idx*3-2:idx*3),'blue','EdgeColor',meshcolor);
                end

                set(h{idx*2+2},'FaceAlpha',abs(ActuationValue))
            end
            plotIndex = idx*2+3;
        end
        for idx=1:size(ActuationElements2,1)
            if ActuationElements2(idx) == 1
                h{idx*2+plotIndex} = patch(defoX(idx*3-2:idx*3),defoY(idx*3-2:idx*3),'white','EdgeColor',meshcolor);
                if ActuationValue2 >=0
                    h{idx*2+plotIndex+1} = patch(defoX(idx*3-2:idx*3),defoY(idx*3-2:idx*3),'red','EdgeColor',meshcolor);
                else
                    h{idx*2+plotIndex+1} = patch(defoX(idx*3-2:idx*3),defoY(idx*3-2:idx*3),'blue','EdgeColor',meshcolor);
                end

                set(h{idx*2+plotIndex+1},'FaceAlpha',abs(ActuationValue2))
            end
            
        end
        
    else
        h = plot(Nodes(:,1)+factor*ux,Nodes(:,2)+factor*uy,'.-','Color', meshcolor, 'Markersize',10);
        c = [];
    end
    
end
        axis equal;
        axis off;
    % Colorbar Setting
    if ~isempty(c)
        SetColorbar
    end
end

function [meshcolor,factor,cameraPos,upVec,c] = parse_inputs(varargin)

defaultColor = 'k';
defaultFactor = 1;
defaultCameraPos = 3;
defaultUpVec = [0;1;0];
defaultComponent = 'U'; % plot norm of displacement
p = inputParser;

addParameter(p,'color',defaultColor)
addParameter(p,'factor',defaultFactor,@(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'cameraPos',defaultCameraPos,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'upVec',defaultUpVec,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'component',defaultComponent);
parse(p,varargin{:});

meshcolor = p.Results.color;
factor = p.Results.factor;
cameraPos = p.Results.cameraPos;
upVec = p.Results.upVec;
c = p.Results.component;

end

function nodes = faceFromsSkinElement(element,faceNumber)
    nodes = [0 0 0; 0 0 0];
    node = faceNumber(1);
    nodes(1,1) = element(node);
    for i = 1:2
        node = node + 1;
        if node == 5
            node = 1;
        end
        nodes(1,i+1) = element(node);
    end
    
    if faceNumber(2) ~= 0
        node = faceNumber(2);
    nodes(2,1) = element(node);
    for i = 1:2
        node = node + 1;
        if node == 5
            node = 1;
        end
        nodes(2,i+1) = element(node);
    end
    end
end
