function h = Plot2MusclesAndConstraints(Nodes,Elements, ...
    ActuationElements,ColorActu1,ActuationElements2,ColorActu2, ...
    ConstrainedElements,ColorConst,varargin)
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
% Additional notes: only tested for TET4 elements
%
% Last modified: 05/05/2024, Mathieu Dubied, ETH Zurich
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
axis equal


%     ux = disp(:,1) ;
%     uy = disp(:,2) ;
%     uz = disp(:,3) ;
%     d = sqrt(ux.^2 + uy.^2 + uz.^2);

if elementdim == 3 % solid in 3D when we simply plot the skin elements 
    originalElements = Elements;
    [skin,~,skinElements,skinElementFaces] = getSkin3D(Elements);      
    skinFaces = skin.';
    nSkinFaces = size(skinFaces,1);      % total number of faces
    nodePerSkinFace = size(skinFaces,2);     % number of nodes per face
end


%     switch c
%         case 'U'
%             c = d;
%         case 'U1'
%             c = ux;
%         case 'U2'
%             c = uy;
%         case 'U3'
%             c = uz;
%     end
% Preparing for plot with patch (Elements contains the faces in 2D)
% X,Y,Z are nnel-by-nel matrices containing x,y,z coordinates for
% each node in the mesh. E.g., X(i,j) contains the x-coordinate of
% the i-th node of the j-th element in the mesh, i.e., the node
% with the index Elements(j,i).

X = Nodes(skinFaces',1); X = reshape(X, nodePerSkinFace, nSkinFaces);
Y = Nodes(skinFaces',2); Y = reshape(Y, nodePerSkinFace, nSkinFaces);
Z = Nodes(skinFaces',3); Z = reshape(Z, nodePerSkinFace, nSkinFaces);

%     UX = ux(skinFaces',1); UX = reshape(UX, nodePerSkinFace, nSkinFaces);
%     UY = uy(skinFaces',1); UY = reshape(UY, nodePerSkinFace, nSkinFaces);
%     UZ = uz(skinFaces',1); UZ = reshape(UZ, nodePerSkinFace, nSkinFaces);
% profile = c(skinFaces',1);
% profile = reshape(profile,[nodePerSkinFace length(profile)/nodePerSkinFace]);

% Plotting the profile of a property on the deformed mesh
%     defoX = X+factor*UX ;
%     defoY = Y+factor*UY ;
%     defoZ = Z+factor*UZ ;
defoX = X;
defoY = Y;
defoZ = Z;

% open View->Camera Toolbar, and View->Properties Inspector (Axes,
% Viewing angles
% settings of deformed mesh w/o animation
% view(-1.5078,39.6338)   
% camproj('orthographic')
% campos([-0.949358102339665,-1.789111161270594,1.491057614878067])
% camtarget([-0.160061975241042,-0.004309698762181,0.011789814756804])
% camup([0.244316050410686,0.552461400881531,0.79692914870002])
% camva(11.0954365);
view(-0.2207,13.2656)   
camproj('orthographic')
campos([-0.927954494978962,-3.313676033697632,2.55113399138335])
camtarget([0.423664179178241,-0.25731872550064,0.017983359980954])
camup([0.2443,0.5525,0.7969])
camva(11.0954365);
rotate3d on
hold on;

% h{1} = patch(defoX,defoY,defoZ,profile,'EdgeColor',meshcolor,...
%     'DisplayName','Deformed Mesh');
h{1} = patch(defoX,defoY,defoZ,'white','EdgeColor',meshcolor,...
'DisplayName','Deformed Mesh');
% h{1} = plot(defoX(:),defoY(:),defoZ(:),'.','Color', meshcolor, 'Markersize',5);
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


        h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),ColorActu1,'EdgeColor',meshcolor);
        hIdx = hIdx + 1;


%             set(h{hIdx},'FaceAlpha',abs(ActuationValue))
%             hIdx = hIdx + 1;
        % possibly a second element in this face being par of the skin
        if faceNodes(2,1) ~= 0
            [i1, ~]=find(skinFaces==faceNodes(2,1));
            [i2, ~]=find(skinFaces==faceNodes(2,2));
            [i3, ~]=find(skinFaces==faceNodes(2,3));
            i12 = intersect(i1,i2);
            ColumnIdx = intersect(i12,i3);

            h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),ColorActu1,'EdgeColor',meshcolor);
            hIdx = hIdx + 1;

%                 h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'white','EdgeColor',meshcolor);
%                 hIdx = hIdx + 1;
%     
%                 if ActuationValue >=0
%                     h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'red','EdgeColor',meshcolor);
%                 else
%                     h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'blue','EdgeColor',meshcolor);
%                 end  
%                 set(h{hIdx},'FaceAlpha',abs(ActuationValue))
%                 hIdx = hIdx + 1;
        end

    end
end

for idx=1:size(ActuationElements2,1)
% check if an element is part of the skin elements and actuation
% elements

    if ActuationElements2(idx) == 1 && skinElements(idx) == 1
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

%             if ActuationValue2 >=0
%                 h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'red','EdgeColor',meshcolor);
%             else
%                 h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'blue','EdgeColor',meshcolor);
%             end  
%             set(h{hIdx},'FaceAlpha',abs(ActuationValue2))
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

%                 if ActuationValue2 >=0
%                     h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'red','EdgeColor',meshcolor);
%                 else
%                     h{hIdx} = patch(defoX(ColumnIdx*3-2:ColumnIdx*3),defoY(ColumnIdx*3-2:ColumnIdx*3),defoZ(ColumnIdx*3-2:ColumnIdx*3),'blue','EdgeColor',meshcolor);
%                 end  
%                 set(h{hIdx},'FaceAlpha',abs(ActuationValue2))
%                 hIdx = hIdx + 1;
        end

    end
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
