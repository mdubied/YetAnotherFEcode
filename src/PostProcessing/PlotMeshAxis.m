function p = PlotMeshAxis(Nodes,Elements,show)
%--------------------------------------------------------------------------
% Purpose:
%         To plot 2D and 3D Finite Element Method Mesh together with the 
%         coordinate axis, currently this function only works for meshes 
%         with a constant number of nodes per element
% Synopsis :
%           PlotMesh(coordinates,nodes)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [X Y Z]
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]
%           show - to dispaly nodal and element numbers
%                  0 (default) - do not display
%                  1           - display
%
% Additional notes: based on PlotMesh. For 3D plots, it shows element
% indexes rather than skinFace indexes as in PlotMesh. Tested for TRI3 and
% TET4 meshes. Maximum 3 faces out of 4 as skin faces in 3D
%
% Last modified: 11/03/2023, Mathieu Dubied, ETH Zurich
%--------------------------------------------------------------------------

if nargin == 2
    show = 0 ;
end

dimension = size(Nodes,2) ;  % Dimension of the mesh
nel = size(Elements,1);      % number of elements
nnode = length(Nodes) ;      % total number of nodes in system
nnel = size(Elements,2) ;    % number of nodes per element

if dimension == 3   % For 3D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:3)));
    
    
    if elementdim == 3 % solid in 3D when we simply plot the skin faces
        [skin,~,~,skinElementFaces] = getSkin3D(Elements);
        skinFaces = skin.';
        nSkinFaces = size(skinFaces,1);          % total number of skin faces 
        nodePerSkinFace = size(skinFaces,2);    % number of nodes per face
    end
    
    X = Nodes(skinFaces',1); X = reshape(X, nodePerSkinFace, nSkinFaces);
    Y = Nodes(skinFaces',2); Y = reshape(Y, nodePerSkinFace, nSkinFaces);
    Z = Nodes(skinFaces',3); Z = reshape(Z, nodePerSkinFace, nSkinFaces);
    
    if elementdim ~= 3 % surface or line in 3D        
        [X,Y,Z] = tune_coordinates(X,Y,Z);
    end
    
    p = patch(X,Y,Z,'w','FaceAlpha',1.0,'EdgeAlpha',1,...
        'EdgeColor','k','LineStyle','-','DisplayName','Mesh');
    hold on;
    % draw coordinate axis
    maxX = max(Nodes(:,1));
    maxY = max(Nodes(:,2));
    maxZ = max(Nodes(:,3));
    quiver3(0,0,0,maxZ*1.4,0,0,'r','LineWidth',1)
    quiver3(0,0,0,0,maxZ*1.4,0,'r','LineWidth',1)
    quiver3(0,0,0,0,0,maxZ*1.4,'r','LineWidth',1)
    text(maxZ*1.4,0,0, sprintf('x'),'Color','r')
    text(0,maxZ*1.4,0, sprintf('y'),'Color','r')
    text(0,0,maxZ*1.4, sprintf('z'),'Color','r')
    
    % open View->Camera Toolbar, and View->Properties Inspector (Axes,
    % Viewing angles
    view(-1.5078,39.6338)   
    camproj('orthographic')
    campos([-0.949358102339665,-1.789111161270594,1.491057614878067])
    camtarget([-0.160061975241042,-0.004309698762181,0.011789814756804])
    camup([0.244316050410686,0.552461400881531,0.79692914870002])
    camva(11.0954365)
    
    
elseif dimension == 2           % For 2D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:2))); % dimension of element
    hold on;
    if elementdim == 2
        X = Nodes(Elements',1); X = reshape(X, nnel, nel);
        Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
        p = patch(X,Y,'w','DisplayName','Mesh');
    else % line
        p=plot(Nodes(:,1),Nodes(:,2),'.-k', 'Markersize',10);
    end
    
end
% display node numbers and (2D or 3D) element numbers
if show ~= 0
    k = 1:nnode ;
    nd = k' ;
    if size(Nodes,2)==2
        for i = 1:nel
            text(X(:,i),Y(:,i),int2str(nd(Elements(i,:))),'fontsize',8,'color','k'); % node number
            text(mean(X(:,i)),mean(Y(:,i)),int2str(i),'fontsize',10,'color','r') ;  % 2D element number
        end
    elseif size(Nodes,2)==3
        % find the elements index associated to each skin face
        elementIdxSkinFace = zeros(nSkinFaces,1);
        skinElementFacesIdx = 1;
        columnCheck = 1;

        for i=1:nSkinFaces
            elementFound = 0;
            while elementFound == 0
                if skinElementFaces(skinElementFacesIdx,columnCheck) ~= 0
                    elementIdxSkinFace(i) = skinElementFacesIdx;
                    columnCheck = columnCheck + 1;
                    if columnCheck == 4
                        skinElementFacesIdx = skinElementFacesIdx + 1;
                        columnCheck = 1;
                    end
                    elementFound = 1;
                else
                    skinElementFacesIdx = skinElementFacesIdx + 1;
                    columnCheck = 1;
                end
            end   
        end
        % print node numbers and element numbers
        for i = 1:nSkinFaces
            
            text(X(:,i),Y(:,i),Z(:,i),int2str(nd(skinFaces(i,:))),'fontsize',8,'color','k'); % node number of the skin element i
            % text(mean(X(:,i)),mean(Y(:,i)),mean(Z(:,i)),int2str(elementIdxSkinFace(i)),'fontsize',10,'color','r') ; % 3D element number
        
        end
        p.FaceAlpha = 0.2;
    end
end
% set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]);
rotate3d on;
axis equal;
axis off;
hold off;
end

function [X,Y,Z] = tune_coordinates(X,Y,Z)
% ONLY for UNDEFORMED mesh:
% to avoid visualization problems when the surface is contained in
% a plane parallel to the main planes (xy,xz,yz):
if sum(X(:)==X(1))==numel(X)
    X(1)=X(1)+1e-15;
elseif sum(Y(:)==Y(1))==numel(Y)
    Y(1)=Y(1)+1e-15;
elseif sum(Z(:)==Z(1))==numel(Z)
    Z(1)=Z(1)+1e-15;
end
end