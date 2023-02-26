function p = PlotMeshandForce(Nodes,Elements,show,F)
%--------------------------------------------------------------------------
% Purpose:
%         To plot 2D and 3D Finite Element Method Mesh with additional
%         force F. Currently this function only works for meshes with a 
%         constant number of nodes per element
% Synopsis :
%           PlotMeshandForce(Nodes,Elements,show,F)
% Variable Description:
%           Nodes - The nodal connectivity of the elements
%           Elements - Set of element with their node
%          
%           show - to dispaly nodal and element numbers
%                  0 (default) - do not display
%                  1           - display
%           F - the force(s) to display
%
% Note: a large part of this code is taken from PlotMesh
%       
% Last modified: 26/02/2023, Mathieu Dubied, ETH Zurich
%--------------------------------------------------------------------------

if nargin == 2
    show = 0 ;
end

dimension = size(Nodes,2) ;  % Dimension of the mesh
nel = size(Elements,1) ;                  % number of elements
nnode = length(Nodes) ;          % total number of nodes in system
nnel = size(Elements,2);                % number of nodes per element

if dimension == 3   % For 3D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:3)));
    
    
    if elementdim == 3 % solid in 3D when we simply plot the skin elements
        faces = getSkin3D(Elements);
        Elements = faces.';
        nel = size(Elements,1);      % total number of faces
        nnel = size(Elements,2);     % number of nodes per face
    end
    
    X = Nodes(Elements',1); X = reshape(X, nnel, nel);
    Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
    Z = Nodes(Elements',3); Z = reshape(Z, nnel, nel);
    
    if elementdim ~= 3 % surface or line in 3D        
        [X,Y,Z] = tune_coordinates(X,Y,Z);
    end
    
    p = patch(X,Y,Z,'w','FaceAlpha',1.0,'EdgeAlpha',1,...
        'EdgeColor','k','LineStyle','-','DisplayName','Mesh');
    view(3)
    hold on
    
elseif dimension == 2           % For 2D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:2))); % dimension of element
    hold on
    if elementdim == 2
        X = Nodes(Elements',1); X = reshape(X, nnel, nel);
        Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
        p = patch(X,Y,'w','DisplayName','Mesh');
    else % line
        p=plot(Nodes(:,1),Nodes(:,2),'.-k', 'Markersize',10);
    end
    
end


set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]);
rotate3d on;
axis equal;
axis off;


% plot force F on nodes
if size(Nodes,2)==2
    for i = 1:nnode
        drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'MaxHeadSize',0.3,'Color','r','LineWidth',1) ;
            x = [Nodes(i,1) Nodes(i,1)+F(i*2-1)/2];
            y = [Nodes(i,2) Nodes(i,2)+F(i*2)/2];
            drawArrow(x,y);
    end
elseif size(Nodes,2)==3
    for i = 1:nnode
        drawArrow = @(x,y,z) quiver3( x(1),y(1),z(1),x(2)-x(1),y(2)-y(1),z(2)-z(1),0,'MaxHeadSize',0.3,'Color','r','LineWidth',1) ;
            x = [Nodes(i,1) Nodes(i,1)+F(i*3-2)/4];
            y = [Nodes(i,2) Nodes(i,2)+F(i*3-1)/4];
            z = [Nodes(i,3) Nodes(i,3)+F(i*3)/4];
            drawArrow(x,y,z);
    end
end


% display Node numbers and Element numbers
if show ~= 0
    k = 1:nnode ;
    nd = k' ;
    if size(Nodes,2)==2
        for i = 1:nel
            text(X(:,i),Y(:,i),int2str(nd(Elements(i,:))),'fontsize',8,'color','k');
            text(mean(X(:,i)),mean(Y(:,i)),int2str(i),'fontsize',10,'color','r') ;
        end
    elseif size(Nodes,2)==3
        for i = 1:nel
            text(X(:,i),Y(:,i),Z(:,i),int2str(nd(Elements(i,:))),'fontsize',8,'color','k');
            text(mean(X(:,i)),mean(Y(:,i)),mean(Z(:,i)),int2str(i),'fontsize',10,'color','r') ;
        end
        p.FaceAlpha = 0.2;
    end
end
set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]);
rotate3d on;
axis equal;
axis off;
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