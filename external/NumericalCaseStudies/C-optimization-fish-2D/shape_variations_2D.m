% shape_variations_2D
%
% Synthax:
% [thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead] 
% = shape_variations_2D(nodes,Lx,Ly)
%
% Description: Defines the shape variations used for the 2D case study
%
% INPUTS:              
% (1) nodes:    nodes and their coordinates
% (2) Lx:       size of the nominal shape in the x-direction
% (3) Ly:       size of the nominal shape in the <-direction
%
% OUTPUTS:
% (1) [...]:    8 shape variations
%     
%
% Last modified: 19/10/2023, Mathieu Dubied, ETH Zurich
function [thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead] = shape_variations_2D(nodes,Lx,Ly)
    % (1) thinner fish (y-axis) 
    nodes_projected = [nodes(:,1), nodes(:,2)*0];   % projection on x-axis
    yDif = (nodes_projected(:,2) - nodes(:,2));     % y-difference projection vs nominal
    thinFish = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    thinFish(2:2:end) = yDif;                       % fill up all y-positions
    
    % (2) shorter fish (x-axis)
    nodes_projected = [nodes(:,1)*0, nodes(:,2)];   % projection on y-axis
    xDif = nodes_projected(:,1) - nodes(:,1);       % x-difference projection vs nominal
    shortFish = zeros(numel(nodes),1);              % create a single long vectors [x1 y1 x2 y2 ...]^T
    shortFish(1:2:end-1) = xDif;                    % fill up all x-positions
    
    % (3) linear tail thinning
    xCross = -Lx*0.25;
    slope = Ly*0.5/abs(min(nodes(:,1))-xCross);
    b = Ly*0.5-slope*xCross;
    nodes_linear = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) <= xCross
            nodes_linear(n,:) = [nodes(n,1), slope*nodes(n,1)+b];    % projection on ellipse
            yDif(n) = abs(nodes_linear(n,2) - max(nodes(:,2))).*-sign(nodes(n,2)).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    linearTail = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    linearTail(2:2:end) = yDif;                       % fill up all y-positions
    
    % (4) long tail
    elCenter = [-Lx*0.25;0];
    a = Lx*0.75;
    b = Ly*0.5;
    nodes_ellipse = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) <= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), sign(nodes(n,2)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
            yDif(n) = (nodes_ellipse(n,2) - max(nodes(:,2)).*sign(nodes(n,2))).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    longTail = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    longTail(2:2:end) = yDif;                       % fill up all y-positions
    
    % (5) short tail
    elCenter = [-Lx*0.75;0];
    a = Lx*0.25;
    b = Ly*0.5;
    nodes_ellipse = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) <= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), sign(nodes(n,2)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
            yDif(n) = (nodes_ellipse(n,2) - max(nodes(:,2)).*sign(nodes(n,2))).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    shortTail = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    shortTail(2:2:end) = yDif;                       % fill up all y-positions
    
    % (6) linear head thinning
    xCross = -Lx*0.75;
    slope = Ly*0.5/(xCross);
    b = 0;
    nodes_linear = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) >= xCross
            nodes_linear(n,:) = [nodes(n,1), slope*nodes(n,1)+b];    % projection on ellipse
            yDif(n) = abs(nodes_linear(n,2) - max(nodes(:,2))).*-sign(nodes(n,2)).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    linearHead = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    linearHead(2:2:end) = yDif;                       % fill up all y-positions
    
    % (7) long head
    elCenter = [-Lx*0.75;0];
    a = Lx*0.75;
    b = Ly*0.5;
    nodes_ellipse = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) >= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), sign(nodes(n,2)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
            yDif(n) = (nodes_ellipse(n,2) - max(nodes(:,2)).*sign(nodes(n,2))).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    longHead = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    longHead(2:2:end) = yDif;                       % fill up all y-positions
    
    % (8) short head
    elCenter = [-Lx*0.25;0];
    a = Lx*0.25;
    b = Ly*0.5;
    nodes_ellipse = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) >= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), sign(nodes(n,2)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
            yDif(n) = (nodes_ellipse(n,2) - max(nodes(:,2)).*sign(nodes(n,2))).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    shortHead = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    shortHead(2:2:end) = yDif;                       % fill up all y-positions
    
end