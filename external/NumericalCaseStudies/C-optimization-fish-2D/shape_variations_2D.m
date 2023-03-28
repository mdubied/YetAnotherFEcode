% shape_variations_2D
%
% Synthax:
% [thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead,smoothTail,smoothHead] 
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
% (1) [...]:    10 shape variations
%     
%
% Last modified: 22/03/2023, Mathieu Dubied, ETH Zurich
function [thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead,smoothTail,smoothHead] = shape_variations_2D(nodes,Lx,Ly)
    % (1) thinner fish (y-axis) 
    nodes_projected = [nodes(:,1), nodes(:,2)*0];   % projection on x-axis
    yDif = (nodes_projected(:,2) - Ly*0.5.*sign(nodes(:,2))).*abs(nodes(:,2))/(Ly*0.5); % y-difference projection vs nominal
    thinFish = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    thinFish(2:2:end) = yDif;                       % fill up all y-positions
    
    % (2) shorter fish (x-axis)
    nodes_projected = [nodes(:,1)*0, nodes(:,2)];   % projection on y-axis
    xDif = nodes_projected(:,1) - nodes(:,1);       % x-difference projection vs nominal
    shortFish = zeros(numel(nodes),1);              % create a single long vectors [x1 y1 x2 y2 ...]^T
    shortFish(1:2:end-1) = xDif;                    % fill up all x-positions
    
    % (3) linear tail thinning
    xCross = Lx*0.25;
    slope = -Ly*0.5/(max(nodes(:,1))-xCross);
    b = Ly*0.5-slope*xCross;
    nodes_linear = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) >= xCross
            nodes_linear(n,:) = [nodes(n,1), slope*nodes(n,1)+b];    % projection on ellipse
            yDif(n) = abs(nodes_linear(n,2) - max(nodes(:,2))).*-sign(nodes(n,2)).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    linearTail = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    linearTail(2:2:end) = yDif;                       % fill up all y-positions
    
    % (4) long tail
    elCenter = [Lx*0.25;0];
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
    longTail = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    longTail(2:2:end) = yDif;                       % fill up all y-positions
    
    % (5) short tail
    elCenter = [Lx*0.75;0];
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
    shortTail = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    shortTail(2:2:end) = yDif;                       % fill up all y-positions
    
    % (6) linear head thinning
    xCross = Lx*0.75;
    slope = Ly*0.5/(xCross);
    b = 0;
    nodes_linear = [nodes(:,1), nodes(:,2)];
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) <= xCross
            nodes_linear(n,:) = [nodes(n,1), slope*nodes(n,1)+b];    % projection on ellipse
            yDif(n) = abs(nodes_linear(n,2) - max(nodes(:,2))).*-sign(nodes(n,2)).*abs(nodes(n,2))/(Ly*0.5);       % y-difference projection vs nominal
        end
    end
    linearHead = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    linearHead(2:2:end) = yDif;                       % fill up all y-positions
    
    % (7) long head
    elCenter = [Lx*0.75;0];
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
    longHead = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    longHead(2:2:end) = yDif;                       % fill up all y-positions
    
    % (8) short head
    elCenter = [Lx*0.25;0];
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
    shortHead = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    shortHead(2:2:end) = yDif;                       % fill up all y-positions
    
    % (9) smooth tail
    elCenter = [Lx*0.75;0];
    a = Lx*0.25;
    b = Ly*0.5;
    nodes_ellipse = [nodes(:,1), nodes(:,2)];
    xDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) >= elCenter(1)
            nodes_ellipse(n,1) = (a/b.*sqrt(b^2-(nodes(n,2)).^2)+elCenter(1));                                  % projection on ellipse
            xDif(n) = (nodes_ellipse(n,1) - max(nodes(:,1))).*abs(nodes(n,2))/a.*(abs(max(nodes(:,1)-nodes(n,1))-a)/a);   % x-difference projection vs nominal
        end
    end
    smoothTail = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    smoothTail(1:2:end) = xDif;                       % fill up all x-positions
    
    % (10) smooth head
    elCenter = [Lx*0.25;0];
    a = Lx*0.25;
    b = Ly*0.5;
    nodes_ellipse = [nodes(:,1), nodes(:,2)];
    xDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) <= elCenter(1)
            nodes_ellipse(n,1) = (-a/b.*sqrt(b^2-(nodes(n,2)).^2)+elCenter(1));                                  % projection on ellipse
            xDif(n) = nodes_ellipse(n,1).*abs(nodes(n,2))/a.*(abs(nodes(n,1)-a)/a);   % x-difference projection vs nominal
        end
    end
    smoothHead = zeros(numel(nodes),1);               % create a single long vectors [x1 y1 x2 y2 ...]^T
    smoothHead(1:2:end) = xDif;                       % fill up all x-positions

end