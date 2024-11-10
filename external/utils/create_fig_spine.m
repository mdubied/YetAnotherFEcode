% ------------------------------------------------------------------------ 
% create_fig_spine.m
% 
% Description: Create a figure showing the spine and tail element of a
% given FE mesh
%
% INPUTS: 
% (1) nodes:    matrix containing the positions of the nodes
% (2) elements: matrix containing the nodes' ID of each element
%
% OUTPUTS:   
% (1) fig: a figure

% Last modified: 10/11/2024, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
function fig = create_fig_spine(elements, nodes, colorSpine, colorTail, figsize)

    fig = figure('units','centimeters','position',figsize);
    
    [~, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TET4(elements,nodes);
    [~, ~, tailElementWeights] = find_tail(elements,nodes,spineElements,nodeIdxPosInElements);
    
    PlotSpineAndTailElements(nodes,elements, ...
        spineElementWeights,colorSpine,tailElementWeights,colorTail);
    

end