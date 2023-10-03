% find_tail_TRI3
%
% Synthax:
% tailElement = find_tail_TRI3(elements, nodes)
%
% Description: 
% Returns the tail element, i.e., an array of 3 element numbers situated at
% the tail of the fish.
%
% INPUTS
%   - elements: a 2D array in which each row contains the nodes' indexes of
%               an element
%   - nodes:    a 2D array containing the x and y positions of each node
%
% OUTPUT:
%   - tailElement:  array of 3 node numbers situated in the element at the 
%                   tail of the fish
%   - tailIndexInElement: 
%   - tailElementBinVec
%
% Additional notes: only work if the tail ends with a single elements
%
% Last modified: 23/09/2023, Mathieu Dubied, ETH Zurich

function [tailElement, tailIndexInElement, tailElementBinVec] = find_tail_TRI3(elements, nodes)
    % find tail element
    [~,xMaxIndex] = max(nodes(:,1));
    [row,~] = find(elements==xMaxIndex);
    tailElement = elements(row,:);
    
    % find tail node index in the tail element (i.e., the position/index of
    % the node at the very end of the tail in the tail element)
    if elements(row,1)==xMaxIndex
        tailIndexInElement = 1;
    elseif elements(row,2)==xMaxIndex
        tailIndexInElement = 2;
    else
        tailIndexInElement = 3;
    end

    % Create binary vector of length nElement, with a one at the index
    % corresponding to the tail element
    tailElementBinVec = zeros(length(elements(:,1)),1);
    tailElementBinVec(row) = 1;
end
