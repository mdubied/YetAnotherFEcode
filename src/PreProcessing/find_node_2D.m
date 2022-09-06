% findnode for the 2D case
% n = findnode(xt,yt,nodes2search)
% 
% Returns the id of a node at coordinates (xt,y).
% INPUTS:   coordinates xt,yt
%           nodes2search: matrix containing nodes' [x y]

function n = find_node_2D(xt,yt, nodes2search)

xy = repmat([xt yt],size(nodes2search,1),1);
sn = abs(nodes2search - xy);
[~,ind]=min(sum(sn,2));
node_labels = 1:size(nodes2search,1);
n = node_labels(ind);

