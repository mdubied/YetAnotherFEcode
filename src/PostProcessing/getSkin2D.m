% -------------------------------------------------------------------------
% getSkin2D(elements): computes the external faces of the mesh. 
% INPUTS:
%   elements: table of elements (Nelements X Ndofs)
% OUTPUTS:
%   skin: table of nodes labels of external faces. Each column contains the
%         nodes of one face.
%   allfaces: as skin, but with all the faces.
%   skinElements: logical vector of size (nElements X 1), 1 meaning the element
%                 has a face which is part of the skin.  
%
% Supported elements: QUAD4
% Note: a large portion of the code is taken from getSkin3D
% Last modified: 18/08/2022, Mathieu Dubied, ETH Zürich
% -------------------------------------------------------------------------
function [skin,allfaces,skinElements] = getSkin2D(elements)

nnel  = size(elements,2); % number of nodes per element

% select faces order
switch nnel
    case 4 % QUAD4
        faces = [1 2; 2 3; 3 4; 4 1];
end

% build face matrix
N = size(elements,1);
FACES = zeros(size(faces,2),N*size(faces,1));
count = 1;
for ii = 1:N
    vertici = elements(ii,:);
    for jj = 1:size(faces,1)
        FACES(:,count) = vertici(faces(jj,:))';
        count = count+1;
    end
end

% function [x1,icc] = OnlyNotRepeatedColumns(x)____________________________
x = FACES;
xs = sort(x,1);
[~,ia,ic] = unique(xs','rows','stable');
ia = ia';
ic = ic';
ic = ia(ic); % a unique index is given to each unique column
icc = hist(ic,1:length(ic));	% counts the occurences 
                                % (occ. after the first appear as zeros)
icc = icc == 1;                 % take only terms occurring ONCE
% x1 = xs(:,icc);__________________________________________________________

% remove all the repeated faces (internal), so to plot only external ones
skin = FACES(:,icc);
allfaces = FACES;

% find elements with faces being part of the skin
skinElements = zeros(N,1);
indexSkinElements = 1;
skinMembers = ismember(allfaces.',skin','rows');
for ii = 1:size(skinMembers,1)
    if skinMembers(ii) == 1
        skinElements(indexSkinElements) = 1;
    end 

    if mod(ii,size(faces,1)) == 0
        indexSkinElements = indexSkinElements + 1;
    end
end 


