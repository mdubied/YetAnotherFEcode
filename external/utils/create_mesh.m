% create_mesh
%
% Synthax:
% [Mesh, nodes, elements, nsetForBC] = create_mesh(filename, myElementConstructor, propRigid)
%
% Description: Creates a FE mesh in Matlab based on an Abaqus input input
% file. Describes a set of nodes used for the boundary conditions (BC)
%
% INPUTS: 
% (1) filename:             name of Abaqus input file              
% (2) myElementConstructor: element constructor type (only tested with Tet4)
% (3) propRigid:            proporion of the fish which is rigid (used for
%                           the boundary conditions)
%
% OUTPUTS:   
% (1) mesh:     mesh object
% (2) nodes:    matrix containing the positions of the nodes
% (3) elements: matrix containing the nodes' ID of each element
% (4) nsetBC:   ID of the nodes subject to boundary conditions (node set)
% (5) eset:     binary vector of size n_elements, 1 (0): element (not) 
%               subject to boundary conditions (element set)
%
% Additional notes: this function is not written in the most generic way.
% It purpose is to facilitate the creation of a mesh in our examples.
%
% Last modified: 02/09/2024, Mathieu Dubied, ETH Zürich
function [mesh, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid)
    
    % read Abaqus mesh
    [nodes, elements, ~, ~] = mesh_ABAQUSread(filename);

    % convert to cm to m and reduce the initial y dimension (specific to
    % our examples)
    nodes = nodes*0.01;
    nodes(:,2) = 0.8*nodes(:,2);
    nodes(:,:) = 0.5*nodes(:,:);

    % create mesh
    mesh = Mesh(nodes);
    mesh.create_elements_table(elements,myElementConstructor);
    
    % set for boundary conditions 
    nel = size(elements,1);
    [Lx, ~, ~] = mesh_dimensions(nodes);
    nsetBC = {};
    esetBC = zeros(nel,1);
    
    for el=1:nel   
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
        if elementCenterX >= -Lx*propRigid
            for n=1:size(elements,2) 
                if  ~any(cat(2, nsetBC{:}) == elements(el,n))
                    nsetBC{end+1} = elements(el,n); 
                end
            end   
            esetBC(el)=1;
        end
    end

end