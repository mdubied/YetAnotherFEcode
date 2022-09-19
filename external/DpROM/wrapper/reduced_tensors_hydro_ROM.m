% reduced_tensors_hydro_ROM
%
% Synthax:
% tensors = reduced_tensors_hydro_ROM(myAssembly, elements, V, skinElements, skinElementFaces, vwater, rho)
%
% Description: This function computes the reduced order hydrodynamic
% tensors at the Assembly level, by combining (i.e., summing), the
% element-level contributions.
%
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - elements: table of the elements
%   - V: Reduced Order Basis (unconstrained)
%   - skinElements:
%   - skinElementFaces
%   - vwater
%   - rho:
% OUTPUT:
%   tensors: a struct variable with the following fields*:
%       .Tr1             
%    	.Tr2u            
%   	.Tr2udot
%       .Tr3uu
%       .Tr3uudot
%       .Tr3udotudot
%      	.time           computational time
%     
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - List of currently supported elements: 
%     TRI3
%
% Last modified: 16/09/2022, Mathieu Dubied, ETH Zürich

function tensors = reduced_tensors_hydro_ROM(myAssembly, elements, V, skinElements, skinElementFaces, vwater, rho)

t0=tic;
% data from myAssembly
nodes    = myAssembly.Mesh.nodes;                   % nodes table
nel      = myAssembly.Mesh.nElements;           	% number of elements
nnodes   = myAssembly.Mesh.nNodes;               	% number of nodes
freedofs = myAssembly.Mesh.EBC.unconstrainedDOFs;   % free DOFs

myMesh = myAssembly.Mesh;
mode = 'ELP';   % element level projection
m = size(V,2);
u0 = zeros( myMesh.nDOFs, 1);    

% create ROM object
RomAssembly = ReducedAssembly(myMesh, V);

% compute reduced tensors
disp(' REDUCED HYDRODYNAMIC TENSORS:')
fprintf(' Assembling %d elements ...', nel)

Tr1 = RomAssembly.vector_skin('T1e', 'weights', skinElements, skinElementFaces, vwater, rho);
Tr2udot = RomAssembly.matrix_skin('T2udote', 'weights', skinElements, skinElementFaces, vwater, rho);

% Q2 = RomAssembly.matrix('tangent_stiffness_and_force', u0);
% Q3 = tensor( RomAssembly.tensor('T2',[m m m],[2 3], mode));
% Q4 = tensor( RomAssembly.tensor('T3',[m m m m],[2 3 4], mode));

time = toc(t0);


fprintf(' %.2f s (%.2f s)\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)
fprintf(' SIZEs: %d \n\n', size(V,2))

tensors.Tr1 = Tr1;           
% tensors.Tr2u            
tensors.Tr2udot = Tr2udot;
% tensors.Tr3uu
% tensors.Tr3uudot
% tensors.Tr3udotudot
% tensors.time = time;

end


