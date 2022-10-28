% tensors_hydro_FOM
%
% Synthax:
% tensors = tensors_hydro_FOM(myAssembly, elements, skinElements, skinElementFaces, vwater, rho)
%
% Description: 
% This function computes the hydrodynamic tensors at the
% Assembly level, by combining the element-level tensors. It uses classical
% element-to-assembly combining method, and results in FOM (large) tensors.
%
%
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - elements: table of the elements
%   - skinElements:
%   - skinElementFaces
%   - vwater
%   - rho:
%
% OUTPUT:
%   tensors: a struct variable with the following fields*:
%       .T1
%       .T2
%       .Tu2
%       .Tu3
%       .Tudot2
%       .Tudot3
%       .Tuu3
%       .Tuu4
%       .Tuudot3
%       .Tuudot4
%    	.Tudotudot3         
%   	.Tudotudot4
%      	.time           computational time
%     
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - List of currently supported elements: 
%     TRI3
%
% Last modified: 21/10/2022, Mathieu Dubied, ETH Zürich

function tensors = tensors_hydro_FOM(myAssembly, elements, skinElements, skinElementFaces, vwater, rho)

t0=tic;
% data from myAssembly
nodes    = myAssembly.Mesh.nodes;                   % nodes table
nel      = myAssembly.Mesh.nElements;           	% number of elements
nnodes   = myAssembly.Mesh.nNodes;               	% number of nodes
freedofs = myAssembly.Mesh.EBC.unconstrainedDOFs;   % free DOFs

myMesh = myAssembly.Mesh;
mode = 'ELP';   % element level projection
u0 = zeros(myMesh.nDOFs, 1);    


% compute FOM tensors
disp(' FOM HYDRODYNAMIC TENSORS:')
fprintf(' Assembling %d elements ...', nel)

T1 = myAssembly.vector_skin('T1', 'weights', skinElements, skinElementFaces, vwater, rho);
T2 = myAssembly.matrix_skin('T2', 'weights', skinElements, skinElementFaces, vwater, rho);
Tu2 = myAssembly.matrix_skin('Tu2', 'weights', skinElements, skinElementFaces, vwater, rho);
Tu3 = myAssembly.tensor_skin('Tu3', [myMesh.nDOFs,myMesh.nDOFs,myMesh.nDOFs], 'weights', skinElements, skinElementFaces, vwater, rho);
Tudot2 = myAssembly.matrix_skin('Tudot2', 'weights', skinElements, skinElementFaces, vwater, rho);



time = toc(t0);

fprintf(' %.2f s (%.2f s)\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)


% outputs
tensors.T1 = T1;
tensors.T2 = T2;
tensors.Tu2 = Tu2;
tensors.Tu3 = Tu3;
tensors.Tudot2 = Tudot2;
% tensors.Tudot3 = Tudot3;
% tensors.Tuu3 = Tuu3;
% tensors.Tuu4 = Tuu4;
% tensors.Tuudot3 = Tuudot3;            
% tensors.Tuudot4 = Tuudot4;
% tensors.Tudotudot3 = Tudotudot3;            
% tensors.Tudotudot4 = Tudotudot4;
tensors.time = time;

end


