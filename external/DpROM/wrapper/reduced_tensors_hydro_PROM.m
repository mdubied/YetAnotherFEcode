% reduced_tensors_hydro_PROM
%
% Synthax:
% tensors = reduced_tensors_hydro_ROM(myAssembly, elements, V, skinElements, skinElementFaces, vwater, rho)
%
% Description: This function computes the reduced order hydrodynamic
% tensors at the Assembly level, by combining (i.e., summing), the
% element-level contributions. The obtained tensors are the ones used for
% the Parametric ROM (PROM).
%
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - elements: table of the elements
%   - V: Reduced Order Basis (unconstrained)
%   - skinElements
%   - skinElementFaces
%   - vwater
%   - rho:
% OUTPUT:
%   tensors: a struct variable with the following fields*:
%       .Tr1
%       .Tr2            specific to PROM
%    	.Tru2
%       .Tru3           specific to PROM
%   	.Trudot2
%   	.Trudot3        specific to PROM
%       .Truu3
%   	.Truu4          specific to PROM
%       .Truudot3
%   	.Truudot4       specific to PROM
%       .Trudotudot3
%   	.Trudotudot4    specific to PROM
%      	.time           computational time
%     
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - List of currently supported elements: 
%     TRI3
%
% Last modified: 17/11/2022, Mathieu Dubied, ETH Zürich

function tensors = reduced_tensors_hydro_PROM(myAssembly, elements, V, U, skinElements, skinElementFaces, vwater, rho)

t0=tic;

% data from myAssembly
nel      = myAssembly.Mesh.nElements;   % number of elements
myMesh = myAssembly.Mesh;
mode = 'ELP';                           % element level projection
m = size(V,2);                          % size of the ROM
md = size(U,2);                         % size of the defect basis

% create ROM object
RomAssembly = ReducedAssembly(myMesh, V);

% compute reduced tensors
disp(' REDUCED HYDRODYNAMIC TENSORS (PROM):')
fprintf(' Assembling %d elements ...\n', nel)

tic;
Tr1 = RomAssembly.vector_skin('Te1', 'weights', skinElements, skinElementFaces, vwater, rho);
Tr2 = RomAssembly.matrix_skin_PROM('Te2', U,'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('   1st order terms in u - Tr1, Tr2: %.2f s\n',toc)

tic;
Tru2 = RomAssembly.matrix_skin('Teu2', 'weights', skinElements, skinElementFaces, vwater, rho);
Tru3 = RomAssembly.tensor_skin_PROM('Teu3', U, [m m md], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho);
Trudot2 = RomAssembly.matrix_skin('Teudot2', 'weights', skinElements, skinElementFaces, vwater, rho);
Trudot3 = RomAssembly.tensor_skin_PROM('Teudot3', U, [m m md], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('   2nd order terms in u - Tru2, Tru3, Trudot2, Trudot3: %.2f s\n',toc)

tic;
Truu3 = 0.5*RomAssembly.tensor_skin('Teuu3', [m m m], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho);
Truudot3 = RomAssembly.tensor_skin('Teuudot3', [m m m], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho);
Trudotudot3 = 0.5*RomAssembly.tensor_skin('Teudotudot3', [m m m], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('   3rd order terms in u - Truu3, Truudot3, Trudotudot3: %.2f s\n',toc)

% display time needed for computation
time = toc(t0);
fprintf(' TOTAL TIME: %.2f s\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)
fprintf(' SIZEs: %d \n\n', size(V,2))

% store outputs
tensors.Tr1 = Tr1;
tensors.Tr2 = Tr2; 
tensors.Tru2 = Tru2;
tensors.Tru3 = Tru3;
tensors.Trudot2 = Trudot2;
tensors.Trudot3 = Trudot3;
tensors.Truu3 = Truu3;
tensors.Truudot3 = Truudot3;
tensors.Trudotudot3 = Trudotudot3;
tensors.time = time;

end


