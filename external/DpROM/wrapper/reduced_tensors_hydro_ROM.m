% reduced_tensors_hydro_ROM
%
% Synthax:
% tensors = reduced_tensors_hydro_ROM(myAssembly, elements, V, skinElements, skinElementFaces, vwater, rho)
%
% Description: This function computes the reduced order hydrodynamic
% tensors at the Assembly level, by combining (i.e., summing), the
% element-level contributions. The obtained tensors are the ones used for
% ROM-n and/or ROM-d, where the shape variation is `fixed'.
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
%    	.Tru2            
%   	.Trudot2
%       .Truu3
%       .Truudot3
%       .Trudotudot3
%      	.time           computational time
%     
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - List of currently supported elements: 
%     TRI3, TET4
%
% Last modified: 12/03/2023, Mathieu Dubied, ETH Zurich

function tensors = reduced_tensors_hydro_ROM(myAssembly, elements, V, skinElements, skinElementFaces, vwater, rho)

t0=tic;

% data from myAssembly
nel      = myAssembly.Mesh.nElements;   % number of elements
myMesh = myAssembly.Mesh;
mode = 'ELP';                           % element level projection
m = size(V,2);                          % size of the ROM

% create ROM object
RomAssembly = ReducedAssembly(myMesh, V);

% compute reduced tensors
disp(' REDUCED HYDRODYNAMIC TENSORS:')
fprintf(' Assembling %d elements ...\n', nel)

tic;
Tr1 = RomAssembly.vector_skin('Te1', 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('   1st order terms - Tr1: %.2f s\n',toc)

tic;
Tru2 = RomAssembly.matrix_skin('Teu2', 'weights', skinElements, skinElementFaces, vwater, rho);
Trudot2 = RomAssembly.matrix_skin('Teudot2', 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('   2nd order terms - Tru2, Trudot2: %.2f s\n',toc)

tic;
Truu3 = 0.5*RomAssembly.tensor_skin('Teuu3',[m m m],[2 3],mode, 'weights', skinElements, skinElementFaces, vwater, rho);
Truudot3 = RomAssembly.tensor_skin('Teuudot3',[m m m],[2 3],mode, 'weights', skinElements, skinElementFaces, vwater, rho);
Trudotudot3 = 0.5*RomAssembly.tensor_skin('Teudotudot3',[m m m],[2 3],mode, 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('   3rd order terms - Truu3, Truudot3, Turudotudot3: %.2f s\n',toc)

% display time needed for computation
time = toc(t0);
fprintf(' TOTAL TIME: %.2f s\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)
fprintf(' SIZEs: %d \n\n', size(V,2))

% store outputs
tensors.Tr1 = Tr1;                     
tensors.Tru2 = Tru2;
tensors.Trudot2 = Trudot2;
tensors.Truu3 = Truu3;
tensors.Truudot3 = Truudot3;
tensors.Trudotudot3 = Trudotudot3;
tensors.time = time;

end


