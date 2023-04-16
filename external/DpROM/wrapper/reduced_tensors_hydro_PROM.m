% reduced_tensors_hydro_PROM
%
% Synthax:
% tensors = reduced_tensors_hydro_ROM(myAssembly, elements, V, U, fourthOrder, skinElements, skinElementFaces, vwater, rho, c)
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
%   - U: Shape variation basis
%   - fourthOrder: 1 or 0, for the computation (1) of 4th order tensors
%   - skinElements
%   - skinElementFaces
%   - vwater
%   - rho:
%   - c: scaling factor for the hydrodynamic thrust force
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
%   - 4th order tensors not using the scaling factor c for now
%
% Last modified: 14/04/2022, Mathieu Dubied, ETH Zurich

function tensors = reduced_tensors_hydro_PROM(myAssembly, elements, V, U, fourthOrder, skinElements, skinElementFaces, vwater, rho, c)

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
fprintf('   0th order in ud:\n')

tic;
Tr1 = RomAssembly.vector_skin('Te1', 'weights', skinElements, skinElementFaces, vwater, rho, c);
fprintf('       1st order terms - Tr1: %.2f s\n',toc)

tic;
Tru2 = RomAssembly.matrix_skin('Teu2', 'weights', skinElements, skinElementFaces, vwater, rho, c);
Trudot2 = RomAssembly.matrix_skin('Teudot2', 'weights', skinElements, skinElementFaces, vwater, rho, c);
fprintf('       2nd order terms - Tru2, Trudot2: %.2f s\n',toc)

tic;
Truu3 = 0.5*RomAssembly.tensor_skin('Teuu3', [m m m], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
Truudot3 = RomAssembly.tensor_skin('Teuudot3', [m m m], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
Trudotudot3 = 0.5*RomAssembly.tensor_skin('Teudotudot3', [m m m], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
fprintf('       3rd order terms - Truu3, Truudot3, Trudotudot3: %.2f s\n',toc)


fprintf('   1st order in ud:\n')

tic;
Tr2 = RomAssembly.matrix_skin_PROM('Te2', U,'weights', skinElements, skinElementFaces, vwater, rho, c);
fprintf('       2nd order terms - Tr2: %.2f s\n',toc)

tic;
Tru3 = RomAssembly.tensor_skin_PROM('Teu3', U, [m m md], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
Trudot3 = RomAssembly.tensor_skin_PROM('Teudot3', U, [m m md], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
fprintf('       3rd order terms - Tru3, Trudot3: %.2f s\n',toc)

if fourthOrder
    tic;
    Truu4 = 0.5*RomAssembly.tensor4_skin_PROM('Teuu4', U, [m m m md], [2 3 4], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
    Truudot4 = RomAssembly.tensor4_skin_PROM('Teuudot4', U,  [m m m md], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
    Trudotudot4 = 0.5*RomAssembly.tensor4_skin_PROM('Teudotudot4', U, [m m m md], [2 3], mode, 'weights', skinElements, skinElementFaces, vwater, rho, c);
    fprintf('       4th order terms - Truu4, Truudot4, Trudotudot4: %.2f s\n',toc)
else
    fprintf('       4th order terms - not computed \n')
end



% display total time needed for computation
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

if fourthOrder
    tensors.Truu4 = Truu4;
    tensors.Truudot4 = Truudot4;
    tensors.Trudotudot4 = Trudotudot4;
end

tensors.time = time;

end


