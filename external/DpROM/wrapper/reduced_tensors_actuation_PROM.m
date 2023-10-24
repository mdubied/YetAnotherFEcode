% reduced_tensors_actuation_PROM
%
% Synthax:
% tensors = reduced_tensors_actuation_PROM(myAssembly, elements, V, actuationElements, actuationDirection)
%
% Description: This function computes the reduced order actuation
% vector and matrices at the Assembly level, by combining (i.e., summing), 
% the element-level contributions. The obtained expressions are the ones 
% used for the Parametric ROM (PROM)..
%
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - V: Reduced Order Basis (unconstrained)
%   - actuationElements: elements subject to actuation
%   - actuationDirection: a global direction in which the actuation takes
%                         place
%
% OUTPUT:
%   tensors: a struct variable with the following fields:
%       .B1                       
%    	.B2 
%    	.B3 
%      	.time           computational time
%     
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - List of currently supported elements: 
%     TRI3
%
% Last modified: 06/03/2023, Mathieu Dubied, ETH Zürich

function tensors = reduced_tensors_actuation_PROM(myAssembly, V, U, actuationElements, actuationDirection)

t0=tic;

% data from myAssembly
nel      = myAssembly.Mesh.nElements;   % number of elements
myMesh = myAssembly.Mesh;

% create ROM object
RomAssembly = ReducedAssembly(myMesh, V);

% compute reduced tensors
disp(' REDUCED ACTUATION TENSORS (PROM):')
fprintf(' Assembling %d elements ...\n', nel)

tic;
B1 = RomAssembly.vector_actuation('B1', 'weights', actuationElements, actuationDirection);
fprintf('   B1: %.2f s\n',toc)

tic;
B2 = RomAssembly.matrix_actuation('B2', 'weights', actuationElements, actuationDirection);
fprintf('   B2: %.2f s\n',toc)

tic;
B3 = RomAssembly.matrix_actuation_PROM('B3', U, 'weights', actuationElements, actuationDirection);
fprintf('   B3: %.2f s\n',toc)


% display time needed for computation
time = toc(t0);
fprintf(' TOTAL TIME: %.2f s\n',toc(t0))
% fprintf(' SPEED: %.1f el/s\n',nel/time)
% fprintf(' SIZEs: %d \n\n', size(V,2))

% store outputs
tensors.B1 = B1; 
tensors.B2 = B2;  
tensors.B3 = B3; 
tensors.time = time;

end
