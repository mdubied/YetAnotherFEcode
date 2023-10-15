% compute_drag_tensors_ROM
%
% Synthax:
% tensors = compute_drag_tensors_ROM(ROMAssembly, skinElements, skinElementFaces, rho)
%
% Description: This function computes the reduced order drag
% tensors at the global level, by combining (i.e., summing), the
% element-level contributions. The obtained tensors are the ones used in
% the ROM.
%
% INPUTS
%   - ROMAssembly:          Reduced assembly from YetAnotherFEcode
%   - skinElements:         logical array of length nElements, with 1 if an
%                           element is a skin element, and 0 else
%   - skinElementFaces:     matrix of size nElements x 2, giving for each
%                           skin element which face (1,2 or 3) is part of
%                           the skin. If two faces are part of the skin, we
%                           use the 2 columns of the matrix.
%   - rho:                  water density
%
% OUTPUTS
%   tensors: a struct variable with the following fields:     
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
%   - List of currently supported elements: TRI3
%
% Last modified: 11/10/2023, Mathieu Dubied, ETH Zurich

function tensors = compute_drag_tensors_ROM(ROMAssembly, skinElements, skinElementFaces, rho) 

    t0=tic;
    
    % data from ROM Assembly
    nel = ROMAssembly.Mesh.nElements;      % number of elements
    V = ROMAssembly.V;
    mode = 'ELP';                           % element level projection
    m = size(V,2);                          % size of the ROM
    
    % compute reduced tensors
    disp(' REDUCED HYDRODYNAMIC TENSORS:')
    fprintf(' Assembling %d elements ...\n', nel)

    tic;
    Tr1 = ROMAssembly.vector_skin('Te1', 'weights', skinElements, skinElementFaces, rho);
    fprintf('   1st order term - Tr1: %.2f s\n',toc)
    
    tic;
    Tru2 = ROMAssembly.matrix_skin('Teu2', 'weights', skinElements, skinElementFaces, rho);
    Trudot2 = ROMAssembly.matrix_skin('Teudot2', 'weights', skinElements, skinElementFaces, rho);
    fprintf('   2nd order terms - Tru2, Trudot2: %.2f s\n',toc)
    
    tic;
    Truu3 = 0.5*ROMAssembly.tensor_skin('Teuu3',[m m m],[2 3],mode, 'weights', skinElements, skinElementFaces, rho);
    Truudot3 = ROMAssembly.tensor_skin('Teuudot3',[m m m],[2 3],mode, 'weights', skinElements, skinElementFaces, rho);
    Trudotudot3 = 0.5*ROMAssembly.tensor_skin('Teudotudot3',[m m m],[2 3],mode, 'weights', skinElements, skinElementFaces, rho);
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


