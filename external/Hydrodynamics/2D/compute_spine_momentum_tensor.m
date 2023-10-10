% compute_spine_momentum_tensor
% WARNING: does not work properly due to sparsify function
%
% Synthax:
% tensor = compute_spine_momentum_tensor(ROMAssembly, spineElementWeights, normalisationFactors, nodeIdxPosInElements, mTilde)
%
% Description: 
% This function computes the 4th order tensor that can be used to compute
% the change of momentum along the fish spine at the Assembly level (FOM)
%
% INPUTS
%   - Assembly:             Assembly from YetAnotherFEcode.
%   - spineElementWeights:  array of length nElements, with 1 if element is
%                           a spine element, 0 else
%   - normalisationFactors: array of size nElements with zeros everywhere,
%                           except at the indexes of the spine elements.
%                           For these indexes, the inverse of the distance
%                           between the two spine nodes is entailed.
%   - nodeIdxPosInelements: matrix of size nElements x 2. For the row
%                           corresponding to spine elements, it gives the
%                           column index (1 to 3) in the element of the 
%                           node close to the tail (1st column) and the 
%                           node close to the head (2nd column).
%   - mTilde:               virtual mass linear density
%
%
% OUTPUTS
%   tensor: the 4th order tensor mentioned above
%       
%
% Additional notes:  WARNING: does not work properly due to sparsify function
%
% Last modified: 10/10/2023, Mathieu Dubied, ETH Zurich
function tensor = compute_spine_momentum_tensor(ROMAssembly, spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde)
    m = size(ROMAssembly.V,2);
    tensor = ROMAssembly.tensor_spine_momentum('spine_momentum_tensor',[m m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde); 
end
