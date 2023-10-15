% compute_spine_momentum_tensors_PROM
%
% Synthax:
% tensors = compute_spine_momentum_tensors_PROM(PROMAssembly, spineElementWeights, normalisationFactors, nodeIdxPosInElements, mTilde)
%
% Description: 
% This function computes the 4th order tensors that can be used to compute
% the change of momentum along the fish spine at the global level of the
% PROM. In particular, it also computes the tensors needed to compute the
% partial derivatives in the PROM (U3, U4, U34, where the number indicates
% the placement of the U matrix multiplications with the tensor T).
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
% Additional notes:  
%
% Last modified: 15/10/2023, Mathieu Dubied, ETH Zurich
function tensors = compute_spine_momentum_tensors_PROM(PROMAssembly, spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde)
    m = size(PROMAssembly.V,2);
    md = size(PROMAssembly.U,2);
    Tr = PROMAssembly.tensor_spine_momentum('spine_momentum_tensor',[m m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde); 
    TrU3 = PROMAssembly.tensor_spine_momentum_U3('spine_momentum_tensor',[m m md m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde); 
    TrU4 = PROMAssembly.tensor_spine_momentum_U4('spine_momentum_tensor',[m m m md], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde); 
    TrU34 = PROMAssembly.tensor_spine_momentum_U34('spine_momentum_tensor',[m m md md], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde); 

    tensors.T = Tr;
    tensors.TU3 = TrU3;
    tensors.TU4 = TrU4;
    tensors.TU34 = TrU34;
end
