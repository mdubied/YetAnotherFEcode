% modal_derivatives
%
% Synthax:
% [MD, names] = modal_derivatives(myAssembly, elements, Phi)
%
% Description: compute Modal Derivatives.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode. It MUST contain the field
%     ".DATA.K" storing the unconstrained linear stiffness matrix K0
%   - elements: table of the elements
%   - Phi: tall matrix containing by columns the selected Vibration Modes
%     (unconstrained)
% OUTPUTS
%   - MD: tall matrix containing by columns all the Modal Derivatives that
%     can be computed from the given Phi (unconstrained MDs are returned).
%   - names: matrix containing the subscripts of the MDs, for convenience.
%
% Additional notes:
%   - this function uses the stiffness_matrix_derivative function 
%     implemented in the Julia module "DpROM.jl". 
%   - as such, this function supports ONLY models meshed with the elements
%     supported by both YetAnotherFEcode AND the DpROM.jl
%   - List of currently supported elements: 
%     Q8, TET10, HEX20, WED15               (in YetAnotherFEcode)
%     Q8, TET10, HEX20, WED15, Q4, HEX8     (in DpROM.jl)
%
% Created: 14 May 2021
% Author: Jacopo Marconi, Politecnico di Milano


function [MD, names] = modal_derivatives(myAssembly, elements, Phi)

n = size(Phi,1);
n_VMs = size(Phi,2);

K0 = myAssembly.DATA.K;
K0 = myAssembly.constrain_matrix( K0 );

MD = zeros(n, n_VMs*(n_VMs+1)/2);
names = zeros(n_VMs*(n_VMs+1)/2, 2);
kk = 1;
for jj = 1 : n_VMs
    
    Phi_j = Phi(:, jj);
    dK_deta_j = stiffness_matrix_derivative(myAssembly, elements, Phi_j);
    dK_deta_j = myAssembly.constrain_matrix( dK_deta_j );
    
    for ii = 1 : n_VMs
        if ii < jj
            continue
        end
        
        Phi_i = myAssembly.constrain_vector( Phi(:, ii) );
        dPhi_i_deta_j = -K0\(dK_deta_j * Phi_i); 
        
        th =  dPhi_i_deta_j / max(abs(dPhi_i_deta_j));
        MD(:,kk) = myAssembly.unconstrain_vector( th );
        names(kk, :) = [ii jj];
        kk = kk + 1;
    end
end
disp(' ')
