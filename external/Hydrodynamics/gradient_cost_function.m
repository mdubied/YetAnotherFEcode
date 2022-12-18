% gradient_cost_function
%
% Synthax:
% nabla_Lr = gradient_cost_function(q,qd,xi,tensors_DpROM)
%
% Description: Implementation of the optimization pipeline 1 presented in
% the paper, based on Newton's method. 
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qdot:             vector of time domain velocities
% (3) tensors_DpROM:    structure array containing reduced tensors of
%                       F_hydro model
%
% OUTPUTS:
% (1) xi_star:          optimal shape parameter(s)
%     
%
% Additional notes: -
%
% Last modified: 17/12/2022, Mathieu Dubied, ETH Zürich
function nablaLr = gradient_cost_function(dr,xi,eta,etad,s,sd,tensors_hydro_PROM)
    N = size(eta,2);
    nablaLr = zeros(size(xi,1),1);
    
    for i=1:N 
        derivative_hydro = DpROM_hydro_derivatives(eta(:,i),etad(:,i),xi,tensors_hydro_PROM);
        dfhydrodeta = derivative_hydro.dfdq; 
        dfhydrodetad = derivative_hydro.dfdqd;
        dfhydrodxi = derivative_hydro.dfdp;
        dfdxi_i = dfhydrodxi + dfhydrodeta*s(:,i) + dfhydrodetad*sd(:,i);
                
        nablaLr = nablaLr - dr.'*dfdxi_i;
    end  

    
end