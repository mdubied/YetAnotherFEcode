% gradient_cost_function
%
% Synthax:
% nabla_Lr = gradient_cost_function(dr,xi,eta,etad,s,sd,tensors_hydro_PROM,FOURTHORDER)
%
% Description:  gradient of the (reduced) cost function Lr. The gradient is
%               analytical and based on the hydrodynamic tensors
%
% INPUTS: 
% (1) dr:                   reduced forward swimming direction vector
% (2) xi:                   vector of shape variation parameters
% (3) eta:                  solution for the reduced state variables
% (4) etad:                 solution for the reduced state derivatives
% (5) s:                    solution for the sensitivity
% (6) sd:                   solution for the sensitivity derivative
% (7) tensors_hydro_PROM:   reduced tensors for the hydrodynamic forces
% (8) FOURTHORDER:          logical value for 4th order tensors (1), or not
%                           (0)
%                       
% OUTPUTS:
% (1) nablaLr:          gradient of the (reduced cost function)
%     
%
% Last modified: 15/03/2023, Mathieu Dubied, ETH Zürich

function nablaLr = gradient_cost_function(dr,xi,eta,etad,s,sd,tensors_hydro_PROM,FOURTHORDER)
    N = size(eta,2);
    nablaLr = zeros(size(xi,1),1);
    
    for i=1:N 
        derivative_hydro = DpROM_hydro_derivatives(eta(:,i),etad(:,i),xi,tensors_hydro_PROM,FOURTHORDER);
        dfhydrodeta = derivative_hydro.dfdq; 
        dfhydrodetad = derivative_hydro.dfdqd;
        dfhydrodxi = derivative_hydro.dfdp;
        dfdxi_i = dfhydrodxi + dfhydrodeta*s(:,i) + dfhydrodetad*sd(:,i);
                
        nablaLr = nablaLr - dr.'*dfdxi_i;
    end  

    
end