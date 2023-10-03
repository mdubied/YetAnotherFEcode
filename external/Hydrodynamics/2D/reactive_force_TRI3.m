% reactive_force_TRI3
%
% Synthax:
% force = reactive_force_TRI3(myAssembly, elements, skinElements, skinElementFaces, vwater, rho)
%
% Description: 
% This function computes the reactive force at the Assembly level,
% using the original nonlinear formulation (not in form of a polynomial).
%
%
% INPUTS:
%   - myAssembly: Assembly from YetAnotherFEcode.
%   -
%
% OUTPUT:
%   force: a struct variable with the following fields:
%       .f              hydrodynamic force at the Assembly level        
%      	.time           computational time
%     
%
% Additional notes:
%
% Last modified: 23/09/2023, Mathieu Dubied, ETH Zurich
function force = reactive_force_TRI3(Assembly, tailElement, tailNodeIndexInElement, mTilde, q, qd)
    u = Assembly.unconstrain_vector(q);
    ud = Assembly.unconstrain_vector(qd);
    force = Assembly.vector_hydro2('reactive_force_full', 'weights', tailElement, tailNodeIndexInElement, mTilde,u,ud); 
    force = Assembly.constrain_vector(force);
end
