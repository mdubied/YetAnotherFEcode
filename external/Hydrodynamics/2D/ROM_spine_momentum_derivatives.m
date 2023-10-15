% ROM_spine_momentum_derivatives
%
% Synthax:
% der = ROM_spine_momentum_derivatives(q,qd,qdd,x0,A,B)
%
% Description: This function returns the partial derivatives of the change
% in the fish momentum, which is understood as a force acting on the spine.
% The partial derivatives are needed to solve the sensitivity ODE as well 
% as the EoMs.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qd:               vector of time domain velocities
% (3) qdd:              vector of time domain acceleration
% (4) x0:               node position projected in ROM
% (5) T:                tensor used to express the spine momentum change
%
% OUTPUTS:
% (1) der:              strucure array containing the partial derivatives 
%     
% Additional notes:
%   - q,qd, qdd should be understood as eta and dot{eta}, ddot{eta}.
%
% Last modified: 13/10/2023, Mathieu Dubied, ETH Zürich
function der = ROM_spine_momentum_derivatives(q,qd,qdd,x0,T)
    T = tensor(T);
    
    % dfdq ________________________________________________________________
    dfdq = ttv(ttv(T,qdd,2),x0+q,3) + ttv(ttv(T,qdd,2),x0+q,2) + ...
           ttv(ttv(T,qd,2),qd,2) + ttv(ttv(T,qd,2),qd,3);
    dfdq = double(dfdq);
    
    % dfdqd _______________________________________________________________
    dfdqd = ttv(ttv(T,qd,3),x0+q,3) + ttv(ttv(T,qd,2),x0+q,3) + ...
            ttv(ttv(T,x0+q,3),qd,3) + ttv(ttv(T,qd,2),x0+q,3);
    dfdqd = double(dfdqd);

    % dfdqdd ______________________________________________________________
    dfdqdd = ttv(ttv(T,x0+q,3),x0+q,3);
    dfdqdd = double(dfdqdd);
    
    % store results in output struct ______________________________________
    der.dfdq = dfdq;
    der.dfdqd = dfdqd;
    der.dfdqdd = dfdqdd;

end
