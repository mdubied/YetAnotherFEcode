% PROM_spine_momentum_derivatives
%
% Synthax:
% der = PROM_spine_momentum_derivatives(q,qd,qdd,eta0,xi,tensors)
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
% (5) xi:               shape variation parameter vector/scalar
% (6) tensors:          tensors (struct) used to express the spine momentum
%                       change
%
% OUTPUTS:
% (1) der:              strucure array containing the partial derivatives 
%     
% Additional notes:
%   - q,qd, qdd should be understood as eta and dot{eta}, ddot{eta}.
%
% Last modified: 13/10/2023, Mathieu Dubied, ETH Zürich
function der = PROM_spine_momentum_derivatives(q,qd,qdd,eta0,xi,tensors)
    
    % access tensors
    T = tensors.T;
    TU3 = tensors.TU3;
    TU4 = tensors.TU4;
    TU34 = tensors.TU34;
    T = tensor(T);
    TU3 = tensor(TU3);
    TU4 = tensor(TU4);
    TU34 = tensor(TU34);
    
    % dfdq ________________________________________________________________
    dfdq = ttv(ttv(T,qdd,2),eta0+q,3) + ttv(ttv(T,qdd,2),eta0+q,2) + ...
           ttv(ttv(T,qd,2),qd,2) + ttv(ttv(T,qd,2),qd,3) + ...
           ttv(ttv(TU4,qd,2),xi,3) + ttv(ttv(TU3,qd,2),xi,2);
    dfdq = double(dfdq);
    
    % dfdqd _______________________________________________________________
    dfdqd = ttv(ttv(T,qd,3),eta0+q,3) + ttv(ttv(T,qd,2),eta0+q,3) + ...
            ttv(ttv(T,eta0+q,3),qd,3) + ttv(ttv(T,qd,2),eta0+q,2) + ...
            ttv(ttv(TU4,qd,3),xi,3) + ttv(ttv(TU4,qd,2),xi,3) + ...
            ttv(ttv(TU3,xi,3),qd,3) + ttv(ttv(TU3,qd,2),xi,2);
    dfdqd = double(dfdqd);

    % dfdqdd ______________________________________________________________
    dfdqdd = ttv(ttv(T,eta0+q,3),eta0+q,3) + ...
             ttv(ttv(TU4,eta0+q,3),xi,3) + ...
             ttv(ttv(TU3,xi,3),eta0+q,3) + ...
             ttv(ttv(TU34,xi,3),xi,3);
    dfdqdd = double(dfdqdd);

    % dfdp ________________________________________________________________
    dfdp = ttv(ttv(TU3,qdd,2),eta0+q,3) + ttv(ttv(TU34,qdd,2),xi,2) + ...
           ttv(ttv(TU4,qdd,2),eta0+q,2) + ttv(ttv(TU34,qdd,2),xi,2) + ...
           ttv(ttv(TU4,qd,2),qd,2) + ttv(ttv(TU3,qd,2),qd,3);
    dfdp = double(dfdp);
    
    % store results in output struct ______________________________________
    der.dfdq = dfdq;
    der.dfdqd = dfdqd;
    der.dfdqdd = dfdqdd;
    der.dfdp = dfdp;

end
