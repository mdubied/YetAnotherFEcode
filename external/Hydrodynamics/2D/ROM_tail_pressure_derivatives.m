% ROM_tail_pressure_derivates
%
% Synthax:
% der = ROM_tail_pressure_derivatives(q,qd,A,B)
%
% Description: This function returns the partial derivatives of the change
% in tail pressure force, which is understood as a force acting on the spine.
% The partial derivatives are needed to solve the sensitivity ODE as well 
% as the EoMs.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qd:               vector of time domain velocities
% (3) A:                matrix used to express the tail pressure force
% (4) B:                matrix used to express the tail pressure force
% (5) R:                rotation matrix
% (6) mTilde:           virtual mass linear density
% (7) x0:               initial node positions of the tail node
% (8) VTail:            part of the ROB matrix corresponding to the tail
%                       element
% OUTPUTS:
% (1) der:              struct containing the partial derivatives
%     
%
% Additional notes:
%   - q and qd should be understood as eta and dot{eta} (ROM).
%   - x0 should be understood as eta_0 (ROM)
%   - p should be understood as xi.
%
% Last modified: 10/10/2023, Mathieu Dubied, ETH Zürich
function der = ROM_tail_pressure_derivatives(q,qd,A,B,R,mTilde,w,x0,VTail)

    % dfdq ________________________________________________________________
    firstTerm1 = A*VTail*qd;                                % vector
    firstTerm2 = dot(A*VTail*qd,R*B*VTail*(x0+q));          % scalar
    firstTerm3 = B*VTail*(x0+q);                            % vector
    firstTerm1RB = firstTerm1.'*R*B*VTail;                  % row vector 
    
    firstTerm = 2*firstTerm2*firstTerm3*firstTerm1RB;     % outer-product between the last two terms to get a matrix
    secondTerm = dot(A*VTail*qd,R*B*VTail*(x0+q))^2*B*VTail;
    
    dfdq = 0.5*mTilde*w^3*VTail.'*(firstTerm + secondTerm);
    
    % dfdqd _______________________________________________________________
    term1RB = R*B*VTail*(x0+q);                                  % vector 
    term1B = (A*VTail).'*term1RB;                             % row vector
    term2 = dot(A*VTail*qd,R*B*VTail*(x0+q));               % scalar
    term3 = B*VTail*(x0+q);                               % vector 
    dfdqd = 0.5*mTilde*w^3*VTail.'*term2*term3*term1B.';    % with outer-product
    
    % store results in output struct ______________________________________
    der.dfdq = dfdq;
    der.dfdqd = dfdqd;

end
