% reduced_cost_function_w_constraints
%
% Synthax:
% Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,dr,A,b,barrierParam)
%
% Description:  Computes the cost function value in the ROB, considering
%               constraints on the shape variation parameters xi. Upper and
%               lower bounds inequality are considered. These inequality
%               constraints are included as (log) barrier functions. 
%
% INPUTS: 
% (1) N:                number of time steps             
% (2) tensors_hydro:    reduced tensors for the hydrodynamic forces
% (3) eta:              solution for the reduced state variables
% (4) etad:             solution for the reduced state derivatives
% (5) dr:               reduced forward swimming direction vector
% (6)-(7) A, b:         constraints on xi of the form Axi<b  
% (8) barrierParam:     parameter to scale (1/barrierParam) the barrier functions      
%                   
%
% OUTPUTS:
% (1) Lr:   reduced cost function value    
%     
%
% Last modified: 31/03/2023, Mathieu Dubied, ETH Zurich

function Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xi,dr,A,b,barrierParam)
    Lr = 0;
    nConstraints = size(b);
   
    LwoB = 0;
    
    for t=1:N-2
        eta_i = eta(:,t);
        etad_i = etad(:,t);

        % hydrodynamic forces needed for cost function
        fhydro = double(tensors_hydro_PROM.Tr1) + ...
        double(tensors_hydro_PROM.Tru2*eta_i) + double(tensors_hydro_PROM.Trudot2*etad_i) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truu3,eta_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad_i,3), etad_i,2));

        % constraints (log barriers) to be included in the cost function
        logBarrierInTimeStep = 0;
        if nConstraints ~= 0
            for i = 1:nConstraints 
                logBarrierInTimeStep = logBarrierInTimeStep - 1/barrierParam*log(-(A(i,:)*xi-b(i)));
            end
        end

        % final cost function at time step t
        Lr = Lr - dr'*fhydro + logBarrierInTimeStep;
        LwoB = LwoB -dr'*fhydro;
    end

    % print cost function without part stemming from barrier functions
    fprintf('Drag minimization: Lr = %.4f\n',LwoB)
end