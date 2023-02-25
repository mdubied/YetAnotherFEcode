% reduced_cost_function
%
% Synthax:
% Lr = reduced_cost_function(N,tensors_hydro_PROM,eta,etad,dr)
%
% Description: Computes the cost function value in the ROB
%
% INPUTS: 
% (1) N:                number of time steps             
% (2) tensors_hydro:    reduced tensors for the hydrodynamic forces
% (3) eta:              solution for the reduced state variables
% (4) etad:             solution for the reduced state derivatives
% (5) dr:               reduced forward swimming direction vector 
%
%
% OUTPUTS:
% (1) Lr:   reduced cost function value    
%     
%
% Additional notes: none
%
% Last modified: 17/12/2022, Mathieu Dubied, ETH Zürich

function Lr = reduced_cost_function(N,tensors_hydro_PROM,eta,etad,dr)
    Lr = 0;
    for i=1:N-2
        eta_i = eta(:,i);
        etad_i = etad(:,i);
        fhydro = double(tensors_hydro_PROM.Tr1) + ...
        double(tensors_hydro_PROM.Tru2*eta_i) + double(tensors_hydro_PROM.Trudot2*etad_i) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truu3,eta_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad_i,3), etad_i,2));
       
        Lr = Lr - dr'*fhydro;
    end
end