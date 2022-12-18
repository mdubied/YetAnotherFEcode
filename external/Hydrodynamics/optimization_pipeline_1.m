% optimization_pipeline_1
%
% Synthax:
% xi_star = optimization_pipeline_1(q,qd,xi,tensors_DpROM)
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
% Additional notes:
%
% Last modified: 17/12/2022, Mathieu Dubied, ETH Zürich
function [xiStar,LrEvo] = optimization_pipeline_1(V,d,tensors_PROM, tensors_hydro_PROM,eta,etad,s,sd)
    xi_k = 0;
    N = size(eta,2);
    dr = reduce_vector(d,V);
    Lr = cost_function(N,tensors_hydro_PROM,eta,etad,dr);
    LrEvo = Lr;
    for i = 1:3
        disp(i)
        nablaLr = gradient_cost_function(dr,xi_k,eta,etad,s,sd,tensors_hydro_PROM);
        LrEvo = [LrEvo, cost_function(N,tensors_hydro_PROM,eta,etad,dr)];
        xi_k = xi_k - 0.8*nablaLr;
    end
    xiStar = xi_k;

end


function dr = reduce_vector(d,V)
    n = size(V,1);
    m = size(V,2);
    dr = zeros(m,1);
    for i=1:6:n
        Ve = V(i:i+5,:);
        de = [d;d;d];
        dr = dr + Ve.'*de;  
    end
end


function Lr = cost_function(N,tensors_hydro_PROM,eta,etad,dr)
    Lr = 0;
    for i=1:N-2
        eta_i = eta(:,i);
        etad_i = etad(:,i);
        fhydro = double(tensors_hydro_PROM.Tr1) + ...
        double(tensors_hydro_PROM.Tru2*eta_i) + double(tensors_hydro_PROM.Trudot2*etad_i) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truu3,eta_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad_i,3), etad_i,2));
       
        Lr = Lr + dr'*fhydro;
    end
end
