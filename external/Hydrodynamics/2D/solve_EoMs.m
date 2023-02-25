% solve_EoMs
%
% Synthax:
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax)
%
% Description: Computes the solutions for eta and dot{eta} for [0,tmax] for 
% a given PROM assembly and time step
%
% INPUTS: 
% (1) V:                    ROB   
% (2) PROM_Assembly:        PROM assembly           
% (3) tensors_hydro_PROM:   (reduced) tensors for the hydrdynamic forces
% (4) h:                    time step for time integration
% (5) tmax:                 simulation for [0,tmax]
%
% OUTPUTS:
% (1) TI_NL_PROM:           struct containing the solutions and related
%                           information
%     
%
% Additional notes:
%
% Last modified: 17/12/2022, Mathieu Dubied, ETH Zürich

function TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax)
    % SIMULATION PARAMETERS AND ICs _______________________________________
    eta0 = zeros(size(V,2),1);
    etad0 = zeros(size(V,2),1);
    etadd0 = zeros(size(V,2),1);

    % HYDRODYNAMIC (EXTERNAL FORCES) ______________________________________
    F_ext = @(t,eta,etad) (double(tensors_hydro_PROM.Tr1) + ...
        double(tensors_hydro_PROM.Tru2*eta) + ...
        double(tensors_hydro_PROM.Trudot2*etad) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truu3,eta,3), eta,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad,3), eta,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad,3), etad,2))); 
    
    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate object
    TI_NL_PROM = ImplicitNewmark('timestep',h,'alpha',0.005);
    
    % modal nonlinear Residual evaluation function handle
    Residual_NL_red = @(eta,etad,etadd,t)residual_reduced_nonlinear_hydro(eta,etad,etadd,t,PROM_Assembly,F_ext);
    
    % time integration
    
    TI_NL_PROM.Integrate(eta0,etad0,etadd0,tmax,Residual_NL_red);
    TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

end 