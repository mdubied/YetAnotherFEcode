% solve_sensitivities
%
% Synthax:
% TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM,tensors_hydro_PROM,etaSol,etadSol,etaddSol,h,tmax,varargin)
%
% Description: Computes the solutions for the sensitivity S for [0,tmax] 
% for a given PROM assembly, xi_k, and time step h
%
% INPUTS: 
% (1) V:                    ROB 
% (2) xi_k:                 current shape variation xi to consider
% (3) PROM_Assembly:        PROM assembly    
% (4) tensors_PROM:         (reduced) tensors for the internal forces  
% (5) tensors_hydro_PROM:   (reduced) tensors for the hydrdynamic forces
% (6) etaSol:               solution for eta to consider
% (7) etadSol:              solution for dot{eta} to consider
% (8) etadSol:              solution for ddot{eta} to consider
% (9) h:                    time step for time integration
% (10) tmax:                simulation for [0,tmax]
%
% OUTPUTS:
% (1) TI_sens:              struct containing the solutions and related
%                           information
%     
%
% Last modified: 26/03/2023, Mathieu Dubied, ETH Zurich

function TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM,tensors_hydro_PROM,etaSol,etadSol,etaddSol,h,tmax,FOURTHORDER,varargin)        

    % PARSE ADDITIONAL INPUTS (ACTUATION) _________________________________
    [ACTUATION,tensors_topMuscle_PROM,tensors_bottomMuscle_PROM] = parse_inputs(varargin{:});

    % SIMULATION PARAMETERS AND ICs _______________________________________
    s0 = zeros(size(V,2),size(xi_k,1));
    sd0 = zeros(size(V,2),size(xi_k,1));
    sdd0 = zeros(size(V,2),size(xi_k,1));
    
    % EVALUATE PARTIAL DERIVATIVE ALONG NOMINAL SOLUTION __________________
    secondOrderDer = 0;
    pd_hydro_PROM = @(eta,etad)DpROM_hydro_derivatives(eta,etad,xi_k,tensors_hydro_PROM,FOURTHORDER,secondOrderDer);
    pd_fint_PROM = @(eta)DpROM_derivatives(eta,tensors_PROM);
    
    if ACTUATION
        pd_factTop_PROM = @(a)PROM_actu_derivatives(tensors_topMuscle_PROM,a);
        pd_factBottom_PROM = @(a)PROM_actu_derivatives(tensors_bottomMuscle_PROM,a);
    end

    % TIME INTEGRATION ____________________________________________________
    
    %instantiate object for time integration
    TI_sens = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true,'sens',true);
    
    %residual function handle
    if ACTUATION
         Residual_sens = @(s,sd,sdd,t)residual_linear_sens(s,sd,sdd,t,PROM_Assembly,etaSol,etadSol,etaddSol,pd_hydro_PROM,pd_fint_PROM,h,...
             'ACTUATION',ACTUATION,'pd_factTop',pd_factTop_PROM,'pd_factBottom',pd_factBottom_PROM);
    else
        Residual_sens = @(s,sd,sdd,t)residual_linear_sens(s,sd,sdd,t,PROM_Assembly,etaSol,etadSol,etaddSol,pd_hydro_PROM,pd_fint_PROM,h);
    end
    
    
    %time integration
    TI_sens.Integrate(s0,sd0,sdd0,tmax,Residual_sens);
    %TI_sens.Solution.s = V * TI_sens.Solution.q; % get full order solution

end

% parse input
function [ACTUATION,topMuscle,bottomMuscle] = parse_inputs(varargin)
defaultACTUATION = 0;
defaultTopMuscle = 0;
defaultBottomMuscle = 0;
p = inputParser;

addParameter(p,'ACTUATION',defaultACTUATION,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'topMuscle',defaultTopMuscle );
addParameter(p,'bottomMuscle',defaultBottomMuscle );

parse(p,varargin{:});

ACTUATION = p.Results.ACTUATION;
topMuscle = p.Results.topMuscle;
bottomMuscle = p.Results.bottomMuscle;

end
