% solve_EoMs
%
% Synthax:
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax,varargin)
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
% (6) varargin:             additional name-value pairs if actuation is
%                           present. See below for the parser function.
%
% OUTPUTS:
% (1) TI_NL_PROM:           struct containing the solutions and related
%                           information
%     
% Last modified: 26/03/2023, Mathieu Dubied, ETH Zurich

function TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax,varargin)

    % PARSE ADDITIONAL INPUTS (ACTUATION) _________________________________
    [ACTUATION,tensors_topMuscle_PROM,tensors_bottomMuscle_PROM] = parse_inputs(varargin{:});

    % SIMULATION PARAMETERS AND ICs _______________________________________
    eta0 = zeros(size(V,2),1);
    etad0 = zeros(size(V,2),1);
    etadd0 = zeros(size(V,2),1);

    % HYDRODYNAMIC AND ACTUATION (EXTERNAL FORCES) ________________________
    if ACTUATION
        fprintf('Solving for a structure with actuation...\n')
        k = 1;
        B1TopMuscle = tensors_topMuscle_PROM.B1;
        B2TopMuscle = tensors_topMuscle_PROM.B2;
        B1BottomMuscle = tensors_bottomMuscle_PROM.B1;
        B2BottomMuscle = tensors_bottomMuscle_PROM.B2;
        % hydrodynamic and actuation forces
        F_ext = @(t,eta,etad) (double(tensors_hydro_PROM.Tr1) + ...
            double(tensors_hydro_PROM.Tru2*eta) + double(tensors_hydro_PROM.Trudot2*etad) + ...
            double(ttv(ttv(tensors_hydro_PROM.Truu3,eta,3), eta,2)) + ...
            double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad,3), eta,2)) + ...
            double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad,3), etad,2)) + ...
            k/2*(1-(1+0.04*sin(t*2*pi)))*(B1TopMuscle+B2TopMuscle*eta) + ...
            k/2*(1-(1-0.04*sin(t*2*pi)))*(B1BottomMuscle+B2BottomMuscle*eta)); % q, qd are reduced order DOFs

    else
        F_ext = @(t,eta,etad) (double(tensors_hydro_PROM.Tr1) + ...
        double(tensors_hydro_PROM.Tru2*eta) + ...
        double(tensors_hydro_PROM.Trudot2*etad) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truu3,eta,3), eta,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad,3), eta,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad,3), etad,2))); 
    end
    
    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate object
    TI_NL_PROM = ImplicitNewmark('timestep',h,'alpha',0.005);
    
    % modal nonlinear Residual evaluation function handle
    Residual_NL_red = @(eta,etad,etadd,t)residual_reduced_nonlinear_hydro(eta,etad,etadd,t,PROM_Assembly,F_ext,tensors_hydro_PROM);
    
    % time integration
    
    TI_NL_PROM.Integrate(eta0,etad0,etadd0,tmax,Residual_NL_red);
    TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

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