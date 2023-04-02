% solve_EoMs_FOM
%
% Synthax:
% TI_NL_PROM = solve_EoMs_FOM(FOM_Assembly,tensors_hydro,h,tmax,varargin)
%
% Description: Computes the solutions for u and dot{u} for [0,tmax] for 
% a given FOM assembly and time step
%
% INPUTS:   
% (1) PROM_Assembly:        FOM assembly           
% (2) tensors_hydro_PROM:   (reduced) tensors for the hydrdynamic forces
% (3) h:                    time step for time integration
% (4) tmax:                 simulation for [0,tmax]
% (5) varargin:             additional name-value pairs if actuation is
%                           present. See below for the parser function.
%
% OUTPUTS:
% (1) TI_NL_PROM:           struct containing the solutions and related
%                           information
%
% Additional note: not tested with actuation yet
%     
% Last modified: 01/04/2023, Mathieu Dubied, ETH Zurich

function TI_NL_FOM = solve_EoMs_FOM(FOMAssembly,tensors_hydro,h,tmax,varargin)

    % PARSE ADDITIONAL INPUTS (ACTUATION) _________________________________
    [ACTUATION,tensors_topMuscle_FOM,tensors_bottomMuscle_FOM] = parse_inputs(varargin{:});

    % SIMULATION PARAMETERS AND ICs _______________________________________
    nUncDOFs = size(FOMAssembly.Mesh.EBC.unconstrainedDOFs,2);
    eta0 = zeros(nUncDOFs,1);
    etad0 = zeros(nUncDOFs,1);
    etadd0 = zeros(nUncDOFs,1);

    % HYDRODYNAMIC AND ACTUATION (EXTERNAL FORCES) ________________________
    if ACTUATION
        fprintf('Solving for a structure with actuation...\n')
        fprintf('Funtion was not tested yet, there is no guarantee...\n')
        k = 1;
        B1TopMuscle = tensors_topMuscle_FOM.B1;
        B2TopMuscle = tensors_topMuscle_FOM.B2;
        B1BottomMuscle = tensors_bottomMuscle_FOM.B1;
        B2BottomMuscle = tensors_bottomMuscle_FOM.B2;
        % hydrodynamic and actuation forces
        T1 = FOMAssembly.constrain_vector(double(tensors_hydro.T1));
        Tu2 = FOMAssembly.constrain_matrix(double(tensors_hydro.Tu2));
        Tudot2 = FOMAssembly.constrain_matrix(double(tensors_hydro.Tudot2));
        Tuu3 = tensor(FOMAssembly.constrain_tensor(double(tensors_hydro.Tuu3)));
        Tuudot3 = tensor(FOMAssembly.constrain_tensor(double(tensors_hydro.Tuudot3)));
        Tudotudot3 = tensor(FOMAssembly.constrain_tensor(double(tensors_hydro.Tudotudot3)));

        F_ext = @(t,q,qd) (T1 + Tu2*q + Tudot2*q + double(ttv(ttv(Tuu3,q,3), q,2)) + ...
                            double(ttv(ttv(Tuudot3,qd,3),q,2))+ double(ttv(ttv(Tudotudot3,qd,3), qd,2))+ ...
                k/2*(1-(1+0.004*sin(t*2*pi/5)))*(B1TopMuscle+B2TopMuscle*eta) + ...
                k/2*(1-(1-0.004*sin(t*2*pi/5)))*(B1BottomMuscle+B2BottomMuscle*eta));

    else
        T1 = FOMAssembly.constrain_vector(double(tensors_hydro.T1));
        Tu2 = FOMAssembly.constrain_matrix(double(tensors_hydro.Tu2));
        Tudot2 = FOMAssembly.constrain_matrix(double(tensors_hydro.Tudot2));
        Tuu3 = tensor(FOMAssembly.constrain_tensor(double(tensors_hydro.Tuu3)));
        Tuudot3 = tensor(FOMAssembly.constrain_tensor(double(tensors_hydro.Tuudot3)));
        Tudotudot3 = tensor(FOMAssembly.constrain_tensor(double(tensors_hydro.Tudotudot3)));
        
        F_ext = @(t,q,qd) (T1 + Tu2*q + Tudot2*q + double(ttv(ttv(Tuu3,q,3), q,2)) + ...
                            double(ttv(ttv(Tuudot3,qd,3),q,2))+ double(ttv(ttv(Tudotudot3,qd,3), qd,2)));
    end
    
    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate object
    TI_NL_FOM = ImplicitNewmark('timestep',h,'alpha',0.005);
    
    % modal nonlinear Residual evaluation function handle
    Residual_NL_red = @(eta,etad,etadd,t)residual_nonlinear_hydro(eta,etad,etadd,t,FOMAssembly,F_ext);
    
    % time integration
    
    TI_NL_FOM.Integrate(eta0,etad0,etadd0,tmax,Residual_NL_red);
    TI_NL_FOM.Solution.u = zeros(FOMAssembly.Mesh.nDOFs,size(TI_NL_FOM.Solution.q,2));
    for t=1:size(TI_NL_FOM.Solution.q,2)
        TI_NL_FOM.Solution.u(:,t) = FOMAssembly.unconstrain_vector(TI_NL_FOM.Solution.q(:,t));
    end

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