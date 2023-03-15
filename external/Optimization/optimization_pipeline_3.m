% optimization_pipeline_3
%
% Synthax:
% [xiStar,LrEvo] = optimization_pipeline_3(MeshNominal,nodes,elements,U,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline 3 presented in
% the paper, based on Newton's method. 
%
% INPUTS: 
% (1) MeshNominal:  nominal mesh converted from Abaqus              
% (2) nodes:        nodes and their coordinates
% (3) elements:     elements and corresponding nodes
% (4) U:            shape variation basis
% (5) d:            forward swimming direction
% (6) h:            time step for time integration
% (7) tmax:         simulation for [0,tmax]
% (9)-(10)          constraints on xi of the form Axi<b
%
% possible additional name-value pair arguments
% (11) maxIteration:maximum number of iterations
% (12) convCrit:    convergence criterium. Norm between two successive
%                   optimal paramter vectors
% (13) FORMULATION: order of the Neumann approximation (N0/N1/N1t)
% (14) VOLUME:      integration over defected (1) or nominal volume (0)
% (15) USEJULIA:    use of JULIA (1) for the computation of internal forces
%                   tensors
%
% OUTPUTS:
% (1) xi_star:      optimal shape parameter(s) (scalar or vector)
% (2) LrEvo:        evolution of the cost function values
%     
%
% Last modified: 15/03/2023, Mathieu Dubied, ETH Zurich

function [xiStar,xiEvo,LrEvo] = optimization_pipeline_3(MeshNominal,nodes,elements,U,d,h,tmax,A,b,varargin)
    
    % parse input
    [maxIteration,convCrit,FORMULATION,VOLUME,USEJULIA,FOURTHORDER] = parse_inputs(varargin{:});
    
    % STEP 1: set xi_0 = 0 ________________________________________________
    fprintf('____________________\n')
    fprintf('STEP 1\n')

    xi_k = 0;
    xiEvo = xi_k;

    % STEP 2: mesh the structure and build a PROM _________________________
    fprintf('____________________\n')
    fprintf('STEP 2\n')

    [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM] = ...
        build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER);
        
    % STEP 3: solve EoMs to get nominal solution eta and dot{eta} _________
    fprintf('____________________\n')
    fprintf('STEP 3\n')
    
    TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);

    % STEP 4: solve sensitivity equation to get S _________________________
    fprintf('____________________\n')
    fprintf('STEP 4\n') 

    TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM, ...
        tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
        TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER);
    
    % STEP 5-13: OPTIMIZATION LOOP ________________________________________
    fprintf('____________________\n')
    fprintf('STEP 5-13\n') 

    eta = TI_NL_PROM.Solution.q;
    etad = TI_NL_PROM.Solution.qd;
    S = TI_sens.Solution.q;
    Sd = TI_sens.Solution.qd;
    eta_k = eta;
    etad_k = etad;

    N = size(eta,2);
    dr = reduced_constant_vector(d,V);
    Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xi_k,dr,A,b);
    LrEvo = Lr;

    for k = 1:maxIteration
        fprintf('Optimization loop iteration: k= %d\n',k-1)
        % step 7
        eta_k = eta_k + S*xi_k; 
        % step 8
        nablaLr = gradient_cost_function_w_constraints(dr,xi_k,eta_k,etad_k,S,Sd,A,b,tensors_hydro_PROM,FOURTHORDER);
        LrEvo = [LrEvo, reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta_k,etad_k,xi_k,dr,A,b)];
        % step 9 and 10
        xi_k = xi_k - 0.8*nablaLr;
        if k>40
            xi_k = xi_k - 0.3*nablaLr;
        end
        xiEvo = [xiEvo,xi_k];

        if norm(xiEvo(end)-xiEvo(end-1))<0.001
            fprintf('Convergence criterion of %.2g fulfilled\n',convCrit)
            break
        end
        if k == maxIteration
            fprintf('Maximum number of %d iterations reached\n',maxIteration)
        end
    end
    
    xiStar = xi_k;

end

% parse input
function [maxIteration,convCrit,FORMULATION,VOLUME,USEJULIA,FOURTHORDER] = parse_inputs(varargin)
defaultMaxIteration = 50;
defaultConvCrit = 0.001;
defaultFORMULATION = 'N1';
defaultVOLUME = 1;
defaultUSEJULIA = 0; 
defaultFOURTHORDER = 0;

p = inputParser;
addParameter(p,'maxIteration',defaultMaxIteration, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'convCrit',defaultConvCrit,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'FORMULATION',defaultFORMULATION,@(x)validateattributes(x, ...
                {'char'},{'nonempty'}))
addParameter(p,'VOLUME',defaultVOLUME,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'USEJULIA',defaultUSEJULIA,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'FOURTHORDER',defaultFOURTHORDER,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );

parse(p,varargin{:});

maxIteration = p.Results.maxIteration;
convCrit = p.Results.convCrit;
FORMULATION = p.Results.FORMULATION;
VOLUME = p.Results.VOLUME;
USEJULIA = p.Results.USEJULIA;
FOURTHORDER = p.Results.FOURTHORDER;

end